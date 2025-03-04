# TODO: save MIUVIG file with pred_genome_type and pred_genome_struc
# TODO: Change to genomad taxonomy
import importlib.resources
import os
import shutil

import click
import pandas as pd
import taxopy

from gb_submitter import utils


def load_segment_db():
    """
    Load the segmented viruses database.

    Returns
    -------
    pandas.DataFrame
        Data frame with the segmented viruses database.
    """
    with (
        importlib.resources.files("gb_submitter.data")
        .joinpath("segmented_viruses.tsv")
        .open("r") as file
    ):
        db = pd.read_csv(file, sep="\t", header=0)
        return db


def load_genome_type_db():
    """
    Load the genome structure database.

    Returns
    -------
    pandas.DataFrame
        Data frame with the genome structure database.
    """
    with (
        importlib.resources.files("gb_submitter.data")
        .joinpath("genome_types.tsv")
        .open("r") as file
    ):
        db = pd.read_csv(file, sep="\t", header=0)
        return db


@click.command(short_help="Assign virus taxonomy to sequences.")
@click.option(
    "-i",
    "--input",
    "fasta_file",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(exists=False),
    help="Output directory",
)
@click.option(
    "-d",
    "--database",
    "database",
    required=True,
    type=click.Path(exists=True),
    help="ICTV MMseqs database path",
)
@click.option(
    "-s",
    "--identity",
    "seqid",
    required=False,
    default=0.7,
    type=float,
    help="Minimum sequence identity for hits to be considered",
)
@click.option(
    "-t",
    "--threads",
    "threads",
    required=False,
    default=4,
    type=int,
    help="Number of threads to use",
)
def taxonomy(fasta_file, database, output_path, seqid, threads):
    """
    This command uses MMseqs2 to assign taxonomy to sequences using protein sequences from ICTV taxa in the nr database.
    """
    if os.path.exists(output_path):
        click.echo(
            f"Warning: Output directory '{output_path}' already exists and will be overwritten."
        )

    os.makedirs(output_path, exist_ok=True)

    taxresult_path = os.path.join(output_path, "taxresults")

    # TODO Add RAM restrictions?
    # TODO Add error handling
    Cmd = "mmseqs easy-taxonomy "
    Cmd += f"{fasta_file} "  # input
    Cmd += f"{database} "  # database
    Cmd += f"{taxresult_path} "  # output
    Cmd += "tmp "  # tmp
    Cmd += "-s 7.5 --blacklist '' --tax-lineage 1 "
    Cmd += f"--threads {threads}"
    utils.Exec(Cmd)

    shutil.rmtree("tmp")

    taxonomy = pd.read_csv(f"{taxresult_path}_lca.tsv", sep="\t", header=None)
    taxonomy.rename(
        {
            0: "query",
            1: "taxid",
            2: "rank",
            3: "name",
            4: "fragments",
            5: "assigned",
            6: "agreement",
            7: "support",
            8: "lineage",
        },
        axis=1,
        inplace=True,
    )

    tophit = pd.read_csv(f"{taxresult_path}_tophit_aln", sep="\t", header=None)
    tophit.rename(
        {
            0: "query",
            1: "target",
            2: "pident",
            3: "len",
            4: "mismatch",
            5: "gapopen",
            6: "qstart",
            7: "qend",
            8: "tstart",
            9: "tend",
            10: "evalue",
            11: "bits",
        },
        axis=1,
        inplace=True,
    )

    # Select top hits
    highest_bits_idx = tophit.groupby("query")["bits"].idxmax()
    top_tophit = tophit.loc[highest_bits_idx]

    merged = pd.merge(taxonomy, top_tophit, on="query", how="left")

    tax_names = []
    for index, row in merged.iterrows():
        if row["rank"] == "no rank":
            click.echo(f"No taxonomy for {row['query']}")
            last_known = "unclassified viruses"
            row["taxid"] = 12429
        elif row["rank"] == "species":  # Fix issue when species contains sp.
            last_known = row["lineage"].split(";")[-2].replace("g_", "")
            last_known += " sp."
        elif (
            row["rank"] == "genus" and row["pident"] < seqid  # TODO check best cutoff
        ):  # if genus rank and sequence identity is lower than 70% (seqid) get family assignment
            last_known = row["lineage"].split(";")[-2].replace("f_", "")
            last_known += " sp."
        else:
            last_known = row["name"].strip()
            last_known += " sp."
        tax_names.append(
            [
                row["query"],
                last_known,
                row["taxid"],
            ]  # taxid is not equal to last_known (eg. species taxid -> last_known is max genus level)
        )

    tax_df = pd.DataFrame(tax_names, columns=["contig", "taxonomy", "taxid"])
    tax_df.to_csv(os.path.join(output_path, "taxonomy.tsv"), sep="\t", index=False)

    # Load database file with segmented viruses
    segment_db = load_segment_db()
    genome_type_db = load_genome_type_db()

    # Load taxonomy database
    taxdb = taxopy.TaxDb(
        nodes_dmp="/lustre1/scratch/337/vsc33750/ictv_db/ictv_taxdump/nodes.dmp",
        names_dmp="/lustre1/scratch/337/vsc33750/ictv_db/ictv_taxdump/names.dmp",
    )  # TODO: Set database path

    results = []  # For segmented virus details from segment_db
    gt_results = []  # For genome type records from genome_type_db
    segmented = False  # Flag to trigger the extra click.echo messages

    for index, row in tax_df.iterrows():
        # Get the lineage for this taxid
        lineage = taxopy.Taxon(row["taxid"], taxdb).name_lineage
        segment_record = None
        genome_type_record = None

        # Loop once over the lineage to get both segment and genome type records.
        for taxa in lineage:
            if segment_record is None and taxa in segment_db["taxon"].values:
                segment_record = (
                    segment_db.loc[segment_db["taxon"] == taxa].iloc[0].to_dict()
                )
                # Record the segment details (for segmented virus info output)
                record = {"contig": row["contig"], **segment_record}
                results.append(record)

                # If a high fraction of viruses are segmented, print extra information
                # (Using >= 25 as a threshold here)
                if float(segment_record["segmented_fraction"]) >= 25:
                    segmented = True
                    click.echo(
                        f"\n{row['contig']} is part of the {segment_record['taxon']} {segment_record['rank']}, "
                        f"{float(segment_record['segmented_fraction']):.2f}% of these are segmented viruses."
                    )
                    if segment_record["min_segment"] != segment_record["max_segment"]:
                        click.echo(
                            f"Most segmented viruses of the {segment_record['rank']} {segment_record['taxon']} have "
                            f"{segment_record['majority_segment']} segments, but it can vary between "
                            f"{segment_record['min_segment']} and {segment_record['max_segment']} depending on the species."
                        )
                    else:
                        click.echo(
                            f"The segmented viruses of the {segment_record['rank']} {segment_record['taxon']} have "
                            f"{segment_record['majority_segment']} segments."
                        )

            if genome_type_record is None and taxa in genome_type_db["taxon"].values:
                genome_type_record = (
                    genome_type_db.loc[genome_type_db["taxon"] == taxa]
                    .iloc[0]
                    .to_dict()
                )

            # If we've found both, we can stop checking further taxa in the lineage.
            if segment_record is not None and genome_type_record is not None:
                break

        # Compute the predicted genome structure based on the segmented_fraction from segment_record.
        # If no segment record was found, we assume non-segmented.
        if segment_record is not None:
            try:
                seg_frac = float(segment_record["segmented_fraction"])
            except (ValueError, TypeError):
                seg_frac = 0.0
            if seg_frac == 100:
                pred_struc = "segmented"
            elif seg_frac > 0:
                pred_struc = "undetermined"
            else:
                pred_struc = "non-segmented"
        else:
            pred_struc = "non-segmented"

        # Add the predicted genome structure to the genome type record.
        # If no genome type record was found, create a minimal record.
        if genome_type_record is not None:
            genome_type_record.pop("taxon")
            genome_type_record["contig"] = row["contig"]
            genome_type_record["pred_genome_struc"] = pred_struc
            gt_results.append(genome_type_record)
        else:
            gt_results.append(
                {
                    "contig": row["contig"],
                    # "taxon": "unclassified viruses",
                    "pred_genome_type": "uncharacterized",
                    "pred_genome_struc": "undetermined",
                }
            )

    # Write the genome type DataFrame (with the added predicted structure) to file.
    genome_type_df = pd.DataFrame(gt_results)

    ordered_columns = ["contig", "pred_genome_type", "pred_genome_struc"]
    genome_type_df = genome_type_df[ordered_columns]

    genome_type_df.to_csv(
        os.path.join(output_path, "miuvig_taxonomy.tsv"),
        sep="\t",
        index=False,
    )

    if segmented:
        click.echo(
            "\nYou might want to look into your data to see if you can identify the missing segments."
        )

    # Write the segmented virus info to file if any records exist.
    if results:
        segmented_df = pd.DataFrame(results)
        segmented_df = segmented_df.sort_values(
            by="segmented_fraction", ascending=False
        )
        segmented_df.to_csv(
            os.path.join(output_path, "segmented_viruses_info.tsv"),
            sep="\t",
            index=False,
        )


if __name__ == "__main__":
    taxonomy()
