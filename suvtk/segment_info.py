#TODO: Rename to genome_info?
import importlib.resources
import os
from pathlib import Path

import click
import pandas as pd
import taxopy

def load_segment_db():
    """
    Load the segmented viruses database.

    Returns
    -------
    pandas.DataFrame
        Data frame with the segmented viruses database.
    """
    with (
        importlib.resources.files("suvtk.data")
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
        importlib.resources.files("suvtk.data")
        .joinpath("genome_types.tsv")
        .open("r") as file
    ):
        db = pd.read_csv(file, sep="\t", header=0)
        return db
    
@click.command(short_help="Get information on possible segmented viruses based on their taxonomy.")
@click.option(
    "--taxonomy",
    required=True,
    type=click.Path(exists=True),
    help="Taxonomy file.",
)
@click.option(
    "-d",
    "--database",
    "database",
    required=True,
    type=click.Path(exists=True),
    help="The suvtk database path.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(exists=False),
    help="Output directory",
)

def segment_info(taxonomy, database, output_path):
    #TODO: generate docstring for long help
    if isinstance(taxonomy, (str, Path)):
        tax_df = pd.read_csv(taxonomy, sep="\t")
    elif isinstance(taxonomy, pd.DataFrame):
        tax_df = taxonomy
    else:
        raise TypeError("taxonomy must be either a file path (str/Path) or a pandas DataFrame")
    
    #if os.path.exists(output_path):
    #    click.echo(
    #        f"Warning: Output directory '{output_path}' already exists and will be overwritten."
    #    )

    os.makedirs(output_path, exist_ok=True)

    # Load database file with segmented viruses
    segment_db = load_segment_db()
    genome_type_db = load_genome_type_db()

    # Load taxonomy database
    taxdb = taxopy.TaxDb(
        nodes_dmp=os.path.join(database, "nodes.dmp"),
        names_dmp=os.path.join(database, "names.dmp"),
    )  # TODO: Set better database path?

    results = []  # For segmented virus details from segment_db
    gt_results = []  # For genome type records from genome_type_db
    segmented = False  # Flag to trigger the extra click.echo messages

    for index, row in tax_df.iterrows():
        if row["taxonomy"] == "unclassified viruses":
            gt_results.append(
                {
                    "contig": row["contig"],
                    "pred_genome_type": "uncharacterized",
                    "pred_genome_struc": "undetermined",
                }
            )
            continue

        # Get the lineage for this taxid
        tax = row["taxonomy"].removesuffix(' sp.')
        try:
            taxid = taxopy.taxid_from_name(tax, taxdb)
            lineage = taxopy.Taxon(taxid[0], taxdb).name_lineage
        except IndexError:
            click.echo(f"Warning: '{tax}' is not part of the official ICTV taxonomy. The info on genome type and structure can therefore not be accessed for {row["contig"]}.")
            gt_results.append(
                {
                    "contig": row["contig"],
                    "pred_genome_type": "uncharacterized",
                    "pred_genome_struc": "undetermined",
                }
            )
            continue
        
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
    segment_info()