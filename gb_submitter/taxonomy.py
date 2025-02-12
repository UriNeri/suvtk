import os
import warnings

import click
import pandas as pd

from gb_submitter import utils


@click.command(help="Assign virus taxonomy to sequences.")
@click.option(
    "-i",
    "--input-file",
    "fasta_file",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file",
)
@click.option(
    "-o",
    "--output-path",
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
    if os.path.exists(output_path):
        warnings.warn(
            f"Warning: Output directory '{output_path}' already exists and may be overwritten."
        )

    os.makedirs(output_path, exist_ok=True)

    # Add RAM restrictions?
    # Add error handling
    Cmd = "mmseqs easy-taxonomy "
    Cmd += f"{fasta_file} "
    Cmd += f"{database} "
    Cmd += f"{output_path}/taxresults "
    Cmd += "tmp "
    Cmd += "-s 7.5 --blacklist '' --tax-lineage 1 "
    Cmd += f"--threads {threads}"
    utils.Exec(Cmd)

    taxonomy = pd.read_csv(f"{output_path}/taxresults_lca.tsv", sep="\t", header=None)
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

    tophit = pd.read_csv(f"{output_path}/taxresults_tophit_aln", sep="\t", header=None)
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
            print(f"No taxonomy for {row['query']}")
            last_known = "unclassified viruses"
            row["taxid"] = 12429
        elif row["rank"] == "species":  # Fix issue when species contains sp.
            last_known = row["lineage"].split(";")[-2].replace("g_", "")
            last_known += " sp."
        else:
            if (
                row["rank"] == "genus"
                and row["pident"] < seqid  # TODO check best cutoff
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
    tax_df.to_csv(f"{output_path}/taxonomy.tsv", sep="\t", index=False)


if __name__ == "__main__":
    taxonomy()
