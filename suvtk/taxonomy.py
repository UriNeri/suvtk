"""
taxonomy.py
===========

This script assigns taxonomy to sequences using MMseqs2 from the ICTV nr database.
It generates taxonomy files and integrates with other modules for further processing.

Functions
---------
taxonomy(fasta_file, database, output_path, seqid, threads)
    Main command to assign taxonomy to sequences.
"""

# TODO: save MIUVIG file with pred_genome_type and pred_genome_struc
# TODO: Change to genomad taxonomy
# TODO: mmseqs overwrite tmp file (mmseqs fails when command was previously aborted) ||| why not use --remove-tmp-files ? shouldn't unless you intend to try and reuse the tmp stuff from a failed run (--force-reuse)
# TODO: add memory size check and then manually set the  --split INT and --split-memory-limit BYTE so that at leas 200mb of the system availalbe memory are not used.

import os
import shutil

import click
import polars as pl

from suvtk import utils, virus_info


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
    help="Path to the suvtk database folder.",
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
    default=utils.get_available_cpus(),
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
    # TODO removing tmp before running mmseqs might be dangerous
    if os.path.exists("tmp"):
        shutil.rmtree("tmp")

    Cmd = "mmseqs easy-taxonomy "
    Cmd += f"{fasta_file} "  # input
    Cmd += os.path.join(database, "ictv_nr_db") + " "  # database
    Cmd += f"{taxresult_path} "  # output
    Cmd += "tmp "  # tmp
    Cmd += "-s 7.5 --blacklist '' --tax-lineage 1 "
    Cmd += f"--threads {threads}"
    utils.Exec(Cmd)

    shutil.rmtree("tmp")

    taxonomy = utils.safe_read_csv(f"{taxresult_path}_lca.tsv", separator="\t", has_header=False)
    taxonomy_columns = {
        "column_1": "query",
        "column_2": "taxid",
        "column_3": "rank",
        "column_4": "name",
        "column_5": "fragments",
        "column_6": "assigned",
        "column_7": "agreement",
        "column_8": "support",
        "column_9": "lineage",
    }
    taxonomy = taxonomy.rename(taxonomy_columns)

    tophit = utils.safe_read_csv(f"{taxresult_path}_tophit_aln", separator="\t", has_header=False)
    tophit_columns = {
        "column_1": "query",
        "column_2": "target",
        "column_3": "pident",
        "column_4": "len",
        "column_5": "mismatch",
        "column_6": "gapopen",
        "column_7": "qstart",
        "column_8": "qend",
        "column_9": "tstart",
        "column_10": "tend",
        "column_11": "evalue",
        "column_12": "bits",
    }
    tophit = tophit.rename(tophit_columns)

    # Select top hits using polars window function
    top_tophit = tophit.filter(
        pl.col("bits") == pl.col("bits").max().over("query")
    )

    merged = taxonomy.join(top_tophit, on="query", how="left")

    tax_names = []
    for row in merged.iter_rows(named=True):
        if row["rank"] == "no rank":
            click.echo(f"No taxonomy for {row['query']}")
            last_known = "unclassified viruses"
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
            ]
        )

    tax_df = pl.DataFrame(tax_names, columns=["contig", "taxonomy"])
    tax_df.write_csv(os.path.join(output_path, "taxonomy.tsv"), separator="\t")

    virus_info.run_segment_info(tax_df, database, output_path)


if __name__ == "__main__":
    taxonomy()
