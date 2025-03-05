# TODO: add date correction option
# TODO: check/display stats file for errors (https://www.ncbi.nlm.nih.gov/genbank/validation/#BioSourceMissing)
# TODO: add missing required miuvig params from src file? -> comments.py?:  collection_date, env_broad_scale, env_local_scale, env_medium, geo_loc_name, ivestigation_type, lat_lon, project_name, seq_meth
# TODO: add option to make genbank file
import os

import click
import pandas as pd

from gb_submitter import utils


def process_comments(src_file, comments_file):
    """
    Processes comments by updating the comments file based on the source file.

    This function reads a source file and a comments file into pandas DataFrames,
    groups the source DataFrame by the "Isolate" column, and processes duplicates.
    For each group of isolates with more than one entry, it updates the comments
    file with the count of isolates, the majority predicted genome type, and sets
    the predicted genome structure to "segmented".

    Args:
        src_file (str): The file path to the source file in tab-separated format.
        comments_file (str): The file path to the comments file in tab-separated format.

    Returns:
        None
    """
    # Read the src_file and comments file into DataFrames
    src_df = pd.read_csv(src_file, sep="\t")
    comments_df = pd.read_csv(comments_file, sep="\t")

    # Group by Isolate and process duplicates
    isolate_groups = src_df.groupby("Isolate")

    for isolate, group in isolate_groups:
        if len(group) > 1:
            # Get the count of Isolate value
            isolate_count = len(group)

            # Get the majority pred_genome_type
            majority_genome_type = group["pred_genome_type"].mode()[0]

            # Update the comments file
            for seqid in group["Sequence_ID"]:
                comments_df.loc[
                    comments_df["Sequence_ID"] == seqid, "number_contig"
                ] = isolate_count
                comments_df.loc[
                    comments_df["Sequence_ID"] == seqid, "pred_genome_type"
                ] = majority_genome_type
                comments_df.loc[
                    comments_df["Sequence_ID"] == seqid, "pred_genome_struc"
                ] = "segmented"

    # Save the updated comments file
    comments_df.to_csv(comments_file, sep="\t", index=False)


@click.command(short_help="Generate .sqn submission for Genbank.")
@click.option(
    "-i",
    "--input",
    "input",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file",
)
@click.option(
    "-o",
    "--output",
    "output",
    required=True,
    type=click.Path(exists=False),
    help="Output prefix",
)
@click.option(
    "-s",
    "--src-file",
    "src_file",
    required=True,
    type=click.Path(exists=True),
    help="File with Source modifiers (.src).",
)
@click.option(
    "-f",
    "--features",
    "features",
    required=True,
    type=click.Path(exists=True),
    help="Feature table file (.tbl).",
)
@click.option(
    "-t",
    "--template",
    "template",
    required=True,
    type=click.Path(exists=True),
    help="Template file with author information (.sbt). See https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
)
@click.option(
    "-c",
    "--comments",
    "comments",
    required=True,
    type=click.Path(exists=True),
    help="Structured comment file (.cmt) with MIUVIG information.",
)
def table2asn(input, output, src_file, features, template, comments):
    """This command generates a .sqn file that you can send to gb-sub@ncbi.nlm.nih.gov"""

    # Process the comments file based on the src_file
    process_comments(src_file, comments)

    Cmd = "table2asn "
    Cmd += f"-i {input} "
    Cmd += f"-o {output}.sqn "
    Cmd += f"-t {template} "
    Cmd += f"-f {features} "
    Cmd += f"-src-file {src_file} "
    Cmd += f"-w {comments} "
    Cmd += "-V vb "  # Check for errors
    Cmd += "-a s"  # allow multifasta file

    utils.Exec(Cmd)

    tag = 0
    error_file = f"{output}.val"
    with open(error_file, "r") as f:
        for line in f.readlines():
            if line.startswith("Warning") or line.startswith("Info"):
                # click.echo(f"{line}")
                continue
            else:
                click.echo(f"UNEXPECTED ERROR -- {line}")
                tag = 1

    if tag == 0:
        print("No major errors reported for Genbank submission")
