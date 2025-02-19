# TODO: add date correction option
# TODO: check/display stats file for errors (https://www.ncbi.nlm.nih.gov/genbank/validation/#BioSourceMissing)
import os

import click

from gb_submitter import utils


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
    Cmd = "table2asn "
    Cmd += f"-i {input} "
    Cmd += f"-o {output}.sqn "
    Cmd += f"-t {template} "
    Cmd += f"-f {features} "
    Cmd += f"-src-file {src_file} "
    Cmd += f"-w {comments} "
    Cmd += "-V v "  # Check for errors
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
