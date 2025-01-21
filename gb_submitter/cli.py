import click
from gb_submitter.gbk2tbl import gbk2tbl

@click.group()
@click.version_option()
def cli():
    """Tool to submit viral sequences to Genbank."""

cli.add_command(gbk2tbl)

if __name__ == "__main__":
    cli()