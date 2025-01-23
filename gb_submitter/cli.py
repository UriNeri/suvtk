import click
from gb_submitter.features import features
from gb_submitter.gbk2tbl import gbk2tbl

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def cli():
    """Tool to submit viral sequences to Genbank."""

cli.add_command(gbk2tbl)
cli.add_command(features)

if __name__ == "__main__":
    cli()