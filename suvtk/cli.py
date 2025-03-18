from gettext import gettext as _

import click

from suvtk import (
    co_occurrence,
    comments,
    features,
    gbk2tbl,
    table2asn,
    taxonomy,
    virus_info,
)


class FullHelpGroup(click.Group):
    """
    Custom Click Group to display commands in the order they were added.

    Methods
    -------
    list_commands(ctx: click.Context)
        Return commands in the order they were added.
    format_commands(ctx: click.Context, formatter: click.HelpFormatter)
        Formats and displays commands in the correct order.
    """

    def list_commands(self, ctx: click.Context):
        """
        Return commands in the order they were added.

        Parameters
        ----------
        ctx : click.Context
            The Click context.

        Returns
        -------
        list
            List of command names in the order they were added.
        """
        return list(self.commands.keys())

    def format_commands(self, ctx: click.Context, formatter: click.HelpFormatter):
        """
        Formats and displays commands in the correct order.

        Parameters
        ----------
        ctx : click.Context
            The Click context.
        formatter : click.HelpFormatter
            The Click help formatter.
        """
        commands = [
            (name, self.get_command(ctx, name))
            for name in self.list_commands(ctx)
            if self.get_command(ctx, name) and not self.get_command(ctx, name).hidden
        ]

        if commands:
            rows = [(name, cmd.short_help or "") for name, cmd in commands]

            with formatter.section(_("Commands")):
                formatter.write_dl(rows)


CONTEXT_SETTINGS = dict(
    help_option_names=["-h", "--help"], show_default=True, max_content_width=120
)


@click.group(context_settings=CONTEXT_SETTINGS, cls=FullHelpGroup)
@click.version_option()
def cli():
    """
    Tool to submit viral sequences to Genbank.

    This CLI tool provides various commands to process and submit viral
    sequences to Genbank, including taxonomy, features, virus information,
    and more.
    """


cli.add_command(taxonomy)
cli.add_command(features)
cli.add_command(virus_info)
cli.add_command(co_occurrence)
cli.add_command(gbk2tbl)
cli.add_command(comments)
cli.add_command(table2asn)


if __name__ == "__main__":
    cli()
