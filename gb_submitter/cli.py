#TODO: format help to show everything
import click
from gb_submitter.features import features
from gb_submitter.gbk2tbl import gbk2tbl
from gb_submitter.taxonomy import taxonomy

from gettext import gettext as _

class FullHelpGroup(click.Group):
    def format_commands(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        """Extra format methods for multi methods that adds all the commands
        after the options.
        """
        commands = []
        for subcommand in self.list_commands(ctx):
            cmd = self.get_command(ctx, subcommand)
            # What is this, the tool lied about a command.  Ignore it
            if cmd is None:
                continue
            if cmd.hidden:
                continue

            commands.append((subcommand, cmd))

        # allow for 3 times the default spacing
        if len(commands):
            limit = formatter.width - 6 - max(len(cmd[0]) for cmd in commands)

            rows = []
            for subcommand, cmd in commands:
                help = cmd.help if cmd.help is not None else ""
                rows.append((subcommand, help))

            if rows:
                with formatter.section(_("Commands")):
                    formatter.write_dl(rows)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], show_default=True)

@click.group(context_settings=CONTEXT_SETTINGS, cls=FullHelpGroup)
@click.version_option()
def cli():
    """Tool to submit viral sequences to Genbank."""

cli.add_command(gbk2tbl)
cli.add_command(features)
cli.add_command(taxonomy)

if __name__ == "__main__":
    cli()