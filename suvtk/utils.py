import os
import subprocess
import sys

import click
import pandas as pd


# Adapted from https://github.com/rcedgar/palm_annot/blob/77ac88ef7454dd3be9e5cbdb55792ce1ed7db95c/py/palm_annot.py#L121-L132
def Exec(CmdLine, fLog=None, capture=False):
    """
    Execute a command line in a shell, logging it to a file if specified.
    If capture is True, suppress printing output to the screen.

    :param CmdLine: The command line to execute
    :type CmdLine: str
    :param fLog: A file object to log the command and results, or None
    :type fLog: file object or None
    :param capture: Whether to capture output instead of printing
    :type capture: bool
    :return: The output of the command if captured, else None
    :rtype: str or None
    """

    def log_or_print(message, is_error=False):
        """Helper to log to file or print to screen, unless capturing."""
        if not capture:  # Only print if capture is False
            if fLog:
                fLog.write(message)
            else:
                output = sys.stderr if is_error else sys.stdout
                output.write(message)

    try:
        # Execute the command
        result = subprocess.run(
            CmdLine,
            shell=True,
            capture_output=capture,
            text=True,
            check=True,
        )

        # Log stdout
        if result.stdout:
            log_or_print(result.stdout)

        # Log stderr
        if result.stderr:
            log_or_print(result.stderr, is_error=True)

        return result.stdout if capture else None  # Return stdout only if capturing

    except subprocess.CalledProcessError as e:
        # Log error details
        if e.stderr:
            log_or_print(e.stderr, is_error=True)
        log_or_print(f"code {e.returncode}\n")
        log_or_print("\n")
        log_or_print(f"{CmdLine}\n")
        log_or_print("\n")
        log_or_print(f"Error code {e.returncode}\n", is_error=True)

        raise  # Re-raise the exception


def safe_read_csv(path, **kwargs):
    """
    Reads a CSV file using ASCII encoding. If a UnicodeDecodeError occurs,
    raises a ClickException showing the offending character.
    """
    try:
        return pd.read_csv(path, encoding="ascii", **kwargs)
    except UnicodeDecodeError as e:
        offending_bytes = e.object[e.start : e.end]
        # Try decoding using UTF-8 to show the offending character
        try:
            offending_char = offending_bytes.decode("utf-8")
        except Exception:
            offending_char = repr(offending_bytes)
        raise click.ClickException(
            f"Only ASCII characters are allowed in file '{path}'. "
            f"Offending character: {offending_char}. Error: {str(e)}"
        )


# Copied from https://github.com/EricDeveaud/genomad/blob/030ab6434654435ce75243347c97be6f40ea175b/genomad/cli.py#L250-L257
def get_available_cpus():
    try:
        # Try to get the number of cores available to this process
        CPUS = len(os.sched_getaffinity(0))
    except AttributeError:
        # Windows / MacOS probably don't have this functionality
        CPUS = os.cpu_count()
    return CPUS
