import sys
import subprocess

def Exec(CmdLine, fLog=None, capture=False):
    """
    Execute a command line in a shell, logging it to a file if specified,
    or printing output to the screen if no log file is given.

    :param CmdLine: The command line to execute
    :type CmdLine: str
    :param fLog: A file object to log the command and results, or None
    :type fLog: file object or None
    :return: The output of the command
    :rtype: str
    """
    def log_or_print(message, is_error=False):
        """Helper to log to file or print to screen."""
        if fLog:
            fLog.write(message)
        else:
            output = sys.stderr if is_error else sys.stdout
            output.write(message)

    try:
        # Execute the command and capture output
        result = subprocess.run(
            CmdLine,
            shell=True,
            capture_output=capture,
            text=True,
            check=True
        )
        # Print or log stdout
        if result.stdout:
            log_or_print(result.stdout)
        # Print or log stderr
        if result.stderr:
            log_or_print(result.stderr, is_error=True)

        return result.stdout  # Return the command's stdout
    except subprocess.CalledProcessError as e:
        # Print or log error details
        if e.stderr:
            log_or_print(e.stderr, is_error=True)
        log_or_print(f"code {e.returncode}\n")
        log_or_print("\n")
        log_or_print(f"{CmdLine}\n")
        log_or_print("\n")
        log_or_print(f"Error code {e.returncode}\n", is_error=True)

        raise  # Re-raise the exception to notify the caller