import os
import subprocess
import rich_click as click
from rich.console import Console

STYLES = {
    "info": "blue",
    "success": "green",
    "error": "bold red",
    "warning": "yellow",
    "debug": "cyan",
    "dim": "dim",
    "highlight": "bold magenta",
}


def resolve_tree_format(
    tree_path: str, specified_format: str | None, console: Console, debug: bool
) -> str:
    """
    Resolves the format of a phylogenetic tree file based on its extension or specified_format.
    Args:
        tree_path (str): Path to the tree file.
        specified_format (str | None): User-specified format ('newick' or 'nexus').
        console (Console): Rich console for output.
        debug (bool): If True, prints debug information.
    Returns:
        str: Resolved format ('newick' or 'nexus').
    Raises:
        click.Abort: If the format cannot be determined and no format is specified.
    """
    resolved_format = specified_format
    if not resolved_format:
        _, ext = os.path.splitext(tree_path)
        ext_lower = ext.lower()
        if ext_lower in [".nwk", ".newick"]:
            resolved_format = "newick"
        elif ext_lower == ".nexus":
            resolved_format = "nexus"
        else:
            if console:
                console.print(
                    f"[{STYLES['error']}]Error: Unknown tree format for file '{tree_path}'. Extension '{ext}' is not recognized.[/{STYLES['error']}]"
                )
                console.print(
                    f"[{STYLES['error']}]Please specify the format using --tree-format ('newick' or 'nexus').[/]"
                )
            raise click.Abort()

    if debug and console:
        console.print(
            f"[{STYLES['warning']}]Resolved tree format for '{tree_path}': {resolved_format}[/]"
        )
    return resolved_format


def run_subprocess_command(
    cmd: list[str],
    console: Console,
    debug: bool,
    success_message: str = "Successfully executed command.",
    error_message_prefix: str = "Error executing command",
) -> bool:
    """
    Runs a subprocess command and handles errors with rich output.
    Args:
        cmd (list[str]): Command to run as a list of strings.
        console (Console): Rich console for output.
        debug (bool): If True, prints debug information.
        success_message (str): Message to print on successful execution.
        error_message_prefix (str): Prefix for error messages.
    Returns:
        bool: True if the command was executed successfully, False otherwise.
    Raises:
        click.Abort: If the command fails or is not found.
    """
    if debug and console:
        console.print(f"[{STYLES['debug']}]Running command: {' '.join(cmd)}[/]")

    try:
        process_result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if debug and console:
            if process_result.stdout:
                console.print(
                    f"[{STYLES['dim']}]{cmd[0]} stdout:\\\\n{process_result.stdout}[/]"
                )
            if process_result.stderr:
                console.print(
                    f"[{STYLES['dim']}]{cmd[0]} stderr:\\\\n{process_result.stderr}[/]"
                )
        if success_message and console:
            console.print(f"[{STYLES['success']}]{success_message}[/]")
        return True
    except FileNotFoundError:
        if console:
            console.print(
                f"[{STYLES['error']}]{error_message_prefix}: {cmd[0]} command not found. Please ensure it is installed and in your PATH.[/]"
            )
        raise click.Abort()
    except subprocess.CalledProcessError as e:
        if console:
            console.print(f"[{STYLES['error']}]{error_message_prefix} {cmd[0]}: {e}[/]")
            if e.stdout:
                console.print(f"[{STYLES['dim']}]{cmd[0]} stdout:\\\\n{e.stdout}[/]")
            if e.stderr:
                console.print(f"[{STYLES['dim']}]{cmd[0]} stderr:\\\\n{e.stderr}[/]")
        raise click.Abort()


def sortFun(x: str) -> int:
    """
    Sort function to extract the numeric position from a mutation string.
    This function is used to sort mutation strings based on their numeric position,
    ignoring the nucleotide identities. It extracts the numeric part of the mutation
    string, which is expected to be in the format 'nuc_position', where 'nuc' is a
    nucleotide identity and 'position' is a numeric value.
    Args:
        x (str): Mutation string in the format 'nuc_position'.
    Returns:
        int: The numeric position extracted from the mutation string.
    """
    # sort based on nuc position, ignoring nuc identities
    return int(x[1 : (len(x) - 1)])
