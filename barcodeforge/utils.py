import os
import re
import subprocess
import rich_click as click
from rich.console import Console
from Bio import SeqIO

console = Console()

STYLES = {
    "info": "blue",
    "success": "green",
    "error": "bold red",
    "warning": "yellow",
    "debug": "cyan",
}


def print_error(msg: str) -> None:
    console.print(f"[ERROR] {msg}", style=STYLES["error"], markup=False)


def print_warning(msg: str) -> None:
    console.print(f"[WARNING] {msg}", style=STYLES["warning"], markup=False)


def print_success(msg: str) -> None:
    console.print(f"[SUCCESS] {msg}", style=STYLES["success"], markup=False)


def print_info(msg: str, bold: bool = False) -> None:
    style = f"bold {STYLES['info']}" if bold else STYLES["info"]
    console.print(f"[INFO] {msg}", style=style, markup=False)


def print_debug(msg: str) -> None:
    console.print(f"[DEBUG] {msg}", style=STYLES["debug"], markup=False)


def resolve_tree_format(
    tree_path: str, specified_format: str | None, debug: bool
) -> str:
    """
    Resolves the format of a phylogenetic tree file based on its extension or specified_format.
    Args:
        tree_path (str): Path to the tree file.
        specified_format (str | None): User-specified format ('newick' or 'nexus').
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
            print_error(
                f"Unknown tree format for file '{tree_path}'. Extension '{ext}' is not recognized."
            )
            print_error(
                "Please specify the format using --tree-format ('newick' or 'nexus')."
            )
            raise click.Abort()

    if debug:
        print_debug(f"Resolved tree format for '{tree_path}': {resolved_format}")
    return resolved_format


def ensure_reference_is_first_in_alignment(
    reference_genome_path: str,
    alignment_path: str,
    output_alignment_path: str,
    debug: bool,
) -> str:
    """
    Ensures the reference genome is the first sequence in the alignment.

    faToVcf treats the first sequence in the alignment as the reference. If the
    reference genome is not already the first sequence, a corrected alignment is
    written to ``output_alignment_path`` with the reference genome prepended (and
    any existing record sharing the reference id removed to avoid duplicates), and
    that path is returned. If the reference is already first, the original
    alignment path is returned unchanged.

    Args:
        reference_genome_path (str): Path to the reference genome FASTA file.
        alignment_path (str): Path to the alignment FASTA file.
        output_alignment_path (str): Path to write the corrected alignment to if
            the reference genome needs to be added.
        debug (bool): If True, prints debug information.
    Returns:
        str: Path to the alignment to use downstream (the original path if no
            change was needed, otherwise ``output_alignment_path``).
    Raises:
        click.Abort: If the alignment is empty, or the reference genome length
            does not match the aligned sequence length.
    """
    try:
        ref = SeqIO.read(reference_genome_path, "fasta")
    except ValueError as e:
        # SeqIO.read raises ValueError if the file has zero or more than one record.
        print_error(
            f"Reference genome '{reference_genome_path}' must contain exactly one "
            f"sequence ({e})."
        )
        raise click.Abort()
    records = list(SeqIO.parse(alignment_path, "fasta"))

    if not records:
        print_error(f"Alignment file '{alignment_path}' contains no sequences.")
        raise click.Abort()

    if records[0].id == ref.id:
        if debug:
            print_debug(
                f"Reference genome '{ref.id}' is already the first sequence in the alignment."
            )
        return alignment_path

    print_warning(
        f"The first sequence in the alignment ('{records[0].id}') is not the "
        f"reference genome ('{ref.id}'). Prepending the reference genome to the alignment."
    )

    alignment_width = len(records[0].seq)
    if len(ref.seq) != alignment_width:
        print_error(
            f"Cannot add the reference genome to the alignment: the reference length "
            f"({len(ref.seq)}) does not match the aligned sequence length "
            f"({alignment_width}). faToVcf requires all sequences to be the same length. "
            f"Please provide a reference that is aligned to the same coordinates as the alignment."
        )
        raise click.Abort()

    # Drop any existing record that shares the reference id to avoid duplicate
    # sample names, then place the reference genome first.
    corrected = [ref] + [rec for rec in records if rec.id != ref.id]
    SeqIO.write(corrected, output_alignment_path, "fasta")
    print_success(
        f"Wrote corrected alignment with reference genome '{ref.id}' as the first "
        f"sequence to {output_alignment_path}"
    )
    return output_alignment_path


def run_subprocess_command(
    cmd: list[str],
    debug: bool,
    success_message: str = "Successfully executed command.",
    error_message_prefix: str = "Error executing command",
) -> bool:
    """
    Runs a subprocess command and handles errors with rich output.
    Args:
        cmd (list[str]): Command to run as a list of strings.
        debug (bool): If True, prints debug information.
        success_message (str): Message to print on successful execution.
        error_message_prefix (str): Prefix for error messages.
    Returns:
        bool: True if the command was executed successfully.
    Raises:
        click.Abort: If the command fails or is not found.
    """
    if debug:
        print_debug(f"Running command: {' '.join(cmd)}")

    try:
        process_result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if debug:
            if process_result.stdout:
                print_debug(f"{cmd[0]} stdout:\n{process_result.stdout}")
            if process_result.stderr:
                print_debug(f"{cmd[0]} stderr:\n{process_result.stderr}")
        if success_message:
            print_success(success_message)
        return True
    except FileNotFoundError:
        print_error(
            f"{error_message_prefix}: {cmd[0]} command not found. Please ensure it is installed and in your PATH."
        )
        raise click.Abort()
    except PermissionError:
        print_error(
            f"{error_message_prefix}: {cmd[0]} is not executable. Please check its file permissions."
        )
        raise click.Abort()
    except subprocess.CalledProcessError as e:
        print_error(f"{error_message_prefix} {cmd[0]}: {e}")
        if debug:
            if e.stdout:
                print_debug(f"{cmd[0]} stdout:\n{e.stdout}")
            if e.stderr:
                print_debug(f"{cmd[0]} stderr:\n{e.stderr}")
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
    match = re.search(r"\d+", x)
    if match is None:
        raise ValueError(f"Cannot extract a numeric position from mutation '{x}'")
    return int(match.group())
