"""Converts a nexus file to a newick file."""

import dendropy
from ete4 import Tree
from pathlib import Path
import re


def _remove_quotes_from_file(file_path: Path | str) -> None:
    """
    Remove all single and double quotes from the contents of a file.
    This function reads the entire contents of the file located at the given
    path, strips out every instance of single quotes (') and double quotes (")
    using a regular expression, and writes the cleaned text back to the same
    file.
    Args:
        file_path (Path or str): Path to the target file, provided either
            as a pathlib.Path object or a string representing the filesystem path.
    Raises:
        FileNotFoundError: If the specified file does not exist.
        OSError: If an error occurs while opening, reading, or writing the file.
    Returns:
        None
    """

    path = Path(file_path)  # Ensure it's a Path object

    with path.open("r") as f:
        content = f.read()

    # Regex matching ' or "
    pattern = re.compile(r"['\"]")
    modified_content = pattern.sub("", content)

    with path.open("w") as f:
        f.write(modified_content)


def convert_nexus_to_newick(
    input_file: Path | str,
    output_file: Path | str,
    input_format: str = "nexus",
    reformat_tree: bool = False,
):
    """
    Convert a phylogenetic tree file from Nexus to Newick format.
    This function reads a tree in the specified input_format (default 'nexus')
    using DendroPy, writes it in Newick format to output_file, and optionally
    reformats the Newick string using ETE4 for cleaner output. Finally, it
    removes any quotes from the taxa labels in the resulting file.
    Args:
        input_file (Path or str): Path to the input tree file.
        output_file (Path or str): Path where the Newick-formatted tree will be written.
        input_format (str, optional): Format of the input tree (e.g., 'nexus', 'newick'). Defaults to 'nexus'.
        reformat_tree (bool, optional): If True, reparse and rewrite the Newick file using ETE4
            to produce a cleaner string. Defaults to False.
    Raises:
        OSError: If the input file cannot be read or the output file cannot be written.
        ValueError: If the input_format is not supported by DendroPy.
    Notes:
        - Taxon labels are case sensitive and underscores are preserved.
        - Internal node taxa are suppressed by default during conversion.
        - Any double quotes around taxon labels are removed in the final output.
    """
    input_file_str = str(input_file)
    output_file_path = Path(output_file)  # Use Path object for operations
    output_file_str = str(output_file_path)

    tree = dendropy.Tree.get_from_path(
        input_file_str,
        input_format,
        suppress_leaf_node_taxa=False,
        suppress_internal_node_taxa=True,
        case_sensitive_taxon_labels=True,
        preserve_underscores=True,
    )
    tree.write_to_path(output_file_str, schema="newick")

    if reformat_tree:
        # ETE4 can sometimes produce a cleaner Newick string
        tree_string_from_file = ""
        with open(output_file_str, "r") as handle:
            for line in handle:
                l = line.strip("\\n")
                if "(" in l:
                    tree_string_start = l.index("(")
                    tree_string_from_file = l[tree_string_start:]
                    break  # Assuming the tree is on one line or the first line with '('

        if tree_string_from_file:
            t = Tree(tree_string_from_file)
            t.write(
                parser=0, outfile=output_file_str
            )  # format=0 is a flexible Newick format

    _remove_quotes_from_file(output_file_path)
