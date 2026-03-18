"""Source https://gist.github.com/huddlej/5d7bd023d3807c698bd18c706974f2db"""

import json
import pandas as pd
import Bio.Phylo
from augur.utils import annotate_parents_for_tree
import click

from .utils import print_error, print_warning, print_success


def json_to_tree(json_dict, root=True, parent_cumulative_branch_length=None):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    Assigns links back to parent nodes for the root of the tree.

    Test opening a JSON from augur export v1.

    >>> import json
    >>> json_fh = open("tests/data/json_tree_to_nexus/flu_h3n2_ha_3y_tree.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> tree.name
    'NODE_0002020'
    >>> len(tree.clades)
    2
    >>> tree.clades[0].name
    'NODE_0001489'
    >>> hasattr(tree, "attr")
    True
    >>> "dTiter" in tree.attr
    True
    >>> tree.clades[0].parent.name
    'NODE_0002020'
    >>> tree.clades[0].branch_length > 0
    True

    Test opening a JSON from augur export v2.

    >>> json_fh = open("tests/data/zika.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> hasattr(tree, "name")
    True
    >>> len(tree.clades) > 0
    True
    >>> tree.clades[0].branch_length > 0
    True

    Branch lengths should be the length of the branch to each node and not the
    length from the root. The cumulative branch length from the root gets its
    own attribute.

    >>> tip = [tip for tip in tree.find_clades(terminal=True) if tip.name == "USA/2016/FLWB042"][0]
    >>> round(tip.cumulative_branch_length, 6)
    0.004747
    >>> round(tip.branch_length, 6)
    0.000186

    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    node.name = json_dict.get("name")
    if node.name is None:  # Fallback for v1 or if "name" is not present
        node.name = json_dict.get("strain")

    # Assign all non-children attributes from the JSON node to the Clade object.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Handle specific attributes like 'num_date', 'div' (cumulative_branch_length),
    # and 'translations' based on JSON version (v1 uses 'attr', v2 uses 'node_attrs').
    if hasattr(node, "attr"):  # v1 style
        node.numdate = node.attr.get("num_date")
        node.cumulative_branch_length = node.attr.get("div")
        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):  # v2 style
        node.cumulative_branch_length = node.node_attrs.get("div")
        # If 'num_date' or 'translations' can also be in v2 'node_attrs', handle them here.
        # For example:
        # node.numdate = node.node_attrs.get("num_date", node.numdate) # If numdate might be elsewhere too
        # if "translations" in node.node_attrs:
        #    node.translations = node.node_attrs["translations"]

    node.branch_length = 0.0
    if (
        parent_cumulative_branch_length is not None
        and hasattr(node, "cumulative_branch_length")
        and node.cumulative_branch_length is not None
    ):
        node.branch_length = (
            node.cumulative_branch_length - parent_cumulative_branch_length
        )

    # Ensure branch_length is non-negative, as small floating point inaccuracies can occur.
    if node.branch_length < 0:
        node.branch_length = 0.0

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [
            json_to_tree(
                child,
                root=False,
                parent_cumulative_branch_length=node.cumulative_branch_length,
            )
            for child in json_dict["children"]
        ]

    if root:
        node = annotate_parents_for_tree(node)

    return node


def process_auspice_json(
    tree_json_path: str,
    output_metadata_path: str | None,
    output_tree_path: str | None,
    include_internal_nodes: bool,
    attributes: list[str] | None,
):
    """
    Converts an Auspice JSON tree to other formats (Newick, metadata table).
    """
    # Load tree from JSON.
    try:
        with open(tree_json_path, "r", encoding="utf-8") as fh:
            tree_json_data = json.load(fh)
    except FileNotFoundError:
        print_error(f"Tree JSON file not found at '{tree_json_path}'")
        raise click.Abort()
    except json.JSONDecodeError:
        print_error(f"Could not decode JSON from '{tree_json_path}'")
        raise click.Abort()

    tree = json_to_tree(tree_json_data)

    # Output the tree in Newick format, if requested.
    if output_tree_path:
        try:
            Bio.Phylo.write(
                tree,
                output_tree_path,
                "newick",
            )
            print_success(f"Tree successfully written to '{output_tree_path}'")
        except Exception as e:
            print_error(f"Error writing Newick tree to '{output_tree_path}': {e}")
            raise click.Abort()

    if output_metadata_path:
        records = []
        attributes_to_export = attributes if attributes else []

        if not attributes_to_export:  # If attributes list is empty or None, auto-detect
            attrs_set = set()
            if hasattr(tree.root, "attr"):  # v1 style
                attrs_set.update(tree.root.attr.keys())
            if hasattr(tree.root, "node_attrs"):  # v2 style
                attrs_set.update(tree.root.node_attrs.keys())
            if hasattr(tree.root, "branch_attrs"):  # v2 style, branch_attrs is optional
                attrs_set.update(getattr(tree.root, "branch_attrs", {}).keys())

            if not attrs_set:
                print_warning(
                    "Could not auto-detect any attributes from the tree root. "
                    "The metadata output might be sparse or only contain 'name'. "
                    "Consider using --attributes to specify columns."
                )
            attributes_to_export = sorted(list(attrs_set))

        for (
            node_obj
        ) in (
            tree.find_clades()
        ):  # Renamed 'node' to 'node_obj' to avoid conflict with 'node' module
            if node_obj.is_terminal() or include_internal_nodes:
                record = {"name": node_obj.name}
                for attr_name in attributes_to_export:
                    value = None
                    # Check v2 style locations first
                    if (
                        hasattr(node_obj, "node_attrs")
                        and attr_name in node_obj.node_attrs
                    ):
                        value = node_obj.node_attrs[attr_name]
                    elif (
                        hasattr(node_obj, "branch_attrs")
                        and attr_name in node_obj.branch_attrs
                    ):
                        value = node_obj.branch_attrs[attr_name]
                    # Check v1 style location
                    elif hasattr(node_obj, "attr") and attr_name in node_obj.attr:
                        value = node_obj.attr[attr_name]
                    # Fallback: check if it's a direct attribute of the node object
                    elif hasattr(node_obj, attr_name):
                        potential_value = getattr(node_obj, attr_name)
                        if not callable(potential_value) and not isinstance(
                            potential_value,
                            (Bio.Phylo.BaseTree.Clade, Bio.Phylo.BaseTree.TreeElement),
                        ):
                            if (
                                isinstance(potential_value, list)
                                and attr_name == "clades"
                            ):
                                pass
                            else:
                                value = potential_value

                    if isinstance(value, dict) and "value" in value:
                        value = value["value"]

                    if isinstance(value, list):
                        value = ";".join(map(str, value))
                    record[attr_name] = value
                records.append(record)

        if records:
            df = pd.DataFrame(records)
            final_columns = ["name"] + [
                attr for attr in attributes_to_export if attr in df.columns
            ]
            # Add any other columns that might have been added but were not in attributes_to_export.
            for col in df.columns:
                if col not in final_columns:
                    final_columns.append(col)
            df = df[final_columns]
        else:
            df = pd.DataFrame(columns=["name"] + attributes_to_export)

        try:
            df.to_csv(
                output_metadata_path,
                sep="\t",
                header=True,
                index=False,
            )
            print_success(f"Metadata successfully written to '{output_metadata_path}'")
        except Exception as e:
            print_error(f"Error writing metadata to '{output_metadata_path}': {e}")
            raise click.Abort()

    if not output_tree_path and not output_metadata_path:
        # This case should be handled by the CLI command to require at least one output.
        # However, a warning here is fine if called directly.
        print_warning(
            "No output requested. Use --output-tree and/or --output-metadata."
        )
