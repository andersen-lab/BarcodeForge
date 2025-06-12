import os
import rich_click as click
from rich.console import Console
from .format_tree import convert_nexus_to_newick
from .utils import (
    resolve_tree_format,
    run_subprocess_command,
    STYLES,
)  # Import new utils
from .ref_muts import process_and_reroot_lineages
from .generate_barcodes import create_barcodes_from_lineage_paths
from .plot_barcode import create_barcode_plot
from .auspice_tree_to_table import process_auspice_json
from . import __version__

console = Console()


@click.group()
@click.version_option(__version__, message="%(prog)s version %(version)s")
@click.option("--debug/--no-debug", default=False, help="Enable debug mode.")
@click.pass_context
def cli(ctx, debug):
    """BarcodeForge CLI"""
    ctx.ensure_object(dict)
    ctx.obj["DEBUG"] = debug
    if debug:
        console.print(f"[{STYLES['debug']}]Debug mode is ON[/{STYLES['debug']}]")


@cli.command()
@click.argument("reference_genome", type=click.Path(exists=True, readable=True))
@click.argument("alignment", type=click.Path(exists=True, readable=True))
@click.argument("tree", type=click.Path(exists=True, readable=True))
@click.argument("lineages", type=click.Path(exists=True, readable=True))
@click.option(
    "--tree_format",
    type=click.Choice(["newick", "nexus"], case_sensitive=False),
    help="Specify the format of the tree file (newick or nexus)",
)
@click.option(
    "--usher-args",
    type=str,
    default="",
    help="Additional arguments to pass to usher (e.g., '-U -l'). Quote multiple arguments.",
)
@click.option(
    "--threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of CPUs/threads to use.",
)
@click.option(
    "--matutils-overlap",
    type=float,
    metavar="FLOAT",
    default="0",
    show_default=True,
    help="Value for --set-overlap in matUtils annotate.",
)
@click.option(
    "--prefix",
    type=str,
    default="",
    show_default=True,
    help="Prefix to add to lineage names in the barcode file.",
)
@click.pass_context
def barcode(
    ctx,
    reference_genome,
    alignment,
    tree,
    lineages,
    tree_format,
    usher_args,
    threads,
    matutils_overlap,
    prefix,
):
    """Process barcode data, including VCF generation, tree formatting, USHER placement, matUtils annotation, and matUtils extraction."""
    is_debug = ctx.obj.get("DEBUG", False)  # More robust way to get debug status

    console.print(
        f"[bold {STYLES['info']}]Reference Genome:[/bold {STYLES['info']}]\t{reference_genome}"
    )
    console.print(
        f"[bold {STYLES['info']}]Alignment:[/bold {STYLES['info']}]\t\t{alignment}"
    )
    console.print(f"[bold {STYLES['info']}]Tree:[/bold {STYLES['info']}]\t\t\t{tree}")
    console.print(
        f"[bold {STYLES['info']}]Lineages:[/bold {STYLES['info']}]\t\t{lineages}"
    )

    # Use utility function to resolve tree format
    resolved_tree_format = resolve_tree_format(tree, tree_format, console, is_debug)

    # Directory for intermediate files
    intermediate_dir = "barcodeforge_workdir"
    os.makedirs(intermediate_dir, exist_ok=True)

    # Run faToVcf command using utility function
    fatovcf_output_vcf = os.path.join(intermediate_dir, "aligned.vcf")
    fatovcf_cmd = ["faToVcf", alignment, fatovcf_output_vcf]
    run_subprocess_command(
        fatovcf_cmd,
        console,
        is_debug,
        success_message=f"Successfully created VCF file: {fatovcf_output_vcf}",
        error_message_prefix="Error running faToVcf",
    )

    output_converted_tree_path = os.path.join(intermediate_dir, "converted_tree.nwk")
    usher_output_pb = os.path.join(intermediate_dir, "tree.pb")

    if resolved_tree_format == "nexus":
        if is_debug:
            console.print(
                f"[{STYLES['info']}]Converting Nexus tree ({tree}) to Newick format at {output_converted_tree_path}...[/{STYLES['info']}]"
            )
        convert_nexus_to_newick(
            input_file=tree,
            output_file=output_converted_tree_path,
            input_format="nexus",
        )
        console.print(
            f"[{STYLES['success']}]Converted tree saved to {output_converted_tree_path}[/{STYLES['success']}]"
        )
    elif resolved_tree_format == "newick":
        if is_debug:
            console.print(
                f"[{STYLES['info']}]Processing Newick tree ({tree}) to {output_converted_tree_path}...[/{STYLES['info']}]"
            )
        convert_nexus_to_newick(
            input_file=tree,
            output_file=output_converted_tree_path,
            input_format="newick",
        )
        console.print(
            f"[{STYLES['success']}]Processed tree saved to {output_converted_tree_path} (if conversion/reformatting occurred)[/{STYLES['success']}]"
        )
    else:
        raise ValueError(
            f"Unsupported tree format: {resolved_tree_format}. Expected 'newick' or 'nexus'."
        )

    # Run usher command
    usher_cmd = ["usher"]
    if usher_args:
        usher_cmd.extend(usher_args.split())  # Split space-separated string of args
    usher_cmd.extend(
        [
            "-t",
            output_converted_tree_path,
            "-v",
            fatovcf_output_vcf,
            "-o",
            usher_output_pb,
            "-T",
            str(threads),
        ]
    )

    run_subprocess_command(
        usher_cmd,
        console,
        is_debug,
        success_message=f"Successfully ran USHER. Output protobuf: {usher_output_pb}",
        error_message_prefix="Error running USHER",
    )

    # Run matUtils annotate command
    matutils_annotate_cmd = ["matUtils", "annotate"]

    # Specific args for annotate
    matutils_annotate_cmd.extend(
        [
            "--set-overlap",
            str(matutils_overlap),
            "-i",
            usher_output_pb,  # Input protobuf from usher
            "-c",
            lineages,  # Input clades/lineages file
            "-o",
            os.path.join(
                intermediate_dir, "annotated_tree.pb"
            ),  # Output annotated protobuf
            "-T",
            str(threads),  # Number of threads
        ]
    )

    run_subprocess_command(
        matutils_annotate_cmd,
        console,
        is_debug,
        success_message=f"Successfully ran matUtils annotate. Output: annotated_tree.pb",
        error_message_prefix="Error running matUtils annotate",
    )

    # Run matUtils extract command
    matutils_extract_cmd = ["matUtils", "extract"]

    matutils_extract_cmd.extend(
        [
            "-i",
            os.path.join(
                intermediate_dir, "annotated_tree.pb"
            ),  # Input annotated protobuf
            "-C",
            os.path.join(intermediate_dir, "lineagePaths.txt"),
            "-S",
            os.path.join(intermediate_dir, "samplePaths.txt"),
            "-j",
            "auspice_tree.json",
            "-T",
            str(threads),  # Number of threads
        ]
    )

    run_subprocess_command(
        matutils_extract_cmd,
        console,
        is_debug,
        success_message=f"Successfully ran matUtils extract. Outputs: lineagePaths.txt, samplePaths.txt, auspice_tree.json",
        error_message_prefix="Error running matUtils extract",
    )

    process_and_reroot_lineages(
        sample_muts_path=os.path.join(intermediate_dir, "samplePaths.txt"),
        reference_fasta_path=reference_genome,
        sequences_fasta_path=alignment,
        input_lineage_paths_path=os.path.join(intermediate_dir, "lineagePaths.txt"),
        output_additional_muts_path=os.path.join(
            intermediate_dir, "additional_mutations.tsv"
        ),
        output_rerooted_lineage_paths_path=os.path.join(
            intermediate_dir, "rerooted_lineage_paths.txt"
        ),
    )

    # Determine base name for output files
    clean_prefix = prefix.strip()
    if clean_prefix:
        console.print(
            f"[{STYLES['info']}]Using prefix '{clean_prefix}' for lineage names in barcodes...[/{STYLES['info']}]"
        )
        base_name = f"{clean_prefix}-barcode"
    else:
        console.print(
            f"[{STYLES['info']}]No prefix provided for lineage names in barcodes.[/{STYLES['info']}]"
        )
        base_name = "barcode"

    csv_path = f"{base_name}.csv"
    barcode_render = f"{base_name}_plot.pdf"

    create_barcodes_from_lineage_paths(
        input_file_path=os.path.join(intermediate_dir, "rerooted_lineage_paths.txt"),
        output_file_path=csv_path,
        prefix=clean_prefix,
    )
    create_barcode_plot(
        input_file_path=csv_path,
        output_file_path=barcode_render,
    )

    console.print(
        f"[{STYLES['success']}]Generated barcodes are saved to '{csv_path}' and plot saved to '{barcode_render}'[/{STYLES['success']}]"
    )


@cli.command()
@click.argument("auspice_json_path", type=click.Path(exists=True, readable=True))
@click.option(
    "--output_metadata_path",
    type=click.Path(),
    default="metadata.csv",
    show_default=True,
    help="Path to save the metadata table (CSV format).",
)
@click.option(
    "--output_tree_path",
    type=click.Path(),
    default="tree.nwk",
    show_default=True,
    help="Path to save the tree in Newick format.",
)
@click.option(
    "--include_internal_nodes",
    is_flag=True,
    default=False,
    show_default=True,
    help="Include internal nodes in the output tree.",
)
@click.option(
    "--attributes",
    type=str,
    multiple=True,
    help="Attributes to include in the metadata table (e.g., 'country', 'date').",
)
@click.pass_context
def extract_auspice_data(
    ctx,
    auspice_json_path,
    output_metadata_path,
    output_tree_path,
    include_internal_nodes,
    attributes,
):
    """
    Extract metadata and tree from an Auspice JSON file.
    Inspired by Dr. John Huddleston's Gist on processing Auspice JSON files.
    Source: https://gist.github.com/huddlej/5d7bd023d3807c698bd18c706974f2db
    """
    is_debug = ctx.obj.get("DEBUG", False)  # More robust way to get debug status
    console.print(
        f"[bold {STYLES['info']}]Processing Auspice JSON file:[/bold {STYLES['info']}]\t{auspice_json_path}"
    )
    process_auspice_json(
        tree_json_path=auspice_json_path,
        output_metadata_path=output_metadata_path,
        output_tree_path=output_tree_path,
        include_internal_nodes=include_internal_nodes,
        attributes=list(attributes) if attributes else None,
        console=console,
    )


def main():
    """Entry point for the CLI application."""
    cli()


if __name__ == "__main__":
    main()
