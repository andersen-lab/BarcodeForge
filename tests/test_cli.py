import pytest
import json
import pandas as pd
from click.testing import CliRunner
from barcodeforge.cli import cli
from barcodeforge import __version__
from pathlib import Path
import shutil
from unittest.mock import call, ANY, MagicMock
from rich.console import Console
from barcodeforge.utils import STYLES


# Helper function to create dummy files
def create_dummy_file(path, content=""):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(content)


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def temp_files(tmp_path):
    ref_genome = tmp_path / "reference.fasta"
    alignment = tmp_path / "alignment.fasta"
    tree = tmp_path / "tree.nwk"
    lineages = tmp_path / "lineages.tsv"

    create_dummy_file(ref_genome, ">ref_genome_id\nAGCTAGCTAGCTAGCT")
    create_dummy_file(
        alignment,
        ">seq1\nAGCTAGCTAGCTAGCT\n>seq2\nAGCTAGCTAGCTCGCT\n>seq3\nAGCTAGCTAGCTAGGT",
    )
    create_dummy_file(tree, "((seq1:0.1,seq2:0.1):0.05,seq3:0.15);")
    create_dummy_file(lineages, "clade\tsequences\nlineageA\tseq1,seq2\nlineageB\tseq3")

    return {
        "ref_genome": str(ref_genome),
        "alignment": str(alignment),
        "tree": str(tree),
        "lineages": str(lineages),
        "tmp_path": tmp_path,
    }


def test_cli_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert result.output.strip() == f"cli version {__version__}"


def test_cli_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output


def test_barcode_command_help(runner):
    result = runner.invoke(cli, ["barcode", "--help"])
    assert result.exit_code == 0
    assert "REFERENCE_GENOME" in result.output


@pytest.mark.skipif(
    not shutil.which("faToVcf")
    or not shutil.which("usher")
    or not shutil.which("matUtils"),
    reason="External dependencies (faToVcf, usher, matUtils) not found. Skipping test_barcode_command_execution.",
)
def test_barcode_command_execution(runner, temp_files, mocker):
    # This test remains skipped as per subtask instructions.
    pass


def test_barcode_command_default_options(runner, temp_files, mocker):
    mock_run_subp = mocker.patch(
        "barcodeforge.cli.run_subprocess_command", return_value=True
    )
    mock_resolve_format = mocker.patch(
        "barcodeforge.cli.resolve_tree_format", return_value="newick"
    )
    mock_convert_tree = mocker.patch("barcodeforge.cli.convert_nexus_to_newick")
    mock_process_reroot = mocker.patch("barcodeforge.cli.process_and_reroot_lineages")
    mock_create_barcodes = mocker.patch(
        "barcodeforge.cli.create_barcodes_from_lineage_paths"
    )
    mock_create_plot = mocker.patch("barcodeforge.cli.create_barcode_plot")
    mock_cli_console = mocker.patch("barcodeforge.cli.console", MagicMock(spec=Console))

    args = [
        "barcode",
        temp_files["ref_genome"],
        temp_files["alignment"],
        temp_files["tree"],
        temp_files["lineages"],
    ]
    result = runner.invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, f"CLI failed: {result.output}"

    intermediate_dir = "barcodeforge_workdir"
    aligned_vcf_fn = f"{intermediate_dir}/aligned.vcf"
    converted_tree_fn = f"{intermediate_dir}/converted_tree.nwk"
    tree_pb_fn = f"{intermediate_dir}/tree.pb"
    annotated_tree_pb_fn = f"{intermediate_dir}/annotated_tree.pb"
    matutils_C_output_fn = f"{intermediate_dir}/lineagePaths.txt"
    matutils_S_output_fn = f"{intermediate_dir}/samplePaths.txt"
    auspice_json_fn = "auspice_tree.json"
    additional_muts_processed_fn = f"{intermediate_dir}/additional_mutations.tsv"
    rerooted_lineage_paths_fn = f"{intermediate_dir}/rerooted_lineage_paths.txt"
    final_barcodes_csv_fn = "barcode.csv"
    final_barcode_plot_fn = "barcode_plot.pdf"

    mock_resolve_format.assert_called_once_with(
        temp_files["tree"], None, mock_cli_console, False
    )
    mock_convert_tree.assert_called_once_with(
        input_file=temp_files["tree"],
        output_file=converted_tree_fn,
        input_format="newick",
    )

    expected_usher_cmd = [
        "usher",
        "-t",
        converted_tree_fn,
        "-v",
        aligned_vcf_fn,
        "-o",
        tree_pb_fn,
        "-T",
        "8",
    ]
    expected_subprocess_calls = [
        call(
            ["faToVcf", temp_files["alignment"], aligned_vcf_fn],
            mock_cli_console,  # console passed
            False,  # debug status passed
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            expected_usher_cmd,
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            [
                "matUtils",
                "annotate",
                "--set-overlap",
                "0.0",  # Default overlap from cli.py
                "-i",
                tree_pb_fn,
                "-c",
                temp_files["lineages"],
                "-o",
                annotated_tree_pb_fn,
                "-T",
                "8",
            ],
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            [
                "matUtils",
                "extract",
                "-i",
                annotated_tree_pb_fn,
                "-C",
                matutils_C_output_fn,
                "-S",
                matutils_S_output_fn,
                "-j",
                auspice_json_fn,
                "-T",
                "8",
            ],
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
    ]
    mock_run_subp.assert_has_calls(expected_subprocess_calls, any_order=False)

    mock_process_reroot.assert_called_once_with(
        sample_muts_path=matutils_S_output_fn,
        reference_fasta_path=temp_files["ref_genome"],
        sequences_fasta_path=temp_files["alignment"],
        input_lineage_paths_path=matutils_C_output_fn,
        output_additional_muts_path=additional_muts_processed_fn,
        output_rerooted_lineage_paths_path=rerooted_lineage_paths_fn,
    )
    mock_create_barcodes.assert_called_once_with(
        input_file_path=rerooted_lineage_paths_fn,
        output_file_path=final_barcodes_csv_fn,
        prefix="",  # Default prefix ""
    )
    mock_create_plot.assert_called_once_with(
        input_file_path=final_barcodes_csv_fn,
        debug=False,  # Default debug status
        output_file_path=final_barcode_plot_fn,
    )


def test_barcode_command_custom_options(runner, temp_files, mocker):
    mock_run_subp = mocker.patch(
        "barcodeforge.cli.run_subprocess_command", return_value=True
    )
    mock_resolve_format = mocker.patch(
        "barcodeforge.cli.resolve_tree_format", return_value="newick"
    )
    mock_convert_tree = mocker.patch("barcodeforge.cli.convert_nexus_to_newick")
    mock_process_reroot = mocker.patch("barcodeforge.cli.process_and_reroot_lineages")
    mock_create_barcodes = mocker.patch(
        "barcodeforge.cli.create_barcodes_from_lineage_paths"
    )
    mock_create_plot = mocker.patch("barcodeforge.cli.create_barcode_plot")
    mock_cli_console = mocker.patch("barcodeforge.cli.console", MagicMock(spec=Console))

    prefix = "MYPREFIX"
    custom_usher_args = "-U -l"
    custom_threads = "4"
    custom_matutils_overlap = "0.5"

    args = [
        "barcode",
        temp_files["ref_genome"],
        temp_files["alignment"],
        temp_files["tree"],
        temp_files["lineages"],
        "--prefix",
        prefix,
        "--usher-args",
        custom_usher_args,
        "--threads",
        custom_threads,
        "--matutils-overlap",
        str(custom_matutils_overlap),
    ]
    result = runner.invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, f"CLI failed: {result.output}"

    intermediate_dir = "barcodeforge_workdir"
    aligned_vcf_fn = f"{intermediate_dir}/aligned.vcf"
    converted_tree_fn = f"{intermediate_dir}/converted_tree.nwk"
    tree_pb_fn = f"{intermediate_dir}/tree.pb"
    annotated_tree_pb_fn = f"{intermediate_dir}/annotated_tree.pb"
    matutils_C_output_fn = f"{intermediate_dir}/lineagePaths.txt"
    matutils_S_output_fn = f"{intermediate_dir}/samplePaths.txt"
    auspice_json_fn = "auspice_tree.json"

    additional_muts_processed_fn = f"{intermediate_dir}/additional_mutations.tsv"
    rerooted_lineage_paths_fn = f"{intermediate_dir}/rerooted_lineage_paths.txt"

    final_barcodes_csv_fn = f"{prefix}-barcode.csv"
    final_barcode_plot_fn = f"{prefix}-barcode_plot.pdf"

    mock_resolve_format.assert_called_once_with(
        temp_files["tree"], None, mock_cli_console, False
    )
    mock_convert_tree.assert_called_once_with(
        input_file=temp_files["tree"],
        output_file=converted_tree_fn,
        input_format="newick",  # Because mock_resolve_format returns "newick"
    )

    expected_usher_cmd = (
        ["usher"]
        + custom_usher_args.split()
        + [
            "-t",
            converted_tree_fn,
            "-v",
            aligned_vcf_fn,
            "-o",
            tree_pb_fn,
            "-T",
            custom_threads,
        ]
    )
    expected_subprocess_calls = [
        call(
            ["faToVcf", temp_files["alignment"], aligned_vcf_fn],
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            expected_usher_cmd,
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            [
                "matUtils",
                "annotate",
                "--set-overlap",
                str(custom_matutils_overlap),
                "-i",
                tree_pb_fn,
                "-c",
                temp_files["lineages"],
                "-o",
                annotated_tree_pb_fn,
                "-T",
                custom_threads,
            ],
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
        call(
            [
                "matUtils",
                "extract",
                "-i",
                annotated_tree_pb_fn,
                "-C",
                matutils_C_output_fn,
                "-S",
                matutils_S_output_fn,
                "-j",
                auspice_json_fn,
                "-T",
                custom_threads,
            ],
            mock_cli_console,
            False,
            success_message=ANY,
            error_message_prefix=ANY,
        ),
    ]
    mock_run_subp.assert_has_calls(expected_subprocess_calls, any_order=False)

    mock_process_reroot.assert_called_once_with(
        sample_muts_path=matutils_S_output_fn,
        reference_fasta_path=temp_files["ref_genome"],
        sequences_fasta_path=temp_files["alignment"],
        input_lineage_paths_path=matutils_C_output_fn,
        output_additional_muts_path=additional_muts_processed_fn,
        output_rerooted_lineage_paths_path=rerooted_lineage_paths_fn,
    )
    mock_create_barcodes.assert_called_once_with(
        input_file_path=rerooted_lineage_paths_fn,
        output_file_path=final_barcodes_csv_fn,
        prefix=prefix,
    )
    mock_create_plot.assert_called_once_with(
        input_file_path=final_barcodes_csv_fn,
        debug=False,  # Default debug status
        output_file_path=final_barcode_plot_fn,
    )


def test_barcode_command_nexus_tree(runner, temp_files, mocker):
    mock_run_subp = mocker.patch(
        "barcodeforge.cli.run_subprocess_command", return_value=True
    )
    mock_resolve_format = mocker.patch(
        "barcodeforge.cli.resolve_tree_format", return_value="nexus"
    )
    mock_convert_tree = mocker.patch("barcodeforge.cli.convert_nexus_to_newick")
    mocker.patch("barcodeforge.cli.process_and_reroot_lineages")
    mocker.patch("barcodeforge.cli.create_barcodes_from_lineage_paths")
    mocker.patch("barcodeforge.cli.create_barcode_plot")
    mock_cli_console = mocker.patch("barcodeforge.cli.console", MagicMock(spec=Console))

    args = [
        "barcode",
        temp_files["ref_genome"],
        temp_files["alignment"],
        temp_files["tree"],
        temp_files["lineages"],
        "--tree_format",
        "nexus",
    ]
    result = runner.invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, f"CLI failed: {result.output}"

    intermediate_dir = "barcodeforge_workdir"
    converted_tree_fn = f"{intermediate_dir}/converted_tree.nwk"
    aligned_vcf_fn = f"{intermediate_dir}/aligned.vcf"
    tree_pb_fn = f"{intermediate_dir}/tree.pb"

    mock_resolve_format.assert_called_once_with(
        temp_files["tree"], "nexus", mock_cli_console, False
    )
    mock_convert_tree.assert_called_once_with(
        input_file=temp_files["tree"],
        output_file=converted_tree_fn,
        input_format="nexus",
    )

    expected_usher_cmd = [
        "usher",
        "-t",
        converted_tree_fn,
        "-v",
        aligned_vcf_fn,
        "-o",
        tree_pb_fn,
        "-T",
        "8",
    ]
    mock_run_subp.assert_any_call(
        expected_usher_cmd,
        mock_cli_console,
        False,
        success_message=ANY,
        error_message_prefix=ANY,
    )


def test_barcode_command_newick_tree_reformat(runner, temp_files, mocker):
    mocker.patch("barcodeforge.cli.run_subprocess_command", return_value=True)
    mock_resolve_format = mocker.patch(
        "barcodeforge.cli.resolve_tree_format", return_value="newick"
    )
    mock_convert_tree = mocker.patch("barcodeforge.cli.convert_nexus_to_newick")
    mocker.patch("barcodeforge.cli.process_and_reroot_lineages")
    mocker.patch("barcodeforge.cli.create_barcodes_from_lineage_paths")
    mocker.patch("barcodeforge.cli.create_barcode_plot")
    mock_cli_console = mocker.patch("barcodeforge.cli.console", MagicMock(spec=Console))

    args = [
        "barcode",
        temp_files["ref_genome"],
        temp_files["alignment"],
        temp_files["tree"],
        temp_files["lineages"],
        "--tree_format",
        "newick",
    ]
    result = runner.invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, f"CLI failed: {result.output}"

    converted_tree_fn = "barcodeforge_workdir/converted_tree.nwk"
    mock_resolve_format.assert_called_once_with(
        temp_files["tree"], "newick", mock_cli_console, False
    )
    mock_convert_tree.assert_called_once_with(
        input_file=temp_files["tree"],
        output_file=converted_tree_fn,
        input_format="newick",
    )


def test_barcode_command_debug_flag(runner, temp_files, mocker):
    mock_run_subp = mocker.patch(
        "barcodeforge.cli.run_subprocess_command", return_value=True
    )
    mock_resolve_format = mocker.patch(
        "barcodeforge.cli.resolve_tree_format", return_value="newick"
    )
    mock_convert_tree = mocker.patch("barcodeforge.cli.convert_nexus_to_newick")
    mock_process_reroot = mocker.patch("barcodeforge.cli.process_and_reroot_lineages")
    mock_create_barcodes = mocker.patch(
        "barcodeforge.cli.create_barcodes_from_lineage_paths"
    )
    mock_create_plot = mocker.patch("barcodeforge.cli.create_barcode_plot")
    mock_cli_console = mocker.patch("barcodeforge.cli.console", MagicMock(spec=Console))

    args = [
        "--debug",  # Main CLI debug flag
        "barcode",
        temp_files["ref_genome"],
        temp_files["alignment"],
        temp_files["tree"],
        temp_files["lineages"],
    ]

    result = runner.invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, f"CLI command failed: {result.output}"

    intermediate_dir = "barcodeforge_workdir"
    matutils_C_output_fn = f"{intermediate_dir}/lineagePaths.txt"
    matutils_S_output_fn = f"{intermediate_dir}/samplePaths.txt"
    additional_muts_processed_fn = f"{intermediate_dir}/additional_mutations.tsv"
    rerooted_lineage_paths_fn = f"{intermediate_dir}/rerooted_lineage_paths.txt"
    final_barcodes_csv_fn = "barcode.csv"
    final_barcode_plot_fn = "barcode_plot.pdf"

    mock_cli_console.print.assert_any_call(
        f"[{STYLES['debug']}]Debug mode is ON[/{STYLES['debug']}]"
    )
    mock_resolve_format.assert_called_once_with(
        temp_files["tree"], None, mock_cli_console, True
    )
    mock_convert_tree.assert_called_once_with(
        input_file=temp_files["tree"],
        output_file="barcodeforge_workdir/converted_tree.nwk",  # Unprefixed
        input_format="newick",
    )

    mock_process_reroot.assert_called_once_with(
        sample_muts_path=matutils_S_output_fn,
        reference_fasta_path=temp_files["ref_genome"],
        sequences_fasta_path=temp_files["alignment"],
        input_lineage_paths_path=matutils_C_output_fn,
        output_additional_muts_path=additional_muts_processed_fn,
        output_rerooted_lineage_paths_path=rerooted_lineage_paths_fn,
    )
    mock_create_barcodes.assert_called_once_with(
        input_file_path=rerooted_lineage_paths_fn,
        output_file_path=final_barcodes_csv_fn,
        prefix="",  # Default prefix ""
    )
    mock_create_plot.assert_called_once_with(
        input_file_path=final_barcodes_csv_fn,
        debug=True,  # Debug status passed
        output_file_path=final_barcode_plot_fn,
    )

    for call_obj in mock_run_subp.call_args_list:
        assert call_obj.args[1] is mock_cli_console
        assert call_obj.args[2] is True


def test_barcode_command_missing_file(runner, temp_files):
    args = [
        "barcode",
        temp_files["ref_genome"],
        "non_existent_alignment.fasta",
        temp_files["tree"],
        temp_files["lineages"],
    ]
    result = runner.invoke(cli, args)
    assert result.exit_code != 0
    assert "Invalid value for 'ALIGNMENT'" in result.output


def test_extract_auspice_data_command(runner, tmp_path):
    json_path = tmp_path / "tree.json"
    json_data = {
        "meta": {},
        "tree": {
            "name": "root",
            "node_attrs": {"div": 0.0, "country": {"value": "USA"}},
            "children": [
                {"name": "A", "node_attrs": {"div": 0.1, "country": {"value": "USA"}}},
                {"name": "B", "node_attrs": {"div": 0.2, "country": {"value": "CAN"}}},
            ],
        },
    }
    json_path.write_text(json.dumps(json_data))

    meta_out = tmp_path / "meta.tsv"
    tree_out = tmp_path / "tree.nwk"

    result = runner.invoke(
        cli,
        [
            "extract-auspice-data",
            str(json_path),
            "--output_metadata_path",
            str(meta_out),
            "--output_tree_path",
            str(tree_out),
            "--include_internal_nodes",
            "--attributes",
            "country",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert meta_out.exists()
    assert tree_out.exists()
    df = pd.read_csv(meta_out, sep="\t")
    assert set(df["name"]) == {"root", "A", "B"}
    assert list(df.columns) == ["name", "country"]
