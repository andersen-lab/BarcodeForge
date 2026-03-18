import pytest
import click
import subprocess
from unittest.mock import MagicMock, patch, call
from rich.console import Console
import barcodeforge.utils
from barcodeforge.utils import (
    sortFun,
    resolve_tree_format,
    run_subprocess_command,
    STYLES,
)


def test_sortFun():
    assert sortFun("A123B") == 123
    assert sortFun("C456D") == 456
    assert sortFun("G789E") == 789


def test_resolve_tree_format_specified_newick():
    assert resolve_tree_format("any.tree", "newick", False) == "newick"


def test_resolve_tree_format_specified_nexus():
    assert resolve_tree_format("any.tree", "nexus", False) == "nexus"


def test_resolve_tree_format_infer_nwk():
    assert resolve_tree_format("test.nwk", None, False) == "newick"


def test_resolve_tree_format_infer_newick():
    assert resolve_tree_format("test.newick", None, False) == "newick"


def test_resolve_tree_format_infer_nexus():
    assert resolve_tree_format("test.nexus", None, False) == "nexus"


def test_resolve_tree_format_unknown_extension():
    with pytest.raises(click.Abort):
        resolve_tree_format("test.txt", None, False)


def test_resolve_tree_format_debug_output():
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        resolve_tree_format("some.nwk", None, debug=True)
    mock_console.print.assert_any_call(
        f"[{STYLES['debug']}][DEBUG] Resolved tree format for 'some.nwk': newick[/{STYLES['debug']}]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_success(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["test_cmd"], returncode=0, stdout="Success output", stderr=""
    )
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        result = run_subprocess_command(
            ["test_cmd"],
            debug=False,
            success_message="Command executed successfully",
        )
    assert result is True
    mock_console.print.assert_called_once_with(
        f"[{STYLES['success']}][SUCCESS] Command executed successfully[/{STYLES['success']}]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_failure_called_process_error(mock_subproc_run):
    mock_subproc_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd=["fail_cmd_cpe"], output="out", stderr="Error output cpe"
    )
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        with pytest.raises(click.Abort):
            run_subprocess_command(
                ["fail_cmd_cpe"],
                debug=False,
                error_message_prefix="Test error CPE",
            )

    mock_console.print.assert_called_once_with(
        f"[{STYLES['error']}][ERROR] Test error CPE fail_cmd_cpe: Command '['fail_cmd_cpe']' returned non-zero exit status 1.[/{STYLES['error']}]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_file_not_found(mock_subproc_run):
    mock_subproc_run.side_effect = FileNotFoundError("Command not found")
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        with pytest.raises(click.Abort):
            run_subprocess_command(
                ["non_existent_cmd"],
                debug=False,
                error_message_prefix="FNF error",
            )

    mock_console.print.assert_called_once_with(
        f"[{STYLES['error']}][ERROR] FNF error: non_existent_cmd command not found. Please ensure it is installed and in your PATH.[/{STYLES['error']}]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_success_debug(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["debug_cmd_success"],
        returncode=0,
        stdout="Debug success output",
        stderr="Debug success stderr",
    )
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        result = run_subprocess_command(
            ["debug_cmd_success", "arg1"],
            debug=True,
            success_message="Debug success",
        )
    assert result is True
    expected_calls = [
        call(
            f"[{STYLES['debug']}][DEBUG] Running command: debug_cmd_success arg1[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['debug']}][DEBUG] debug_cmd_success stdout:\nDebug success output[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['debug']}][DEBUG] debug_cmd_success stderr:\nDebug success stderr[/{STYLES['debug']}]"
        ),
        call(f"[{STYLES['success']}][SUCCESS] Debug success[/{STYLES['success']}]"),
    ]
    mock_console.print.assert_has_calls(expected_calls, any_order=False)


@patch("subprocess.run")
def test_run_subprocess_command_success_debug_empty_stderr(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["debug_cmd_success_no_stderr"],
        returncode=0,
        stdout="Debug success output",
        stderr="",
    )
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        result = run_subprocess_command(
            ["debug_cmd_success_no_stderr"],
            debug=True,
            success_message="Debug success no stderr",
        )
    assert result is True
    expected_calls = [
        call(
            f"[{STYLES['debug']}][DEBUG] Running command: debug_cmd_success_no_stderr[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['debug']}][DEBUG] debug_cmd_success_no_stderr stdout:\nDebug success output[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['success']}][SUCCESS] Debug success no stderr[/{STYLES['success']}]"
        ),
    ]
    mock_console.print.assert_has_calls(expected_calls, any_order=False)
    # Verify no stderr output line was printed (empty stderr is skipped)
    assert not any(" stderr:\n" in c.args[0] for c in mock_console.print.call_args_list)


@patch("subprocess.run")
def test_run_subprocess_command_failure_debug(mock_subproc_run):
    cmd_list = ["debug_cmd_fail", "arg_fail"]
    mock_subproc_run.side_effect = subprocess.CalledProcessError(
        returncode=1,
        cmd=cmd_list,
        output="Debug fail stdout",
        stderr="Debug fail stderr",
    )
    mock_console = MagicMock(spec=Console)
    with patch.object(barcodeforge.utils, "console", mock_console):
        with pytest.raises(click.Abort):
            run_subprocess_command(
                cmd_list, debug=True, error_message_prefix="Debug fail error"
            )

    expected_calls_in_order = [
        call(
            f"[{STYLES['debug']}][DEBUG] Running command: {' '.join(cmd_list)}[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['error']}][ERROR] Debug fail error {cmd_list[0]}: Command '{cmd_list}' returned non-zero exit status 1.[/{STYLES['error']}]"
        ),
        call(
            f"[{STYLES['debug']}][DEBUG] {cmd_list[0]} stdout:\nDebug fail stdout[/{STYLES['debug']}]"
        ),
        call(
            f"[{STYLES['debug']}][DEBUG] {cmd_list[0]} stderr:\nDebug fail stderr[/{STYLES['debug']}]"
        ),
    ]
    mock_console.print.assert_has_calls(expected_calls_in_order, any_order=False)
