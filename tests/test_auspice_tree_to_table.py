import json
import pandas as pd
import pytest
from rich.console import Console
from unittest.mock import patch
import click
import barcodeforge.utils
from barcodeforge.auspice_tree_to_table import json_to_tree, process_auspice_json


@pytest.fixture
def sample_auspice_json(tmp_path):
    data = {
        "meta": {},
        "tree": {
            "name": "root",
            "node_attrs": {"div": 0.0, "country": {"value": "USA"}},
            "children": [
                {
                    "name": "A",
                    "node_attrs": {"div": 0.1, "country": {"value": "USA"}},
                },
                {
                    "name": "B",
                    "node_attrs": {"div": 0.2, "country": {"value": "CAN"}},
                },
            ],
        },
    }
    json_path = tmp_path / "tree.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    return json_path


def test_json_to_tree_basic(sample_auspice_json):
    with open(sample_auspice_json, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    tree = json_to_tree(data)
    assert tree.name == "root"
    assert len(tree.clades) == 2
    assert tree.clades[0].name == "A"
    assert pytest.approx(tree.clades[0].branch_length, rel=1e-6) == 0.1
    assert tree.clades[0].parent is tree


def test_process_auspice_json(tmp_path, sample_auspice_json):
    meta_out = tmp_path / "meta.tsv"
    tree_out = tmp_path / "tree.nwk"
    process_auspice_json(
        tree_json_path=str(sample_auspice_json),
        output_metadata_path=str(meta_out),
        output_tree_path=str(tree_out),
        include_internal_nodes=False,
        attributes=None,
    )
    assert meta_out.exists()
    assert tree_out.exists()
    df = pd.read_csv(meta_out, sep="\t")
    assert list(df["name"]) == ["A", "B"]
    assert "country" in df.columns


def test_process_auspice_json_include_internal(tmp_path, sample_auspice_json):
    meta_out = tmp_path / "meta.tsv"
    tree_out = tmp_path / "tree.nwk"
    process_auspice_json(
        tree_json_path=str(sample_auspice_json),
        output_metadata_path=str(meta_out),
        output_tree_path=str(tree_out),
        include_internal_nodes=True,
        attributes=["country"],
    )
    df = pd.read_csv(meta_out, sep="\t")
    assert set(df["name"]) == {"root", "A", "B"}
    assert list(df.columns) == ["name", "country"]


def test_process_auspice_json_missing_file(tmp_path):
    recording_console = Console(record=True)
    with patch.object(barcodeforge.utils, "console", recording_console):
        with pytest.raises(click.Abort):
            process_auspice_json(
                tree_json_path=str(tmp_path / "no.json"),
                output_metadata_path=str(tmp_path / "meta.tsv"),
                output_tree_path=None,
                include_internal_nodes=False,
                attributes=None,
            )
    assert "Tree JSON file not found" in recording_console.export_text()


def test_process_auspice_json_bad_json(tmp_path):
    bad_path = tmp_path / "bad.json"
    bad_path.write_text("not valid")
    recording_console = Console(record=True)
    with patch.object(barcodeforge.utils, "console", recording_console):
        with pytest.raises(click.Abort):
            process_auspice_json(
                tree_json_path=str(bad_path),
                output_metadata_path=None,
                output_tree_path=None,
                include_internal_nodes=False,
                attributes=None,
            )
    assert "Could not decode JSON" in recording_console.export_text()
