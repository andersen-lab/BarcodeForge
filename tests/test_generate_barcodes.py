import pytest
import pandas as pd
from barcodeforge.generate_barcodes import (
    parse_tree_paths,
    convert_to_barcodes,
    reversion_checking,
    check_no_flip_pairs,
    identify_chains,
    check_mutation_chain,
    replace_underscore_with_dash,
    create_barcodes_from_lineage_paths,
)
from barcodeforge.utils import sortFun  # Assuming sortFun is in utils


# Sample data for testing
@pytest.fixture
def sample_lineage_data():
    data = {"clade": ["A", "B"], "from_tree_root": [">T123C>G456A", ">C789T"]}
    return pd.DataFrame(data)


@pytest.fixture
def sample_barcode_data():
    data = {
        "T123C": [1, 0],
        "G456A": [1, 0],
        "C789T": [0, 1],
        "A123T": [0, 0],  # For reversion check
    }
    df = pd.DataFrame(data, index=["A", "B"])
    return df


def test_parse_tree_paths(sample_lineage_data):
    parsed_df = parse_tree_paths(sample_lineage_data)
    assert not parsed_df.empty
    assert "from_tree_root" in parsed_df.columns
    assert isinstance(parsed_df.loc["A", "from_tree_root"], list)
    assert parsed_df.loc["A", "from_tree_root"] == ["T123C", "G456A"]


def test_convert_to_barcodes(sample_lineage_data):
    parsed_df = parse_tree_paths(sample_lineage_data)
    barcodes_df = convert_to_barcodes(parsed_df)
    assert not barcodes_df.empty
    assert "T123C" in barcodes_df.columns
    assert barcodes_df.loc["A", "T123C"] == 1
    assert barcodes_df.loc["B", "C789T"] == 1


def test_reversion_checking(sample_barcode_data):
    # Add a reversion pair
    data_reversion = {
        "T123C": [1, 0, 1],
        "C123T": [0, 0, 1],  # Reversion of T123C
        "G456A": [1, 1, 0],
    }
    df_reversion = pd.DataFrame(data_reversion, index=["A", "B", "C"])
    reversion_checked_df = reversion_checking(
        df_reversion.copy()
    )  # Use copy to avoid modifying fixture
    # After reversion checking, one of the pair should be 0 or dropped
    # This assertion depends on the specific logic of reversion_checking
    # For example, if it zeros out the minimum of the pair:
    assert (
        reversion_checked_df.loc["C", "T123C"] == 0
        or reversion_checked_df.loc["C", "C123T"] == 0
    )
    # Or if it drops columns that sum to 0 after adjustment
    assert not (
        (
            "T123C" in reversion_checked_df.columns
            and "C123T" in reversion_checked_df.columns
        )
        and (reversion_checked_df[["T123C", "C123T"]].sum(axis=1) > 0).all()
    )


def test_identify_chains(sample_barcode_data):
    # This test requires a more complex setup to reliably test chains
    # For now, a basic check that it runs and returns a list
    chains = identify_chains(sample_barcode_data)
    assert isinstance(chains, list)


def test_check_mutation_chain(sample_barcode_data):
    # Similar to identify_chains, needs careful setup
    # Basic check that it runs and returns a DataFrame
    chained_df = check_mutation_chain(sample_barcode_data.copy())  # Use copy
    assert isinstance(chained_df, pd.DataFrame)


def test_replace_underscore_with_dash():
    data = {"value": [1, 2]}
    df = pd.DataFrame(data, index=["lineage_A", "lineage_B"])
    replaced_df = replace_underscore_with_dash(df)
    assert "lineage-A" in replaced_df.index
    assert "lineage-B" in replaced_df.index


@pytest.fixture
def temp_barcode_file(tmp_path):
    file_path = tmp_path / "test_barcodes.csv"
    data = {"mut1": [1, 0], "mut2": [0, 1]}
    df = pd.DataFrame(data, index=["L1", "L2"])
    df.to_csv(file_path)
    return file_path


def test_test_no_flip_pairs_no_flips(temp_barcode_file):
    # This should pass without raising an exception
    try:
        check_no_flip_pairs(str(temp_barcode_file))  # Renamed from test_no_flip_pairs
    except Exception as e:
        pytest.fail(f"check_no_flip_pairs raised an exception unexpectedly: {e}")


@pytest.fixture
def temp_barcode_file_with_flips(tmp_path):
    file_path = tmp_path / "test_barcodes_flips.csv"
    data = {"A123T": [1, 0], "T123A": [1, 0]}  # Flip pair
    df = pd.DataFrame(data, index=["L1", "L2"])
    df.to_csv(file_path)
    return file_path


def test_test_no_flip_pairs_with_flips(temp_barcode_file_with_flips):
    with pytest.raises(Exception, match=r"FAIL: flip pairs found"):
        check_no_flip_pairs(
            str(temp_barcode_file_with_flips)
        )  # Renamed from test_no_flip_pairs


@pytest.fixture
def temp_lineage_paths_file(tmp_path):
    file_path = tmp_path / "lineage_paths.tsv"
    content = "clade\tfrom_tree_root\nLineageA\t>A123T>T456C\nLineageB\t>G789A"
    with open(file_path, "w") as f:
        f.write(content)
    return file_path


def test_create_barcodes_from_lineage_paths(temp_lineage_paths_file, tmp_path):
    output_file = tmp_path / "output_barcodes.csv"
    create_barcodes_from_lineage_paths(False, str(temp_lineage_paths_file), str(output_file))
    assert output_file.exists()
    df = pd.read_csv(output_file)
    assert not df.empty
    # Add more specific assertions based on expected output
    assert (
        "LineageA" in df.iloc[:, 0].values
    )  # Check if LineageA is in the first column (index)
