import pytest
from pathlib import Path
from barcodeforge.format_tree import convert_nexus_to_newick, _remove_quotes_from_file
import tempfile

# Sample Nexus content
SAMPLE_NEXUS_CONTENT = """
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=3;
    TAXLABELS A B C;
END;
BEGIN TREES;
    TRANSLATE
        1 A,
        2 B,
        3 C;
    TREE tree1 = (1,(2,3));
END;
"""

# Expected Newick output (may vary slightly based on DendroPy's exact output)
# This is a simplified expectation. DendroPy might add branch lengths or other info.
EXPECTED_NEWICK_SIMPLE = "(A,(B,C));"

# Sample Newick content for testing reformatting
SAMPLE_NEWICK_FOR_REFORMAT_CONTENT = "(A:0.1,(B:0.2,C:0.3):0.4);"


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def test_convert_nexus_to_newick_simple(temp_dir):
    nexus_file = temp_dir / "test.nexus"
    newick_file = temp_dir / "test.nwk"

    with open(nexus_file, "w") as f:
        f.write(SAMPLE_NEXUS_CONTENT)

    convert_nexus_to_newick(
        nexus_file, newick_file, input_format="nexus", reformat_tree=False
    )

    assert newick_file.exists()
    with open(newick_file, "r") as f:
        content = f.read().strip()
        content_normalized = content.replace(":1.0", "")
        assert content_normalized.endswith(EXPECTED_NEWICK_SIMPLE)


def test_convert_nexus_to_newick_reformat(temp_dir):
    nexus_file = temp_dir / "test_reformat.nexus"
    newick_file = temp_dir / "test_reformat.nwk"

    with open(nexus_file, "w") as f:
        f.write(SAMPLE_NEXUS_CONTENT)  # Using nexus here, reformat should still work

    convert_nexus_to_newick(
        nexus_file, newick_file, input_format="nexus", reformat_tree=True
    )

    assert newick_file.exists()
    with open(newick_file, "r") as f:
        content = f.read().strip()
        assert "A" in content
        assert "B" in content
        assert "C" in content
        assert content.startswith("(") and content.endswith(");")


def test_convert_newick_input_to_newick_output(temp_dir):
    input_nwk_file = temp_dir / "input.nwk"
    output_nwk_file = temp_dir / "output.nwk"

    with open(input_nwk_file, "w") as f:
        f.write(SAMPLE_NEWICK_FOR_REFORMAT_CONTENT)

    convert_nexus_to_newick(
        input_nwk_file, output_nwk_file, input_format="newick", reformat_tree=False
    )

    assert output_nwk_file.exists()
    with open(output_nwk_file, "r") as f:
        content = f.read().strip()
        assert content == SAMPLE_NEWICK_FOR_REFORMAT_CONTENT


def test_convert_newick_input_to_newick_output_with_reformat(temp_dir):
    input_nwk_file = temp_dir / "input_reformat.nwk"
    output_nwk_file = temp_dir / "output_reformat.nwk"

    with open(input_nwk_file, "w") as f:
        f.write(SAMPLE_NEWICK_FOR_REFORMAT_CONTENT)

    convert_nexus_to_newick(
        input_nwk_file, output_nwk_file, input_format="newick", reformat_tree=True
    )

    assert output_nwk_file.exists()
    with open(output_nwk_file, "r") as f:
        content = f.read().strip()
        # Check basic structure and content
        assert "A" in content
        assert "B" in content
        assert "C" in content
        assert content.startswith("(") and content.endswith(");")


def test_input_file_not_found(temp_dir):
    non_existent_file = temp_dir / "non_existent.nexus"
    output_file = temp_dir / "output.nwk"
    with pytest.raises(FileNotFoundError):
        convert_nexus_to_newick(non_existent_file, output_file)


# Tests for _remove_quotes_from_file


def test_remove_quotes_single_and_double(temp_dir):
    file_path = temp_dir / "single_double_quotes.txt"
    original_content = '"\'A\'B"C"D"'
    expected_content = "ABCD"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_only_single(temp_dir):
    file_path = temp_dir / "only_single_quotes.txt"
    original_content = "'A'B'C'"
    expected_content = "ABC"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_only_double(temp_dir):
    file_path = temp_dir / "only_double_quotes.txt"
    original_content = '"A"B"C"'
    expected_content = "ABC"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_no_quotes(temp_dir):
    file_path = temp_dir / "no_quotes.txt"
    original_content = "ABC"
    expected_content = "ABC"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_empty_file(temp_dir):
    file_path = temp_dir / "empty_file.txt"
    original_content = ""
    expected_content = ""

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_mixed_content(temp_dir):
    file_path = temp_dir / "mixed_content.txt"
    original_content = "A'B\"C'D\"E"
    expected_content = "ABCDE"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_at_beginning_end(temp_dir):
    file_path = temp_dir / "quotes_at_edges.txt"
    original_content = "\"'ABC'\""
    expected_content = "ABC"

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content


def test_remove_quotes_only_quotes_file(temp_dir):
    file_path = temp_dir / "only_quotes_file.txt"
    original_content = "\"''\"'"
    expected_content = ""

    file_path.write_text(original_content)
    _remove_quotes_from_file(file_path)
    processed_content = file_path.read_text()
    assert processed_content == expected_content
