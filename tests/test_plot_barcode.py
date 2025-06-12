import pytest
import pandas as pd
from pathlib import Path
from barcodeforge.plot_barcode import plot_barcode_altair, create_barcode_plot
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

@pytest.fixture
def df_refposalt_long_for_plotter():
    """
    Provides a long-format DataFrame where 'Mutation' is in 'RefPosAlt' format,
    suitable for direct input to the plot_barcode_altair (matplotlib) function.
    Includes cases that will be filtered (z=0) and kept.
    """
    data = {
        "Lineage": ["L1", "L1", "L2", "L2", "L3", "L3", "L1", "L2"],
        "Mutation": ["A1T", "C22G", "A1T", "G333C", "C22G", "G333C", "T4A", "T4A"],
        "z": [1, 0, 1, 1, 0, 1, 1, 0],  # z=0 for L1/C22G and L3/C22G, L2/T4A
    }
    return pd.DataFrame(data)


@pytest.fixture
def csv_generic_mutations_file(tmp_path: Path) -> Path:
    """
    Provides a CSV file with generic mutation names (e.g., "M1", "M2").
    This will cause issues in the current plot_barcode_altair if not handled,
    as str.extract will not find RefPosAlt.
    """
    file_path = tmp_path / "sample_barcodes_generic.csv"
    data = {
        "M1": {"L1": 1, "L2": 1, "L3": 0},
        "M2": {"L1": 0, "L2": 0, "L3": 1},
        "M3": {"L1": 0, "L2": 1, "L3": 1},
    }
    df = pd.DataFrame(data)
    df.index.name = "Lineage"
    df.to_csv(file_path)
    return file_path


@pytest.fixture
def csv_refposalt_mutations_file(tmp_path: Path) -> Path:
    """
    Provides a CSV file where column headers (mutations) are in 'RefPosAlt' format.
    Suitable for e2e testing create_barcode_plot with the matplotlib plotter.
    """
    file_path = tmp_path / "sample_barcodes_refposalt.csv"
    data = {
        "A1T": {"L1": 1, "L2": 1, "L3": 0},
        "C22G": {"L1": 1, "L2": 0, "L3": 1}, # L1/C22G is present
        "G333C": {"L1": 0, "L2": 1, "L3": 1},
    }
    df = pd.DataFrame(data)
    df.index.name = "Lineage"
    df.to_csv(file_path)
    return file_path


def test_plot_barcode_matplotlib_creates_image_file(df_refposalt_long_for_plotter: pd.DataFrame, tmp_path: Path):
    """
    Tests if plot_barcode_altair (matplotlib version) creates an output image file
    when given data it can parse.
    """
    output_file = tmp_path / "test_plot.png"
    plot_barcode_altair(df_refposalt_long_for_plotter, str(output_file))
    assert output_file.exists()
    assert output_file.stat().st_size > 0  # Check if file is not empty


def test_plot_barcode_matplotlib_uses_correct_calls(df_refposalt_long_for_plotter: pd.DataFrame, tmp_path: Path, mocker):
    """
    Tests if plot_barcode_altair (matplotlib version) makes the expected calls
    to matplotlib and seaborn.
    """
    output_file = tmp_path / "test_plot_structure.png"

    mock_subplots = mocker.patch("matplotlib.pyplot.subplots")
    mock_heatmap = mocker.patch("seaborn.heatmap")
    mock_savefig = mocker.patch("matplotlib.pyplot.savefig")
    mock_close = mocker.patch("matplotlib.pyplot.close")

    plot_barcode_altair(df_refposalt_long_for_plotter, str(output_file))

    mock_subplots.assert_called_once()
    mock_heatmap.assert_called_once()
    mock_savefig.assert_called_once_with(str(output_file), bbox_inches="tight")
    mock_close.assert_called_once()

    heatmap_args, heatmap_kwargs = mock_heatmap.call_args
    heatmap_data = heatmap_args[0] # First positional argument is the data
    assert isinstance(heatmap_data, pd.DataFrame)
