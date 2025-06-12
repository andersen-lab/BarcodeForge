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


def test_plot_barcode_matplotlib_creates_image_file(
    df_refposalt_long_for_plotter: pd.DataFrame, tmp_path: Path
):
    """
    Tests if plot_barcode_altair (matplotlib version) creates an output image file
    when given data it can parse.
    """
    output_file = tmp_path / "test_plot.png"
    plot_barcode_altair(df_refposalt_long_for_plotter, str(output_file))
    assert output_file.exists()
    assert output_file.stat().st_size > 0  # Check if file is not empty
