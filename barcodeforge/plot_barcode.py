"""Plot barcode from CSV file."""

import pandas as pd
from rich.console import Console

from .utils import sortFun, STYLES
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

console = Console()


def create_barcode_visualization(
    barcode_df_long: pd.DataFrame, output_path: str
) -> None:
    """
    Generates and saves an seaborn plot from a long-format barcode DataFrame.

    Args:
        barcode_df_long: Pandas DataFrame in long format with columns ["Lineage", "Mutation", "z"].
        output_path: Path to save the generated plot.
    """
    # Filter out rows with zero z and extract details from Mutation
    df = barcode_df_long.loc[barcode_df_long["z"] != 0].copy()
    df[["Reference", "pos", "alt"]] = df["Mutation"].str.extract(
        r"([A-Za-z]+)(\d+)([A-Za-z]+)"
    )
    df["pos"] = df["pos"].astype(int)

    # Ensure unique alt per (pos, Lineage)
    assert not df.duplicated(
        subset=["pos", "Lineage"]
    ).any(), "Each [pos, Lineage] pair must have a single alt value."

    # Pivot to wide format and reshape for plotting
    wide_df = (
        df.drop(columns=["Mutation", "z"])
        .pivot(index=["pos", "Reference"], columns="Lineage", values="alt")
        .reset_index()
    )

    plot_df = wide_df.melt(
        id_vars=["pos"], var_name="Lineage", value_name="alt"
    ).fillna({"alt": "Unchanged"})

    palette = {
        "A": "#CC3311",
        "C": "#33BBEE",
        "G": "#EE7733",
        "T": "#009988",
        "Unchanged": "#BBBBBB",  # Gray for unchanged bases
    }

    # Dynamic sizing based on unique lineages and positions
    n_positions = plot_df["pos"].nunique()
    n_lineages = plot_df["Lineage"].nunique()
    max_lineage_length = plot_df["Lineage"].str.len().max()

    console.print(
        f"[{STYLES['info']}]Creating barcode plot with {n_positions} positions and {n_lineages} lineages...[/{STYLES['info']}]"
    )

    fig, ax = plt.subplots(
        figsize=(n_positions * 0.4 + max_lineage_length * 0.2, n_lineages * 0.25)
    )

    # Create heatmap data and mapping
    matrix_df = plot_df.pivot(index="pos", columns="Lineage", values="alt")
    base_list = sorted(plot_df["alt"].unique())
    code_map = {base: idx for idx, base in enumerate(base_list)}
    heat_data = matrix_df.replace(code_map).T

    # Build a discrete colormap following the palette order from base_list
    cmap = ListedColormap([palette[base] for base in base_list])

    # Create annotations by reverting the numeric codes back to bases, omitting "Unchanged"
    annot = heat_data.replace({code: base for base, code in code_map.items()})
    annot = annot.map(lambda x: x if x != "Unchanged" else "")

    # Reorder rows so that "Reference" is at the top if present
    if "Reference" in heat_data.index:
        new_order = ["Reference"] + [
            lineage for lineage in heat_data.index if lineage != "Reference"
        ]
        heat_data = heat_data.loc[new_order]
        annot = annot.loc[new_order]

    sns.heatmap(
        heat_data,
        cmap=cmap,
        linewidths=2,
        xticklabels=True,
        yticklabels=True,
        ax=ax,
        cbar=False,
        annot=annot,
        fmt="",
    )

    lineage_palette = ["#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA"]

    # Draw vertical lines using a compact loop
    for i in range(0, heat_data.shape[1]):
        for j in range(1, heat_data.shape[0]):
            ax.vlines(
                x=i + 0.0325,
                ymin=j + 0.05,
                ymax=j + 0.95,
                color=lineage_palette[j % len(lineage_palette)],
                linewidth=2,
            )

    # set the y-tick colors to match the lineage colors
    for i in range(1, heat_data.shape[0]):
        ax.get_yticklabels()[i].set_color(lineage_palette[i % len(lineage_palette)])

    # Add extra space between "Reference" and the other lineages by drawing a thicker white line
    if "Reference" in heat_data.index and len(heat_data.index) > 1:
        ax.hlines(y=1, xmin=-0.5, xmax=heat_data.shape[1], color="white", linewidth=4)

    ax.set_xlabel("Genome Position")
    ax.set_ylabel("Lineage")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def create_barcode_plot(
    input_file_path: str, debug: bool, output_file_path: str
) -> None:
    """
    Reads a barcode CSV file, transforms it to long format, and generates a plot.

    Args:
        input_file_path: Path to the input barcode CSV file.
        output_file_path: Path to save the generated plot.
    """
    if debug:
        console.print(
            f"[{STYLES['info']}]Reading barcode data from: {input_file_path}[/{STYLES['info']}]"
        )
    barcode_df = pd.read_csv(input_file_path, header=0, index_col=0)

    if debug:
        console.print(
            f"[{STYLES['debug']}]Barcode DataFrame shape: {barcode_df.shape}[/{STYLES['debug']}]"
        )
        console.print(
            f"[{STYLES['debug']}]Barcode DataFrame columns: {barcode_df.columns.tolist()}[/{STYLES['debug']}]"
        )
        console.print(
            f"[{STYLES['info']}]Transforming barcode data to long format...[/{STYLES['info']}]"
        )
    barcode_df_long = barcode_df.stack().reset_index()
    barcode_df_long.columns = ["Lineage", "Mutation", "z"]

    create_barcode_visualization(barcode_df_long, output_file_path)
