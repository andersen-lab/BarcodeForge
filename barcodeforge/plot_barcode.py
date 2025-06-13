"""Plot barcode from CSV file."""

import pandas as pd
from rich.console import Console

from .utils import sortFun, STYLES
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import math

console = Console()


def create_barcode_visualization(
    barcode_df_long: pd.DataFrame, chunk_size: int, output_path: str
) -> None:
    """
    Generates and saves an optimized seaborn plot from a long-format barcode DataFrame.
    """
    # Filter, extract, pivot, and reshape
    df = barcode_df_long[barcode_df_long.z.ne(0)].copy()
    df[["Reference", "pos", "alt"]] = df.Mutation.str.extract(
        r"([A-Za-z]+)(\d+)([A-Za-z]+)"
    )
    df.pos = df.pos.astype(int)
    wide = (
        df.drop(columns=["Mutation", "z"])
        .pivot(index=["pos", "Reference"], columns="Lineage", values="alt")
        .reset_index()
    )
    plot_df = wide.melt(id_vars=["pos"], var_name="Lineage", value_name="alt").fillna(
        {"alt": "Unchanged"}
    )

    # Build color palettes once
    base_palette = {"A": "#CC3311", "C": "#33BBEE", "G": "#EE7733", "T": "#009988"}
    # add any new or "Unchanged"
    for b in plot_df.alt.unique():
        base_palette.setdefault(b, "#555555")
    base_palette["Unchanged"] = "#BBBBBB"
    base_list = sorted(base_palette)
    cmap = ListedColormap([base_palette[b] for b in base_list])
    lineage_colors = ["#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA"]

    all_pos = sorted(plot_df.pos.unique())
    CHUNK = chunk_size if chunk_size > 0 else len(all_pos)
    num_chunks = math.ceil(len(all_pos) / CHUNK)

    # Prepare chunks
    chunks = []
    dims = []
    cat_dtype = pd.CategoricalDtype(categories=base_list, ordered=True)
    for i in range(num_chunks):
        ps = all_pos[i * CHUNK : (i + 1) * CHUNK]
        sub = plot_df[plot_df.pos.isin(ps)]
        mat = sub.pivot(index="pos", columns="Lineage", values="alt")

        # fast code array via pd.Categorical
        flat = mat.values.ravel()
        codes = pd.Categorical(flat, categories=base_list, ordered=True).codes
        code_arr = codes.reshape(mat.shape).astype(float)
        heat = pd.DataFrame(code_arr, index=mat.index, columns=mat.columns).T

        annot = mat.T.where(lambda x: x != "Unchanged", "")
        if "Reference" in heat.index:
            order = ["Reference"] + [ln for ln in heat.index if ln != "Reference"]
            heat = heat.loc[order]
            annot = annot.loc[order]

        # size
        w = len(ps) * 0.4 + mat.columns.str.len().max() * 0.2
        h = mat.columns.size * 0.25
        dims.append((w, h))
        chunks.append((heat, annot, ps))

    max_w = max(w for w, _ in dims)
    total_h = sum(h for _, h in dims) + 0.5 * (num_chunks - 1)
    fig, axes = plt.subplots(nrows=num_chunks, figsize=(max_w, total_h))
    if num_chunks == 1:
        axes = [axes]

    # Plot
    for ax, (heat, annot, ps) in zip(axes, chunks):
        lo, hi = min(ps), max(ps)
        ax.set_title(f"Positions {lo} to {hi}", loc="left", fontsize=20)
        sns.heatmap(
            heat,
            cmap=cmap,
            annot=annot,
            fmt="",
            cbar=False,
            linewidths=2,
            xticklabels=True,
            yticklabels=True,
            ax=ax,
        )

        rows, cols = heat.shape
        # grid & lineage colors
        for i in range(0, cols):
            for j in range(1, rows):
                ax.vlines(
                    x=i + 0.065,
                    ymin=j + 0.05,
                    ymax=j + 0.95,
                    color=lineage_colors[j % len(lineage_colors)],
                    linewidth=2,
                )
        for idx, lbl in enumerate(ax.get_yticklabels()):
            if idx:
                lbl.set_color(lineage_colors[idx % len(lineage_colors)])
        if "Reference" in heat.index and rows > 1:
            ax.hlines(1, -0.5, cols, color="white", linewidth=4)

        ax.set_xlabel("Genome Position")
        ax.set_ylabel("Lineage")
        # adjust tick label appearance
        ax.tick_params(axis="x", labelsize=14, labelrotation=45)
        for lbl in ax.get_xticklabels():
            lbl.set_ha("right")
        ax.tick_params(axis="y", labelsize=14)

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def create_barcode_plot(
    debug: bool, input_file_path: str, chunk_size: int, output_file_path: str
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

    create_barcode_visualization(barcode_df_long, chunk_size, output_file_path)
