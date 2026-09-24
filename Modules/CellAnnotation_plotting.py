import hashlib
import math
import warnings
from logging import getLogger
from pathlib import Path
from typing import Dict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from matplotlib.colors import ListedColormap

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def get_dynamic_palette(adata: AnnData, column_name: str) -> Dict[str, str]:
    """Generates a consistent color palette based on string hashing."""
    categories = adata.obs[column_name].unique()
    base_colors = list(sc.pl.palettes.default_102)
    palette_dict = {}

    for cat in categories:
        if pd.isna(cat):
            palette_dict[cat] = "#d3d3d3"
            continue

        cat_str = str(cat)
        if cat_str.lower() in ["unknown", "low confidence", "nan"]:
            palette_dict[cat] = "#d3d3d3"
        else:
            hash_int = int(hashlib.md5(cat_str.encode("utf-8")).hexdigest(), 16)
            palette_dict[cat] = base_colors[hash_int % len(base_colors)]

    return palette_dict


def plot_celltypecomp(
    adata, cluster_col, celltype_col, output_file, min_frac_threshold=0.05
):
    """
    Helpful plot that shows you what individual cells in each cluster were assigned

    - min_frac_threshold = number between 0-1, cell types that take up less than this
    threshold in a cluster won't be included. Default set to 5% (0.05)

    - cluster_col: column for cluster
    - celltype_col: column for individual cell type annotation (CellTypist has this, need
    to check with ScType)
    """

    cell_level_preds = adata.obs[
        celltype_col
    ]  # this is the Cell Type for individual-cell cell type
    clusters = adata.obs[cluster_col]

    # Build composition matrix
    composition_matrix = pd.crosstab(clusters, cell_level_preds, normalize="index")

    # Filter out values < threshold
    filtered_matrix = composition_matrix.where(
        composition_matrix >= min_frac_threshold, 0
    )

    # Renormalise rows
    filtered_matrix = filtered_matrix.div(filtered_matrix.sum(axis=1), axis=0)

    # Drop columns that are now all zero (removes from legend too)
    filtered_matrix = filtered_matrix.loc[:, (filtered_matrix > 0).any(axis=0)]

    # Generate enough distinct colors (UPDATED for Matplotlib 3.7+)
    n_types = filtered_matrix.shape[1]
    cmap = mpl.colormaps["tab20"]  # Modern way to call the colormap
    colors = cmap(np.linspace(0, 1, n_types))

    # Plot
    _ = filtered_matrix.plot(kind="bar", stacked=True, figsize=(14, 6), color=colors)

    plt.legend(
        bbox_to_anchor=(1.01, 1), loc="upper left", title="Cell-Level Predictions"
    )

    plt.title(
        f"Distribution of Cell-Level Predictions within Clusters (>={min_frac_threshold * 100}%)"
    )
    plt.xlabel(f"Cluster ({cluster_col})")
    plt.ylabel("Fraction of Cells")
    plt.tight_layout()
    plt.savefig(output_file)


def summary_celltypecomp(
    adata, cluster_col, celltype_col, output_file, min_frac_threshold=0.05
):
    """
    Helpful text summary that shows you what individual cells in each cluster were assigned
    Useful when theres so many cell types the plot looks confusing

    - min_frac_threshold = number between 0-1, cell types that take up less than this
    threshold in a cluster won't be included. Default set to 5% (0.05)

    - cluster_col: column for cluster
    - celltype_col: column for individual cell type annotation (CellTypist has this, need
    to check with ScType)
    """
    cell_level_preds = adata.obs[celltype_col]
    clusters = adata.obs[cluster_col]

    # Proportion matrix (clusters × cell types)
    composition_matrix = pd.crosstab(clusters, cell_level_preds, normalize="index")

    # Open the file in write mode
    with open(output_file, "w", encoding="utf-8") as f:
        for cluster in composition_matrix.index:
            row = composition_matrix.loc[cluster]

            # majority cell type
            majority = row.idxmax()
            majority_frac = row.max()

            # filter by threshold
            filtered = row[row >= min_frac_threshold].sort_values(ascending=False)

            # Redirect prints to the file using file=f
            print("\n" + "=" * 60, file=f)
            print(f"Cluster: {cluster}", file=f)
            print(f"Majority cell type: {majority} ({majority_frac:.3f})", file=f)
            print("-" * 60, file=f)
            print(f"Cell types ≥ {min_frac_threshold} fraction:", file=f)

            if filtered.empty:
                print("  None above threshold", file=f)
            else:
                for ct, frac in filtered.items():
                    print(f"  - {ct}: {frac:.3f}", file=f)

        print("\n" + "=" * 60, file=f)
        print("Done", file=f)


def plot_spatialplotsplit(
    adata: AnnData,
    spatial_key: str,
    sample_key: str,
    celltype_col: str,
    output_file: str,
    n_cols: int = 4,
) -> None:
    """Creates split spatial plots per celltype."""
    out_path = Path(output_file)

    for sample in adata.obs[sample_key].unique():
        adata_sample = adata[adata.obs[sample_key] == sample].copy()
        cell_types = sorted(adata.obs[celltype_col].astype(str).unique())

        n_types = len(cell_types)
        n_rows = math.ceil(n_types / n_cols)

        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(n_cols * 5, n_rows * 2.5), dpi=300
        )
        axes = np.atleast_1d(axes).flatten()  # Ensure safety for 1x1 subplots

        for i, cell_type in enumerate(cell_types):
            ax = axes[i]
            adata_sample.obs["highlight"] = np.where(
                adata_sample.obs[celltype_col].astype(str) == cell_type,
                cell_type,
                "other",
            )
            adata_sample.obs["highlight"] = adata_sample.obs["highlight"].astype(
                "category"
            )
            adata_sample.obs["highlight"] = adata_sample.obs[
                "highlight"
            ].cat.set_categories([cell_type, "other"], ordered=True)

            if "highlight_colors" in adata_sample.uns:
                del adata_sample.uns["highlight_colors"]

            sq.pl.spatial_scatter(
                adata_sample,
                color="highlight",
                spatial_key=spatial_key,
                shape=None,
                size=2,
                frameon=False,
                img=False,
                outline=False,
                palette=ListedColormap(["red", "lightgrey"]),
                ax=ax,
                title=f"{cell_type} (Sample: {sample})",
            )

        # Safely hide empty subplots without IndexError
        for i in range(n_types, len(axes)):
            if i < len(axes):
                axes[i].set_visible(False)
                axes[i].axis("off")

        plt.tight_layout()
        dynamic_out_file = (
            out_path.parent / f"{out_path.stem}_{sample}{out_path.suffix}"
        )
        plt.savefig(dynamic_out_file)
        plt.close(fig)

        logger.info(f"Saved spatial split plot for {sample} to {dynamic_out_file}")


def plot_umapspatialscatter(
    adata, spatial_key, sample_key, celltype_col, umap_col, output_path
):
    """
    This plots a spatial scatter plot next to a UMAP in the same figure.
    Creates ONE separate image file per sample.
    """

    if adata.obs[celltype_col].isna().any():
        logger.info(
            f"Filling NaN values in '{celltype_col}' with 'Unknown' to prevent plotting errors."
        )

        # Add "Unknown" to categories if it doesn't exist
        if "Unknown" not in adata.obs[celltype_col].cat.categories:
            adata.obs[celltype_col] = adata.obs[celltype_col].cat.add_categories(
                "Unknown"
            )

        # Replace the NaN values
        adata.obs[celltype_col] = adata.obs[celltype_col].fillna("Unknown")

    color_key = f"{celltype_col}_colors"
    if color_key not in adata.uns:
        logger.info(
            f"No colors found for {celltype_col}. Generating dynamic palette..."
        )
        palette_map = get_dynamic_palette(adata, celltype_col)
        ordered_categories = adata.obs[celltype_col].cat.categories
        adata.uns[color_key] = [palette_map[cat] for cat in ordered_categories]
    else:
        logger.info(f"Using existing colors for {celltype_col}")

    out_path = Path(output_path)

    logger.info(f"Generating side-by-side plots for {celltype_col} per sample...")

    for sample in adata.obs[sample_key].unique():
        adata_sample = adata[adata.obs[sample_key] == sample].copy()

        # Create the figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10), dpi=300)

        # 1. Spatial Scatter Plot (Just for this specific sample)
        sq.pl.spatial_scatter(
            adata_sample,
            color=celltype_col,
            spatial_key=spatial_key,
            ax=ax1,
            shape=None,
            img=False,
            size=1,
            alpha=0.8,
            frameon=False,
            legend_fontsize="x-small",
            title=f"Spatial: {celltype_col} ({sample})",
        )

        # 2. UMAP Plot (We use the GLOBAL adata here so we see all cells for context)
        sc.pl.embedding(
            adata,
            basis=umap_col,
            color=celltype_col,
            ax=ax2,
            show=False,
            frameon=False,
            title=f"Global UMAP: {celltype_col}",
        )

        plt.tight_layout()

        # Dynamic filename per sample
        dynamic_out_file = (
            out_path.parent / f"{out_path.stem}_{sample}{out_path.suffix}"
        )
        plt.savefig(dynamic_out_file, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved combined plots to: {out_path.parent}")


def run_CellType_plotting(
    method,
    datatype,
    module_dir,
    adata,
    sample_key,
    cluster_col,
    celltype_col,
    umap_col,
    indivcellanno_col=None,
):
    """
    Makes 4 figures related to cell type annotation
    """

    if datatype == "CosMx":
        spatial_key = "global"
    elif datatype == "Xenium":
        spatial_key = "spatial"

    fig_dir = Path(module_dir) / method / cluster_col
    fig_dir.mkdir(parents=True, exist_ok=True)

    if indivcellanno_col:
        plot_celltypecomp(
            adata,
            cluster_col,
            indivcellanno_col,
            output_file=f"{fig_dir}/celltypecomp_{celltype_col}.png",
        )
        summary_celltypecomp(
            adata,
            cluster_col,
            indivcellanno_col,
            output_file=f"{fig_dir}/celltypecomp_{celltype_col}.txt",
        )
    # temporary disabling spatial split plots for time save
    # plot_spatialplotsplit(
    #     adata,
    #     spatial_key,
    #     sample_key,
    #     celltype_col,
    #     output_file=f"{fig_dir}/spatialplotsplit_{celltype_col}.png",
    # )
    plot_umapspatialscatter(
        adata,
        spatial_key,
        sample_key,
        celltype_col,
        umap_col,
        output_path=f"{fig_dir}/umapspatialscatter_{celltype_col}.png",
    )
