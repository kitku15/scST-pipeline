import warnings
from logging import getLogger

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import numpy as np
import squidpy as sq
import math
import hashlib
import scanpy as sc
from pathlib import Path

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def get_dynamic_palette(adata, column_name):
    """
    Automatically assigns consistent colors to categories in adata.obs[column_name].
    Ensures 'Unknown' is always grey, and colors remain consistent across datasets.
    """
    # Get unique categories (handling pandas Categorical safely)
    categories = adata.obs[column_name].unique()

    # Use ONE large palette so the color pool never changes
    base_colors = list(sc.pl.palettes.default_102)
    palette_dict = {}

    for cat in categories:
        # Handle actual NaN values safely
        if pd.isna(cat):
            palette_dict[cat] = "#d3d3d3"
            continue

        cat_str = str(cat)

        if cat_str.lower() in ["unknown", "low confidence", "nan"]:
            palette_dict[cat] = "#d3d3d3"  # Light Grey
        else:
            # Create a deterministic integer from the category name string
            # We use MD5 because Python's built-in hash() changes every time you restart Python
            hash_int = int(hashlib.md5(cat_str.encode("utf-8")).hexdigest(), 16)

            # Map that integer to our color palette
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


def plot_spatialplotsplit(adata, spatial_key, celltype_col, output_file, n_cols=4):
    """
    Makes a spatial plot for all celltypes highlighting regions of that specific cell type
    Makes it easy to see spatial distribution of Cell Types
    Can adjust number of columns
    celltype_col -> choose which cell type column to color by
    """
    # 1. Get all unique cell types from the CellTypist column
    cell_types = adata.obs[celltype_col].astype(str).unique()
    cell_types = sorted(cell_types)

    # 2. Set up the matplotlib grid
    n_types = len(cell_types)
    n_rows = math.ceil(n_types / n_cols)

    # Create the figure and axes
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(n_cols * 5, n_rows * 2.5), dpi=300
    )
    axes = axes.flatten()

    # 3. Loop through each cell type and plot
    for i, cell_type in enumerate(cell_types):
        ax = axes[i]

        # Check against "majority_voting"
        adata.obs["highlight"] = np.where(
            adata.obs[celltype_col].astype(str) == cell_type, cell_type, "other"
        )

        adata.obs["highlight"] = adata.obs["highlight"].astype("category")
        adata.obs["highlight"] = adata.obs["highlight"].cat.reorder_categories(
            [cell_type, "other"], ordered=True
        )

        if "highlight_colors" in adata.uns:
            del adata.uns["highlight_colors"]

        # Plot on the specific axis
        sq.pl.spatial_scatter(
            adata,
            color="highlight",
            spatial_key=spatial_key,
            shape=None,
            size=2,
            frameon=False,
            img=False,
            outline=False,
            palette=ListedColormap(["red", "lightgrey"]),
            ax=ax,
            title=f"{cell_type}",
        )

    # 4. Clean up any empty subplots
    for i in range(n_types, len(axes)):
        axes[i].set_visible(False)
        axes[i].axis("off")

    # 5. Finalize and show
    plt.tight_layout()
    plt.savefig(output_file)

    # adata.obs.drop(columns=["highlight"], inplace=True)
    # adata.obs.drop(columns=["highlight_colors"], inplace=True)


def plot_umapspatialscatter(adata, spatial_key, celltype_col, umap_col, output_path):
    """
    This plots a spatial scatter plot next to a UMAP in the same figure
    """
    # Apply color mappping that ensures in every single plot each cell type is the same color
    # Generate the palette based on the results
    palette_map = get_dynamic_palette(adata, celltype_col)

    # Apply to adata.uns so Scanpy/Squidpy uses it automatically
    ordered_categories = adata.obs[celltype_col].cat.categories
    adata.uns[f"{celltype_col}_colors"] = [
        palette_map[cat] for cat in ordered_categories
    ]

    logger.info(f"Generating side-by-side plots for {celltype_col}...")

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10), dpi=300)

    # 1. Spatial Scatter Plot (Squidpy)
    # Note: 'show' argument is REMOVED here
    sq.pl.spatial_scatter(
        adata,
        color=celltype_col,
        spatial_key=spatial_key,
        library_id=None,  # This helps with the WARNING you saw
        ax=ax1,
        shape=None,
        img=False,
        size=1,
        alpha=0.8,
        frameon=False,
        legend_fontsize="x-small",
        title=f"Spatial: {celltype_col}",
    )

    # 2. UMAP Plot (Scanpy)
    # Note: Scanpy DOES use 'show=False'
    sc.pl.embedding(
        adata,
        basis=umap_col,
        color=celltype_col,
        ax=ax2,
        show=False,
        frameon=False,
        title=f"UMAP: {celltype_col}",
    )

    # Final touches
    plt.tight_layout()

    # Save the combined figure
    # cant use the same syntax as scpy but have to write full path from root dir
    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"🎉 Saved combined plot to: {output_path}")


def run_CellType_plotting(
    method,
    datatype,
    module_dir,
    adata,
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

    plot_spatialplotsplit(
        adata,
        spatial_key,
        celltype_col,
        output_file=f"{fig_dir}/spatialplotsplit_{celltype_col}.png",
    )
    plot_umapspatialscatter(
        adata,
        spatial_key,
        celltype_col,
        umap_col,
        output_path=f"{fig_dir}/umapspatialscatter_{celltype_col}.png",
    )
