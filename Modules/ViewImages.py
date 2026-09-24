"""Image viewing module."""

import math
import warnings
from logging import getLogger
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from matplotlib import patches

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def plot_embedding_with_legend(adata, embedding_key, color_col, module_dir):
    """Plot an embedding colored by a single column with the legend outside the axes."""
    is_categorical = (
        color_col in adata.obs.columns
        and not pd.api.types.is_numeric_dtype(adata.obs[color_col])
    )

    if not is_categorical:
        fig, ax = plt.subplots(figsize=(7, 6), facecolor="white")
        sc.pl.embedding(
            adata,
            basis=embedding_key,
            color=color_col,
            ax=ax,
            show=False,
            frameon=False,
        )
    else:
        n_categories = adata.obs[color_col].astype("category").cat.categories.size

        # Keep at most ~25 legend rows; add legend columns beyond that
        max_rows = 25
        ncol = max(1, math.ceil(n_categories / max_rows))
        n_rows = math.ceil(n_categories / ncol)
        fontsize = 8 if n_categories > 30 else 10

        # Legend sits outside the axes; bbox_inches="tight" expands the saved
        # image to include it. Only grow the height for very long legends.
        side = max(6, n_rows * 0.28)
        fig, ax = plt.subplots(figsize=(side, side), facecolor="white")

        sc.pl.embedding(
            adata,
            basis=embedding_key,
            color=color_col,
            ax=ax,
            show=False,
            frameon=False,
            legend_loc=None,
        )

        categories = adata.obs[color_col].cat.categories
        colors = adata.uns.get(f"{color_col}_colors", [])
        handles = [
            plt.Line2D(
                [], [], marker="o", linestyle="", color=color, markersize=6, label=cat
            )
            for cat, color in zip(categories, colors)
        ]
        ax.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            ncol=ncol,
            frameon=False,
            fontsize=fontsize,
            title=color_col,
            title_fontsize=fontsize + 1,
            borderaxespad=0,
        )

    umap_out = module_dir / f"embedding_{embedding_key}_{color_col}.png"
    fig.savefig(umap_out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved embedding plot to {umap_out}")


def run_view_images(
    data_type: str,
    input_adata_path: Union[str, Path],
    sample_key: str,
    module_dir: Path,
    gene_list: List[str],
    cluster_name: str,
    n_grid_x: int = 10,
    n_grid_y: int = 10,
    embedding_key: Optional[str] = None,
    umap_color_columns: Optional[Union[str, List[str]]] = None,
) -> Optional[Path]:
    """Run the image viewing module robustly."""

    spatial_key = "global" if data_type == "CosMx" else "spatial"
    grid_csv_path = None

    module_dir.mkdir(exist_ok=True)
    sc.settings.figdir = module_dir

    logger.info("Loading data...")
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    # 1. Embedding plots
    if embedding_key is not None and embedding_key in adata.obsm:
        if umap_color_columns is None:
            umap_color_columns = []
        elif isinstance(umap_color_columns, str):
            umap_color_columns = [umap_color_columns]

        # One figure per color column so each legend has room to fit
        for color_col in umap_color_columns:
            if color_col not in adata.obs.columns and color_col not in adata.var_names:
                logger.warning(
                    f"'{color_col}' not found in adata.obs or adata.var_names. Skipping."
                )
                continue

            logger.info(f"Plotting embedding: {embedding_key} colored by {color_col}")
            plot_embedding_with_legend(adata, embedding_key, color_col, module_dir)
    elif embedding_key is not None:
        logger.warning(
            f"Embedding '{embedding_key}' not found in adata.obsm. Skipping embedding plot."
        )

    # 2. View plots (Tissue Clusters)
    logger.info("Visualize clusters on tissue...")

    sq.pl.spatial_scatter(
        adata,
        spatial_key=spatial_key,
        library_key=sample_key,
        shape=None,
        outline=False,
        color=[cluster_name, "total_counts"],
        wspace=0.4,
        size=1,
        dpi=300,
        figsize=(10, 6),
    )

    plt.savefig(
        module_dir / "leiden_clusters_all_samples.png",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close()
    logger.info(
        f"Saved leiden clusters plot to {module_dir / 'leiden_clusters_all_samples.png'}"
    )

    # 3. View plots (Leiden Clusters with FOV Bounding Boxes)
    fov_col = "fov" if "fov" in adata.obs.columns else None

    if fov_col:
        logger.info("Plotting clusters with FOV bounding boxes...")

        # Loop over samples to draw rectangles properly
        for sample in adata.obs[sample_key].unique():
            adata_sample = adata[adata.obs[sample_key] == sample].copy()

            fig_fov, ax_fov = plt.subplots(figsize=(10, 10), facecolor="white")

            sq.pl.spatial_scatter(
                adata_sample,
                spatial_key=spatial_key,
                shape=None,
                outline=False,
                color=cluster_name,
                size=1,
                dpi=300,
                ax=ax_fov,
                fig=fig_fov,
            )

            # Extract spatial coordinates and create a dataframe for easy grouping
            coords = adata_sample.obsm[spatial_key]
            df_coords = pd.DataFrame(
                {
                    "x": coords[:, 0],
                    "y": coords[:, 1],
                    "fov": adata_sample.obs[fov_col].values,
                }
            )

            # Calculate min/max for each FOV and draw boxes
            for fov_id, group in df_coords.groupby("fov"):
                x_min, x_max = group["x"].min(), group["x"].max()
                y_min, y_max = group["y"].min(), group["y"].max()

                width = x_max - x_min
                height = y_max - y_min

                rect = patches.Rectangle(
                    (x_min, y_min),
                    width,
                    height,
                    linewidth=1.2,
                    edgecolor="black",
                    facecolor="none",
                    linestyle="solid",
                )
                ax_fov.add_patch(rect)

                ax_fov.text(
                    x_min + (width / 2),
                    y_min + (height / 2),
                    str(fov_id),
                    color="black",
                    fontsize=12,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    bbox=dict(
                        boxstyle="circle,pad=0.3",
                        facecolor="white",
                        edgecolor="black",
                        alpha=0.8,
                    ),
                )

            fig_fov.savefig(
                module_dir / f"FOV_mapping_{sample}.png",
                dpi=300,
                facecolor="white",
                bbox_inches="tight",
            )
            plt.close(fig_fov)
            logger.info(
                f"Saved FOV clusters plot to {module_dir / f'FOV_mapping_{sample}.png'}"
            )
    else:
        logger.warning(
            "No 'fov' or 'FOV' column found in adata.obs. Skipping FOV bounding boxes plot."
        )

    # 4. ROI Grid generation (Xenium)
    if data_type == "Xenium":
        logger.info("Generating spatial ROI grid for Xenium data...")

        for sample in adata.obs[sample_key].unique():
            adata_sample = adata[adata.obs[sample_key] == sample].copy()
            coords = adata_sample.obsm[spatial_key]

            # Defensive calculation against NaNs
            x_min_global = np.nanmin(coords[:, 0])
            x_max_global = np.nanmax(coords[:, 0])
            y_min_global = np.nanmin(coords[:, 1])
            y_max_global = np.nanmax(coords[:, 1])

            x_step = (x_max_global - x_min_global) / n_grid_x
            y_step = (y_max_global - y_min_global) / n_grid_y

            box_data = []
            box_id = 1

            for row in range(n_grid_y):
                for col in range(n_grid_x):
                    x_min = x_min_global + col * x_step
                    x_max = x_min_global + (col + 1) * x_step
                    y_max = y_max_global - row * y_step
                    y_min = y_max_global - (row + 1) * y_step

                    box_data.append(
                        {
                            "Box_ID": box_id,
                            "x_min": x_min,
                            "x_max": x_max,
                            "y_min": y_min,
                            "y_max": y_max,
                        }
                    )
                    box_id += 1

            fig_grid, ax_grid = plt.subplots(figsize=(10, 10), facecolor="white")

            sq.pl.spatial_scatter(
                adata_sample,
                spatial_key=spatial_key,
                shape=None,
                outline=False,
                color=cluster_name,
                size=1,
                dpi=300,
                ax=ax_grid,
                fig=fig_grid,
            )

            for box in box_data:
                rect = patches.Rectangle(
                    (box["x_min"], box["y_min"]),
                    box["x_max"] - box["x_min"],
                    box["y_max"] - box["y_min"],
                    linewidth=1.2,
                    edgecolor="black",
                    facecolor="none",
                    linestyle="solid",
                )
                ax_grid.add_patch(rect)

                ax_grid.text(
                    (box["x_min"] + box["x_max"]) / 2,
                    (box["y_min"] + box["y_max"]) / 2,
                    str(box["Box_ID"]),
                    color="black",
                    fontsize=8,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    bbox=dict(
                        boxstyle="circle,pad=0.3",
                        facecolor="white",
                        edgecolor="black",
                        alpha=0.8,
                    ),
                )

            grid_plot_path = module_dir / f"Xenium_ROI_grid_{sample}.png"
            fig_grid.savefig(
                grid_plot_path, dpi=300, facecolor="white", bbox_inches="tight"
            )
            plt.close(fig_grid)

            # Save coordinates to CSV PER SLIDE so MuSpAn can find them
            df_grid = pd.DataFrame(box_data).round(2)
            grid_csv_path = module_dir / f"Xenium_ROI_grid_coordinates_{sample}.csv"
            df_grid.to_csv(grid_csv_path, index=False)
            logger.info(f"Saved Xenium ROI grid coordinates to {grid_csv_path}")

    # 5. View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    for sample in adata.obs[sample_key].unique():
        logger.info(f"Plotting genes for sample: {sample}...")

        adata_sample = adata[adata.obs[sample_key] == sample].copy()
        n_genes = len(gene_list)

        ncols = math.ceil(math.sqrt(n_genes))
        nrows = math.ceil(n_genes / ncols)

        fig, axes = plt.subplots(
            nrows, ncols, figsize=(5 * ncols, 5 * nrows), facecolor="white"
        )

        if n_genes == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for ax in axes[n_genes:]:
            ax.remove()

        axes = axes[:n_genes]

        sq.pl.spatial_scatter(
            adata_sample,
            spatial_key=spatial_key,
            color=gene_list,
            use_raw=True,
            shape=None,
            size=2,
            img=False,
            ax=axes,
            fig=fig,
        )

        out_filename = module_dir / f"gene_expression_{sample}.png"
        fig.savefig(out_filename, dpi=300, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Saved gene expression plot to {out_filename}")

    # Save anndata object
    out_path = module_dir / input_adata_path.name
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
    logger.info("Imaging module completed successfully.")

    return grid_csv_path
