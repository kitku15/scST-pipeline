"""Image viewing module."""

import math
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq
from matplotlib import patches

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_view_images(
    data_type,
    input_adata_path,
    sample_key,
    module_dir,
    gene_list,
    cluster_name,
    n_grid_x=10,
    n_grid_y=10,
    embedding_key=None,
    umap_color_columns=None,
):
    """Run the image viewing module."""

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    grid_csv_path = None

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    sc.settings.figdir = module_dir  # set the figures dir to not be figures

    # Import data
    logger.info("Loading data...")
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    if embedding_key in adata.obsm:
        logger.info(
            f"Plotting embedding: {embedding_key} colored by {umap_color_columns}"
        )

        # We use scanpy's embedding plot which handles a list of colors automatically
        # by creating a grid of subplots.
        sc.pl.embedding(
            adata,
            basis=embedding_key,
            color=umap_color_columns,
            show=False,
            frameon=False,
            wspace=0.3,
        )

        umap_out = module_dir / f"embedding_{embedding_key}.png"
        plt.savefig(umap_out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved embedding plot to {umap_out}")
    else:
        logger.warning(
            f"Embedding '{embedding_key}' not found in adata.obsm. Skipping UMAP plot."
        )

    # View plots
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

    # 2. View plots (Leiden Clusters with FOV Bounding Boxes)
    fov_col = "fov" if "fov" in adata.obs.columns else None

    if fov_col:
        logger.info("Plotting clusters with FOV bounding boxes...")

        # LOOP OVER SAMPLES TO DRAW RECTANGLES PROPERLY
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

    if data_type == "Xenium":
        logger.info("Generating spatial ROI grid for Xenium data...")

        for sample in adata.obs[sample_key].unique():
            adata_sample = adata[adata.obs[sample_key] == sample].copy()
            coords = adata_sample.obsm[spatial_key]

            x_min_global, x_max_global = coords[:, 0].min(), coords[:, 0].max()
            y_min_global, y_max_global = coords[:, 1].min(), coords[:, 1].max()

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
            df_grid = pd.DataFrame(box_data)
            df_grid = df_grid.round(2)
            grid_csv_path = module_dir / f"Xenium_ROI_grid_coordinates_{sample}.csv"
            df_grid.to_csv(grid_csv_path, index=False)
            logger.info(f"Saved Xenium ROI grid coordinates to {grid_csv_path}")

    # 3. View specific gene expression
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
