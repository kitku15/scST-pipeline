"""Image viewing module."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import scanpy as sc
import squidpy as sq
import math

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_view_images(
    data_type,
    prev_module_dir,
    module_dir,
    gene_list,
    cluster_name,
    n_grid_x=10,
    n_grid_y=10,
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
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # View plots
    logger.info("Visualize clusters on tissue...")
    fig, axes = plt.subplots(1, 2, figsize=(20, 6), facecolor="white")

    sq.pl.spatial_scatter(
        adata,
        spatial_key=spatial_key,
        shape=None,
        outline=False,
        color=[cluster_name, "total_counts"],
        wspace=0.4,
        size=1,
        dpi=300,
        ax=axes,
        fig=fig,
    )
    fig.savefig(
        module_dir / "leiden_clusters.png",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)
    logger.info(f"Saved leiden clusters plot to {module_dir / 'leiden_clusters.png'}")

    # 2. View plots (Leiden Clusters with FOV Bounding Boxes)

    # Determine the column name for FOV (usually 'fov' or 'FOV')
    fov_col = (
        "fov" if "fov" in adata.obs.columns else None
    )  # this will be None for Xenium

    if fov_col:
        logger.info("Plotting clusters with FOV bounding boxes...")
        fig_fov, ax_fov = plt.subplots(figsize=(10, 10), facecolor="white")

        sq.pl.spatial_scatter(
            adata,
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
        coords = adata.obsm[spatial_key]
        df_coords = pd.DataFrame(
            {"x": coords[:, 0], "y": coords[:, 1], "fov": adata.obs[fov_col].values}
        )

        # Calculate min/max for each FOV and draw boxes
        for fov_id, group in df_coords.groupby("fov"):
            x_min, x_max = group["x"].min(), group["x"].max()
            y_min, y_max = group["y"].min(), group["y"].max()

            width = x_max - x_min
            height = y_max - y_min

            # Create a Rectangle patch (Updated to match Xenium style)
            rect = patches.Rectangle(
                (x_min, y_min),
                width,
                height,
                linewidth=1.2,  # Changed from 0.8 to 1.2
                edgecolor="black",
                facecolor="none",
                linestyle="solid",  # Changed from "dashed" to "solid"
            )
            ax_fov.add_patch(rect)

            # Add FOV label in the center of the box (Updated to match Xenium style)
            ax_fov.text(
                x_min + (width / 2),
                y_min + (height / 2),
                str(fov_id),
                color="black",
                fontsize=12,  # Changed from 10 to 12
                fontweight="bold",  # Added bold text
                ha="center",
                va="center",
                # Changed boxstyle to circle and added a black edgecolor
                bbox=dict(
                    boxstyle="circle,pad=0.3",
                    facecolor="white",
                    edgecolor="black",
                    alpha=0.8,
                ),
            )

        fig_fov.savefig(
            module_dir / "FOV_mapping.png",
            dpi=300,
            facecolor="white",
            bbox_inches="tight",
        )
        plt.close(fig_fov)
        logger.info(f"Saved FOV clusters plot to {module_dir / 'FOV_mapping.png'}")
    else:
        logger.warning(
            "No 'fov' or 'FOV' column found in adata.obs. Skipping FOV bounding boxes plot."
        )

    if data_type == "Xenium":
        # script that divides into specific regions of interest
        # based on spatial coordinates (e.g. bounding boxes)
        # just a nice grid, numbered 1 to n for each box,
        # also saves csv box number and coordinates of bounding boxes for later use in selection csv module
        logger.info("Generating spatial ROI grid for Xenium data...")

        coords = adata.obsm[spatial_key]
        x_min_global, x_max_global = coords[:, 0].min(), coords[:, 0].max()
        y_min_global, y_max_global = coords[:, 1].min(), coords[:, 1].max()

        x_step = (x_max_global - x_min_global) / n_grid_x
        y_step = (y_max_global - y_min_global) / n_grid_y

        box_data = []
        box_id = 1

        # Build grid coordinates (top-to-bottom, left-to-right)
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

        # Plot the grid on top of spatial data
        fig_grid, ax_grid = plt.subplots(figsize=(10, 10), facecolor="white")

        sq.pl.spatial_scatter(
            adata,
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

            # Center text for the Box ID
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

        grid_plot_path = module_dir / "Xenium_ROI_grid.png"
        fig_grid.savefig(
            grid_plot_path, dpi=300, facecolor="white", bbox_inches="tight"
        )
        plt.close(fig_grid)
        logger.info(f"Saved Xenium ROI grid plot to {grid_plot_path}")

        # Save coordinates to CSV
        df_grid = pd.DataFrame(box_data)
        df_grid = df_grid.round(2)  # round to 2 decimals
        grid_csv_path = module_dir / "Xenium_ROI_grid_coordinates.csv"
        df_grid.to_csv(grid_csv_path, index=False)
        logger.info(f"Saved Xenium ROI grid coordinates to {grid_csv_path}")

    # 3. View specific gene expression
    logger.info("Plotting genes of interest on tissue...")

    # dynamic number of cols and rows that adjust to no. of genes
    n_genes = len(gene_list)

    ncols = math.ceil(math.sqrt(n_genes))
    nrows = math.ceil(n_genes / ncols)

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 5 * nrows), facecolor="white"
    )

    if n_genes == 1:  # Handle the case where there's only 1 gene
        axes = [axes]
    else:
        axes = axes.flatten()

    for ax in axes[n_genes:]:
        ax.remove()

    axes = axes[:n_genes]

    sq.pl.spatial_scatter(
        adata,
        spatial_key=spatial_key,
        color=gene_list,
        shape=None,
        size=2,
        img=False,
        ax=axes,
        fig=fig,
    )
    fig.savefig(
        module_dir / "gene_expression.png",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)
    logger.info(f"Saved gene expression plot to {module_dir / 'gene_expression.png'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Imaging module completed successfully.")

    return grid_csv_path
