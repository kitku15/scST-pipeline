"""Dimension reduction module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
from config import settings, get_module
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_dimension_reduction(
    data_type,
    prev_module_dir,
    module_dir,
    module_name,
    n_comps,
    n_neighbors_list,
    resolution_list,
    cluster_name,
):
    """Run dimension reduction on data."""

    # Ensure inputs are lists for iteration
    if not isinstance(n_neighbors_list, list):
        n_neighbors_list = [n_neighbors_list]
    if not isinstance(resolution_list, list):
        resolution_list = [resolution_list]

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    cluster_palette_25 = [
        "#e6194b",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#46f0f0",
        "#f032e6",
        "#bcf60c",
        "#fabebe",
        "#008080",
        "#e6beff",
        "#9a6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000075",
        "#808080",
        "#ff4500",
        "#2e8b57",
        "#1e90ff",
        "#ff1493",
        "#d1008f",
    ]

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir
    sc.set_figure_params(
        facecolor="white", transparent=False, dpi=300, figsize=(12, 12)
    )

    # Import data
    logger.info("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Perform dimension reduction analysis
    logger.info("Compute PCA...")
    sc.pp.pca(adata, n_comps=n_comps)  # compute principal components
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=50,
        show=False,
        save="PCA.png",
    )
    logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    for n_neighbors in n_neighbors_list:
        logger.info(f"Compute neighbors for n={n_neighbors}...")

        # 1. Create a unique key for this neighbor graph
        neighbors_key = f"neighbors_n{n_neighbors}"
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=30,
            key_added=neighbors_key,  # Save graph to unique key
        )

        logger.info(f"Create UMAPs and cluster cells for n={n_neighbors}...")

        # Tell UMAP to use the specific neighbors graph
        sc.tl.umap(adata, neighbors_key=neighbors_key)

        # Scanpy saves the UMAP to 'X_umap' by default.
        # We must copy it to a unique name so the next loop doesn't overwrite it
        custom_umap_basis = f"umap_n{n_neighbors}"
        adata.obsm[f"X_{custom_umap_basis}"] = adata.obsm["X_umap"].copy()

        for resolution in resolution_list:
            current_cluster_name = f"{cluster_name}_n{n_neighbors}_r{resolution}"

            # combination specific subfolder
            combo_dir = module_dir / f"n{n_neighbors}_r{resolution}"
            combo_dir.mkdir(parents=True, exist_ok=True)

            # Point Scanpy/Squidpy to save figures in this subfolder
            sc.settings.figdir = combo_dir

            # Tell Leiden to use the specific neighbors graph

            logger.info(f"Running Leiden clustering for {current_cluster_name}...")
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=current_cluster_name,
                neighbors_key=neighbors_key,  # use correct graph
            )

            n_clusters = adata.obs[current_cluster_name].nunique()

            # handle palettes when there are > 25 clusters
            if n_clusters <= len(cluster_palette_25):
                adata.uns[f"{current_cluster_name}_colors"] = cluster_palette_25[
                    :n_clusters
                ]
            else:
                # Loop the palette so Scanpy doesn't crash from missing colors
                repeated_palette = cluster_palette_25 * (
                    (n_clusters // len(cluster_palette_25)) + 1
                )
                adata.uns[f"{current_cluster_name}_colors"] = repeated_palette[
                    :n_clusters
                ]

            # plot UMAP
            logger.info(f"Plotting UMAPs for {current_cluster_name}...")
            sc.pl.embedding(
                adata,
                basis=custom_umap_basis,  # Tell plot to use our uniquely saved UMAP
                color=[
                    "total_counts",
                    "n_genes_by_counts",
                    current_cluster_name,
                ],
                wspace=0.4,
                show=False,
                save=f"_{current_cluster_name}.png",
                frameon=False,
            )

            logger.info(f"Plotting Spatial Scatter for {current_cluster_name}...")
            fig, ax = plt.subplots(figsize=(15, 15), facecolor="white")
            ax.set_facecolor("white")

            sq.pl.spatial_scatter(
                adata,
                color=[current_cluster_name],
                spatial_key=spatial_key,
                shape=None,
                facecolor="white",
                size=2,
                frameon=False,
                img=False,
                ax=ax,
                outline=False,
                dpi=300,
            )

            fig.savefig(
                combo_dir / f"{current_cluster_name}_spatial.png",
                dpi=300,
                facecolor="white",
                bbox_inches="tight",
            )
            plt.close(fig)

    # Reset global figdir
    sc.settings.figdir = module_dir

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")


if __name__ == "__main__":
    data_type = settings["project"]["data_type"]

    module_1_name, module_1_dir = get_module(1)
    module_2_name, module_2_dir = get_module(2)

    n_comps = settings["modules"]["DimensionReduction"]["n_comps"]
    n_neighbors = settings["modules"]["DimensionReduction"]["n_neighbors"]
    resolution = settings["modules"]["DimensionReduction"]["resolution"]
    cluster_name = settings["modules"]["DimensionReduction"]["cluster_name"]

    run_dimension_reduction(
        data_type,
        module_1_dir,
        module_2_dir,
        module_2_name,
        n_comps,
        n_neighbors,
        resolution,
        cluster_name,
    )
