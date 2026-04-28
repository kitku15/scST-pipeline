"""Dimension reduction module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
from config import settings, get_module
import matplotlib

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
        save=f"_{module_name}.png",
    )
    logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    plotted_spatial_resolutions = set()

    for n_neighbors in n_neighbors_list:
        logger.info("Compute neighbors...")
        sc.pp.neighbors(
            adata, n_neighbors=n_neighbors, n_pcs=15
        )  # compute a neighborhood graph

        logger.info("Create UMAPs and cluster cells..")
        sc.tl.umap(adata)  # calculate umap
        for resolution in resolution_list:
            # Create a unique name for this specific parameter combination
            current_cluster_name = f"{cluster_name}_n{n_neighbors}_r{resolution}"
            sc.tl.leiden(
                adata,
                resolution=resolution,  # choose resolution for clustering
                key_added=current_cluster_name,
            )  # name leiden clusters

            # get number of clusters actually present
            n_clusters = adata.obs[current_cluster_name].nunique()

            # take only as many colors as needed
            adata.uns[f"{current_cluster_name}_colors"] = cluster_palette_25[
                :n_clusters
            ]

            # plot UMAP
            logger.info(f"Plotting UMAPs for {current_cluster_name}...")
            sc.pl.umap(
                adata,
                color=[
                    "total_counts",
                    "n_genes_by_counts",
                    current_cluster_name,
                ],
                wspace=0.4,
                show=False,
                save=f"_{current_cluster_name}_{module_name}.png",  # save the figure with the module name
                frameon=False,
            )
            logger.info(f"UMAP plot saved to {sc.settings.figdir}")

            if resolution not in plotted_spatial_resolutions:
                # plot visualization of leiden clusters
                logger.info(f"Plotting {current_cluster_name} clusters...")
                # Create a plot where each FOV is its own panel

                sq.pl.spatial_scatter(
                    adata,
                    color=[current_cluster_name],
                    spatial_key=spatial_key,
                    shape=None,
                    facecolor="white",
                    size=2,  # Use a very small size for the full slide
                    # alpha=0.6,
                    frameon=False,
                    img=False,
                    outline=False,
                    figsize=(15, 15),
                    save=f"{resolution}_full_stitched.png",
                    dpi=300,
                )

                plotted_spatial_resolutions.add(resolution)

            # if data_type == "CosMx":
            #     sq.pl.spatial_scatter(
            #         adata,
            #         color=[current_cluster_name],
            #         library_key="fov",  # Use the 'fov' column from your adata.obs
            #         ncols=4,  # Arrange in 4 columns
            #         shape=None,  # Circles
            #         size=1,  # Adjust size if dots are too big/small
            #         img=False,  # Keep False until we confirm coordinates are right
            #         outline=False,
            #         save=f"{current_cluster_name}_by_fov.png",
            #         dpi=300
            #     )

            # logger.info(f"{current_cluster_name} spatial scatter plot saved to {module_dir}")

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
