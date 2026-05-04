"""Spatial statistics module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
from config import settings, get_module

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_spatial_statistics(module_dir, prev_module_dir, cluster_name):
    """Run spatial statistics."""

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    sc.settings.figdir = module_dir  # set the figures dir to not be figures

    # Import data
    logger.info("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Make sure cluster column is categorical and get the number of unique clusters
    if adata.obs[cluster_name].dtype.name != "category":
        adata.obs[cluster_name] = adata.obs[cluster_name].astype("category")

    valid_clusters = list(adata.obs[cluster_name].cat.categories)
    num_clusters = len(valid_clusters)
    logger.info(f"Detected {num_clusters} unique clusters. Adjusting plot sizes...")

    # --- DYNAMIC FIGURE SIZE CALCULATIONS ---
    # Centrality: 3 panels horizontally. Scales width heavily, height slightly.
    cent_width = max(16.0, num_clusters * 1.5)
    cent_figsize = (cent_width, max(5.0, num_clusters * 0.3))

    # Co-occurrence: Grid plot or plot with big legend. Scales width and height.
    co_size = max(10.0, num_clusters * 4)
    co_figsize = (co_size, 10.0)

    # Neighborhood enrichment: Heatmap. Scales width and height equally.
    nhood_size = max(8.0, num_clusters * 0.6)
    nhood_figsize = (nhood_size, nhood_size)
    # ----------------------------------------

    # Calculate spatial statistics
    logger.info("Building spatial neighborhood graph...")
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", delaunay=True
    )  # compute connectivity

    logger.info("Computing and plotting centrality scores...")
    sq.gr.centrality_scores(adata, cluster_key=cluster_name)
    sq.pl.centrality_scores(
        adata,
        cluster_key=cluster_name,
        figsize=cent_figsize,
        save="centrality_scores.png",
    )
    logger.info(
        f"Centrality scores plot saved to {module_dir / 'centrality_scores.png'}"
    )

    # Compute co-occurrence probability
    logger.info("Computing co-occurrence probability...")
    # Create subset table layer
    adata_subsample = sc.pp.subsample(
        adata, fraction=0.5, copy=True
    )  # subsample to speed up computation

    # Visualize co-occurrence
    sq.gr.co_occurrence(
        adata_subsample,
        cluster_key=cluster_name,
    )

    # Ensure subsampled data has exact same valid categories for plotting
    valid_clusters_sub = list(adata_subsample.obs[cluster_name].cat.categories)

    sq.pl.co_occurrence(
        adata_subsample,
        cluster_key=cluster_name,
        clusters=valid_clusters_sub,
        figsize=co_figsize,
        save="co_occurrence.png",
    )
    logger.info(f"Co-occurrence plot saved to {module_dir / 'co_occurrence.png'}")

    # Neighborhood enrichment analysis
    logger.info("Performing neighborhood enrichment analysis...")
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_name)

    # Plot neighborhood enrichment
    sq.pl.nhood_enrichment(
        adata,
        cluster_key=cluster_name,
        figsize=nhood_figsize,
        title="Neighborhood enrichment",
        save="nhood_enrichment.png",
    )
    logger.info(
        f"Neighborhood enrichment plot saved to {module_dir / 'nhood_enrichment.png'}"
    )

    # Moran's I
    logger.info("Calculating Moran's I...")

    # Build spatial neighborhood graph on a subsample dataset
    sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

    # Calculate Moran's I for spatial autocorrelation on subsample data
    sq.gr.spatial_autocorr(
        adata_subsample,
        mode="moran",
        n_perms=100,
        n_jobs=1,
    )

    # Save Moran's I results
    adata_subsample.uns["moranI"].to_csv(module_dir / "moranI_results.csv", index=True)
    logger.info(f"Moran's I results saved to {module_dir / 'moranI_results.csv'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Spatial statistics module completed successfully.")


if __name__ == "__main__":
    module_4_name, module_4_dir = get_module(4)
    module_5_name, module_5_dir = get_module(5)

    cluster_labels = settings["modules"]["Squidpy"]["cluster_labels"]

    run_spatial_statistics(module_5_dir, module_4_dir, cluster_labels)
