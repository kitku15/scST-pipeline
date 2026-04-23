"""Spatial statistics module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
from config import get_module

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_spatial_statistics(module_dir, prev_module_dir):
    """Run spatial statistics."""

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    sc.settings.figdir = module_dir  # set the figures dir to not be figures

    # Import data
    logger.info("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Calculate spatial statistics
    logger.info("Building spatial neighborhood graph...")
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", delaunay=True
    )  # compute connectivity

    logger.info("Computing and plotting centrality scores...")
    sq.gr.centrality_scores(adata, cluster_key="leiden")
    sq.pl.centrality_scores(
        adata,
        cluster_key="leiden",
        figsize=(16, 5),
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
        cluster_key="leiden",
    )

    valid_clusters = list(adata_subsample.obs["leiden"].cat.categories)

    sq.pl.co_occurrence(
        adata_subsample,
        cluster_key="leiden",
        clusters=valid_clusters,
        figsize=(10, 10),
        save="co_occurrence.png",
    )
    logger.info(f"Co-occurrence plot saved to {module_dir / 'co_occurrence.png'}")

    # Neighborhood enrichment analysis
    logger.info("Performing neighborhood enrichment analysis...")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")

    # Plot neighborhood enrichment
    sq.pl.nhood_enrichment(
        adata,
        cluster_key="leiden",
        figsize=(8, 8),
        title="Neighborhood enrichment adata",
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

    run_spatial_statistics(module_5_dir, module_4_dir)
