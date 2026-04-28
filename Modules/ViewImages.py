"""Image viewing module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
from config import settings, get_module

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_view_images(data_type, prev_module_dir, module_dir, gene_list, cluster_name):
    """Run the image viewing module."""

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    sc.settings.figdir = module_dir  # set the figures dir to not be figures

    # Import data
    logger.info("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # View plots
    logger.info("Visualize clusters on tissue...")
    sq.pl.spatial_scatter(
        adata,
        spatial_key=spatial_key,
        shape=None,
        outline=False,
        color=[cluster_name, "total_counts"],
        wspace=0.4,
        size=1,
        save="leiden_clusters.png",
        dpi=300,
    )
    logger.info(f"Saved leiden clusters plot to {module_dir / 'leiden_clusters.png'}")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        spatial_key=spatial_key,
        color=gene_list,
        shape=None,
        size=2,
        img=False,
        save="gene_expression.png",
    )
    logger.info(f"Saved gene expression plot to {module_dir / 'gene_expression.png'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Imaging module completed successfully.")


if __name__ == "__main__":
    data_type = settings["project"]["data_type"]
    module_3_name, module_3_dir = get_module(3)
    module_4_name, module_4_dir = get_module(4)
    gene_list = settings["modules"]["ViewImages"]["gene_list"]

    run_view_images(data_type, module_3_dir, module_4_dir, gene_list)
