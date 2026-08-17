"""MuSpan Xenium handling module."""

import logging
import warnings
from logging import getLogger

import matplotlib.pyplot as plt

try:
    import muspan as ms
except ModuleNotFoundError as err:
    raise ModuleNotFoundError(
        "Could not load necessary MuSpAn package. "
        "Please install it before running this script."
    ) from err

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)
logger = getLogger(__name__)


def map_cell_types_to_domain_xenium(adata, domain, adata_cell_id, cluster_labels):
    """Helper function to map AnnData cell types to a Muspan domain (Used for Xenium)"""
    domain_cell_ids_ordered = [
        str(cell_id) for cell_id in domain.labels["Cell ID"]["labels"]
    ]
    domain_cell_ids_unique = set(domain_cell_ids_ordered)

    logger.info(f"Number of unique cells in the domain: {len(domain_cell_ids_unique)}")

    filt_adata = adata[adata.obs[adata_cell_id].isin(domain_cell_ids_unique)]

    logger.info(f"Filtered adata from {adata.n_obs} to {filt_adata.n_obs} cells")
    logger.info("Adding cell_type IDs to domain with cluster labels")

    cell_id_to_type = dict(
        zip(filt_adata.obs[adata_cell_id], filt_adata.obs[cluster_labels])
    )
    cell_types_ordered = [
        cell_id_to_type.get(cell_id, "Unknown") for cell_id in domain_cell_ids_ordered
    ]

    domain.add_labels(label_name=cluster_labels, labels=cell_types_ordered)
    return domain


def xenium_initial_plotting(adata, domain, cluster_labels, out_dir, color_dict=None):
    # Map labels
    cellid_col = "cell_id"
    domain = map_cell_types_to_domain_xenium(adata, domain, cellid_col, cluster_labels)

    if color_dict is not None:
        try:
            domain.update_colors(
                color_dict, colors_to_update="labels", label_name=cluster_labels
            )
            logger.info("Successfully synced Squidpy colors to MuSpAn Xenium Domain.")
        except Exception as e:
            logger.warning(f"Failed to apply custom colors to MuSpAn domain: {e}")

    # Queries
    qCells = ms.query.query(domain, ("Collection",), "is", "Cell boundaries")
    qTrans = ms.query.query(domain, ("Collection",), "is", "Transcripts")
    qNuc = ms.query.query(domain, ("Collection",), "is", "Nucleus boundaries")

    # Plotting 1: Overviews (Xenium)
    fig, ax = plt.subplots(figsize=(20, 15), nrows=2, ncols=2)
    ms.visualise.visualise(domain, ax=ax[0, 0], marker_size=0.05)
    ax[0, 0].set_title("All objects")

    ms.visualise.visualise(
        domain,
        color_by=("label", cluster_labels),
        ax=ax[0, 1],
        objects_to_plot=qCells,
    )
    ax[0, 1].set_title("Cell type")

    ms.visualise.visualise(
        domain,
        color_by=("label", "Transcript ID"),
        ax=ax[1, 0],
        objects_to_plot=qTrans,
        marker_size=5,
    )
    ax[1, 0].set_title("Transcripts")

    ms.visualise.visualise(
        domain,
        color_by=("label", "Nucleus Area"),
        ax=ax[1, 1],
        objects_to_plot=qNuc,
    )
    ax[1, 1].set_title("Nuclei")

    plt.tight_layout()
    plt.savefig(out_dir / "muspan_domain_visualization.png")

    # Convert boundaries to centroids
    logger.info("Convert cell boundaries to cell centres (centroids)")
    domain.convert_objects(
        population=("Collection", "Cell boundaries"),
        object_type="point",
        conversion_method="centroids",
        collection_name="Cell centroids",
        inherit_collections=False,
    )

    # Plotting 2: Centroids (Xenium)
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by=cluster_labels,
        ax=plt.gca(),
    )
    plt.tight_layout()
    plt.savefig(out_dir / "muspan_cell_centroids.png")

    # Plotting 3: Overlay (Xenium)
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain, objects_to_plot=qCells, color_by=cluster_labels, ax=plt.gca()
    )
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by=("constant", "black"),
        ax=plt.gca(),
        marker_size=1,
        add_cbar=False,
    )
    plt.tight_layout()
    plt.savefig(out_dir / "muspan_cell_centroids_n_boundaries.png")
    plt.close("all")
