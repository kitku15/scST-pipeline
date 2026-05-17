"""MuSpan handling module."""

import warnings
import logging
from logging import getLogger
from pathlib import Path
import scanpy as sc

import spatialdata as sd
from MuSpan_CosMxHandling import CosMx_to_domain, cosmx_initial_plotting
from MuSpan_XeniumHandling import xenium_initial_plotting

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


def run_muspan(
    dataset_type: str,
    module_dir: str,
    prev_module_dir: str,
    domain_name: str,
    cluster_labels: str,  # chosen cluster
    transcripts_of_interest: list = None,
    cell_selection_csv: str = None,  # used for both cosmx or xenium
    # CosMx specific kwargs
    zarr_path: str = None,
    flat_files_dir: str = None,
    # Xenium specific kwargs
    xenium_dir: str = None,
):
    if dataset_type not in ["CosMx", "Xenium"]:
        print(dataset_type)
        raise ValueError("dataset_type must be either 'CosMx' or 'Xenium'")

    out_dir = Path(module_dir) / domain_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load AnnData
    adata_path = Path(prev_module_dir) / "adata.h5ad"
    logger.info(f"Loading Annotated AnnData from {adata_path}...")
    adata = sc.read_h5ad(adata_path)

    # Make CosMx domain, and a figure with 2 plots
    # Selected transcripts (left)
    # Cell boundaries colored by cluster labels (right)
    if dataset_type == "CosMx":
        logger.info("Processing CosMx data...")
        domain = ms.domain(domain_name)
        sdata = sd.read_zarr(zarr_path)

        # Create MuSpAn domain and add transcripts and cell boundaries based on cell selection
        # basically custom version of the xenium to domain function but for CosMx
        domain, qTrans, qCells, qBoundaries = CosMx_to_domain(
            sdata,
            domain,
            transcripts_of_interest,
            cluster_labels,
            cell_selection_csv,
            flat_files_dir,
            out_dir,
        )
        # Make cosmx figure
        logger.info(f"Visualizing the MuSpAn domain: {domain_name}")
        cosmx_initial_plotting(
            domain, qTrans, qCells, qBoundaries, cluster_labels, out_dir
        )

    # Make Xenium domain, and a figure
    elif dataset_type == "Xenium":
        logger.info("Processing Xenium data...")

        # Load Xenium Domain directly via MuSpan io function
        domain = ms.io.xenium_to_domain(
            path_to_xenium_data=str(xenium_dir),
            cells_from_selection_csv=str(cell_selection_csv),
            domain_name=domain_name,
            load_transcripts=True,
            selected_transcripts=transcripts_of_interest,
            load_nuclei=True,
            load_cells_as_shapes=True,
            exclude_no_nuclei_cells=True,
        )
        # Make xenium figure
        xenium_initial_plotting(adata, domain, cluster_labels, out_dir)

    # Save domain (do we need to this takes ages?)=========
    # ms.io.save_domain(
    #     domain,
    #     name_of_file=f"muspan_object_{domain_name.replace(' ', '_')}",
    #     path_to_save=str(out_dir),
    # )
    # logger.info(f"Domain saved successfully to {out_dir}")
    # =======================================================

    return domain
