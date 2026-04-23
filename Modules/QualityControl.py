"""Quality control module."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import spatialdata as sd
from config import settings, get_module
import gc

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_qc(data_type, module_dir, zarr_path, min_counts, min_cells, min_dapi):
    if data_type == "CosMx":
        area_col = "Area"
        DAPI_filter = True
    elif data_type == "Xenium":
        area_col = "cell_area"
        DAPI_filter = False

    # spatialdata object
    sdata = sd.read_zarr(zarr_path)

    # Save anndata object (stored in spatialdata.tables layer)
    adata = sdata.tables["table"]

    del sdata  # Free up memory by deleting the spatialdata object
    gc.collect()

    logger.info("logger.infoing adata obs collumns-----")
    logger.info(adata.obs.columns)

    # exit()

    # 1. Identify the control probes and codewords in the gene list
    # CosMx typically labels these with "Negative" or "SystemControl"
    adata.var["is_control_probe"] = adata.var_names.str.contains(
        "Negative", case=False, na=False
    )
    adata.var["is_control_codeword"] = adata.var_names.str.contains(
        "SystemControl", case=False, na=False
    )

    # 2. Calculate QC metrics, instructing scanpy to isolate our controls
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[
            "is_control_probe",
            "is_control_codeword",
        ],  # Tells scanpy to sum these up
        percent_top=(10, 20, 50, 150),
        inplace=True,
    )
    # logger.info(adata.obs.columns)

    # 3. Use the new columns scanpy just generated for us
    cprobes = (
        adata.obs["total_counts_is_control_probe"].sum()
        / adata.obs["total_counts"].sum()
        * 100
    )
    cwords = (
        adata.obs["total_counts_is_control_codeword"].sum()
        / adata.obs["total_counts"].sum()
        * 100
    )

    logger.info(f"Negative DNA probe count % : {cprobes:.4f}%")
    logger.info(f"Negative decoding count % : {cwords:.4f}%")

    # Calculate averages
    avg_total_counts = np.mean(adata.obs["total_counts"])
    logger.info(f"Average number of transcripts per cell: {avg_total_counts:.2f}")

    avg_total_unique_counts = np.mean(adata.obs["n_genes_by_counts"])
    logger.info(f"Average unique transcripts per cell: {avg_total_unique_counts:.2f}")

    area_max = np.max(adata.obs[area_col])
    area_min = np.min(adata.obs[area_col])

    logger.info(f"Max cell area: {area_max}")
    logger.info(f"Min cell area: {area_min}")

    # plot raw data
    plot_metrics(module_dir, adata, area_col, DAPI_filter)

    # $ QC data #

    # Filter cells
    logger.info("Filtering cells and genes...")

    if DAPI_filter:
        logger.info("Applying DAPI filter...")
        # Filter out the 'empty' cells seen in your DAPI histogram
        adata = adata[adata.obs["Mean.DAPI"] > min_dapi].copy()

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize data
    logger.info("Normalize data...")
    adata.layers["counts"] = adata.X.copy()  # make copy of raw data
    sc.pp.normalize_total(adata, inplace=True)  # normalize data
    sc.pp.log1p(adata)  # Log transform data

    # Identify the genes that actually matter for distinguishing cell types
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # Save the raw state for plotting later
    adata.raw = adata

    # Scale the data so highly expressed genes don't overpower the PCA
    sc.pp.scale(
        adata, max_value=10
    )  # Sparse to Dense Matrix Conversion (storing every single zero as a physical number in RAM)

    # Save data
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Quality control completed successfully.")


def plot_metrics(module_dir, adata, area_col, DAPI_filter):
    module_dir.mkdir(parents=True, exist_ok=True)

    # We will plot 4 metrics, using DAPI as the 4th
    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    # 1. Total transcripts
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], color="blue")

    # 2. Unique transcripts
    axs[1].set_title("Unique genes per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], color="green")

    # 3. Cell Area
    axs[2].set_title("Cell Area (Total)")
    sns.histplot(adata.obs[area_col], kde=False, ax=axs[2], color="orange")

    # 4. Nucleus / DAPI plot
    if DAPI_filter:
        axs[3].set_title("Mean DAPI (Nucleus Signal)")
        sns.histplot(adata.obs["Mean.DAPI"], kde=False, ax=axs[3], color="purple")
    else:
        axs[3].set_title("Nucleus ratio")
        sns.histplot(
            adata.obs["nucleus_area"] / adata.obs["cell_area"],
            kde=False,
            ax=axs[3],
        )

    plt.tight_layout()
    out_file = module_dir / "cell_summary_histograms.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved plots to {out_file.absolute()}")


# Use pathlib.Path so the '/' operator in the function works correctly
if __name__ == "__main__":
    module_1_name, module_1_dir = get_module(1)
    zarr_path = settings["io"]["zarr_dir"]
    data_type = settings["project"]["data_type"]

    min_counts = settings["modules"]["QualityControl"]["min_counts"]
    min_cells = settings["modules"]["QualityControl"]["min_cells"]
    min_dapi = settings["modules"]["QualityControl"]["min_dapi"]

    run_qc(data_type, module_1_dir, zarr_path, min_counts, min_cells, min_dapi)
