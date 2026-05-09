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
import squidpy as sq
import matplotlib

matplotlib.use("Agg")


warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_qc(data_type, module_dir, zarr_path, min_counts, min_cells, min_dapi):
    # 1. Platform-specific configuration
    if data_type == "CosMx":
        cfg = {
            "area_col": "Area",
            "has_dapi": True,
            "control_pattern": "^NegPrb|^SystemControl",
            "spatial_key": "global",
            "nucleus_col": "Mean.DAPI",
        }
    elif data_type == "Xenium":
        cfg = {
            "area_col": "cell_area",
            "has_dapi": False,
            "control_pattern": "control_probe|control_codeword",
            "spatial_key": "spatial",
            "nucleus_col": "nucleus_area",
        }

    # 2. Load data
    sdata = sd.read_zarr(zarr_path)
    adata = sdata.tables["table"]
    del sdata
    gc.collect()

    # print(adata.obs.columns)
    # print(adata.var.columns)

    # 3. Harmonize CosMx to Xenium logic
    if data_type == "CosMx":
        spatial_key = cfg["spatial_key"]

        # Copy to avoid any view/SettingWithCopy warnings
        coords = adata.obsm[spatial_key].copy()

        # Invert the Y-axis (column index 1) while keeping values positive
        coords[:, 1] = np.max(coords[:, 1]) - coords[:, 1]

        # Reassign the updated coordinates back to the AnnData object
        adata.obsm[spatial_key] = coords
        logger.info(f"Permanently flipped Y-axis coordinates for {data_type}.")

        # Create control counts manually for CosMx so logging/metrics look the same as Xenium
        neg_probes = adata.var_names.str.contains("^NegPrb", case=False)
        sys_controls = adata.var_names.str.contains("^SystemControl", case=False)

        print(
            f"CosMx: Found {neg_probes.sum()} negative probes and {sys_controls.sum()} system controls in the data."
        )

        adata.obs["control_probe_counts"] = np.array(
            adata[:, neg_probes].X.sum(axis=1)
        ).flatten()
        adata.obs["control_codeword_counts"] = np.array(
            adata[:, sys_controls].X.sum(axis=1)
        ).flatten()
        # Proxy for nucleus plot (CosMx doesn't have a nucleus_area col usually)
        adata.obs["nucleus_area"] = adata.obs["Area"]

    # 4. Calculate QC metrics
    # Tag control genes in var for Scanpy metrics
    adata.var["control"] = adata.var_names.str.contains(
        cfg["control_pattern"], case=False, na=False
    )

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["control"],
        percent_top=(10, 20, 50, 150),
        inplace=True,
    )

    # 5. Logging (Logic is now identical for both)
    cprobes = (
        adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
    )
    cwords = (
        adata.obs["control_codeword_counts"].sum()
        / adata.obs["total_counts"].sum()
        * 100
    )

    logger.info(f"Negative DNA probe count % : {cprobes:.4f}%")
    logger.info(f"Negative decoding count % : {cwords:.4f}%")

    avg_total_counts = np.mean(adata.obs["total_counts"])
    logger.info(f"Average number of transcripts per cell: {avg_total_counts:.2f}")

    avg_total_unique_counts = np.mean(adata.obs["n_genes_by_counts"])
    logger.info(f"Average unique transcripts per cell: {avg_total_unique_counts:.2f}")

    area_max = np.max(adata.obs[cfg["area_col"]])
    area_min = np.min(adata.obs[cfg["area_col"]])
    logger.info(f"Max cell area: {area_max}")
    logger.info(f"Min cell area: {area_min}")

    # Plotting
    plot_metrics(module_dir, adata, cfg, min_counts, min_dapi)
    plot_spatial_qc(module_dir, adata, cfg)

    # 6. Filtering and Normalization
    logger.info("Filtering cells and genes...")

    if cfg["has_dapi"]:
        logger.info("Applying DAPI filter...")
        adata = adata[adata.obs["Mean.DAPI"] > min_dapi].copy()

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    logger.info("Normalize data...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata

    sc.pp.scale(adata, max_value=10)

    # Save
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Quality control completed successfully.")


def plot_metrics(module_dir, adata, cfg, min_counts, min_dapi):
    module_dir.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    # 1. Total transcripts
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], color="blue")
    axs[0].axvline(min_counts, color="red", linestyle="--", linewidth=2)

    # 2. Unique transcripts
    axs[1].set_title("Unique genes per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], color="green")

    # 3. Cell Area
    axs[2].set_title("Cell Area (Total)")
    sns.histplot(adata.obs[cfg["area_col"]], ax=axs[2], color="orange")

    # 4. Nucleus / DAPI plot
    if cfg["has_dapi"]:
        axs[3].set_title("Mean DAPI (Nucleus Signal)")
        sns.histplot(adata.obs["Mean.DAPI"], ax=axs[3], color="purple")
        axs[3].axvline(min_dapi, color="red", linestyle="--", linewidth=2)
    else:
        axs[3].set_title("Nucleus ratio")
        # Ensure division by zero doesn't happen if area is missing
        sns.histplot(adata.obs["nucleus_area"] / adata.obs[cfg["area_col"]], ax=axs[3])

    plt.tight_layout()
    out_file = module_dir / "cell_summary_histograms.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved plots to {out_file.absolute()}")


def plot_spatial_qc(module_dir, adata, cfg):
    module_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = module_dir

    logger.info(f"Visualize {str(cfg['nucleus_col'])} on tissue...")

    fig, ax = plt.subplots(figsize=(6, 6), facecolor="white")
    ax.set_facecolor("white")

    sq.pl.spatial_scatter(
        adata,
        spatial_key=cfg["spatial_key"],
        color=cfg["nucleus_col"],
        shape=None,
        outline=False,
        wspace=0.4,
        size=1,
        dpi=300,
        img=False,
        ax=ax,
    )

    fig.savefig(
        module_dir / f"{cfg['nucleus_col']}_scatter.png",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)

    logger.info("Visualize cell area on tissue...")
    fig, ax = plt.subplots(figsize=(6, 6), facecolor="white")
    ax.set_facecolor("white")

    sq.pl.spatial_scatter(
        adata,
        spatial_key=cfg["spatial_key"],
        color=cfg["area_col"],
        shape=None,
        outline=False,
        wspace=0.4,
        size=1,
        dpi=300,
        img=False,
        ax=ax,
    )

    fig.savefig(
        module_dir / "Area_scatter.png", dpi=300, facecolor="white", bbox_inches="tight"
    )
    plt.close(fig)


if __name__ == "__main__":
    module_1_name, module_1_dir = get_module(1)
    zarr_path = settings["io"]["zarr_dir"]
    data_type = settings["project"]["data_type"]

    min_counts = settings["modules"]["QualityControl"]["min_counts"]
    min_cells = settings["modules"]["QualityControl"]["min_cells"]
    min_dapi = settings["modules"]["QualityControl"]["min_dapi"]

    run_qc(data_type, module_1_dir, zarr_path, min_counts, min_cells, min_dapi)
