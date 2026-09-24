import gc
import json
import logging
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import spatialdata as sd
import squidpy as sq

matplotlib.use("Agg")
warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# ==========================================
# 1. Core Logic Functions
# ==========================================


def calculate_qc_metrics(adata, cfg):
    """Calculates probe counts and base QC metrics safely without dense RAM bloat."""
    neg_probes = adata.var_names.str.contains(cfg["control_pattern"], case=False)

    if sp.issparse(adata.X):
        adata.obs["control_counts"] = adata[:, neg_probes].X.sum(axis=1).A1
    else:
        adata.obs["control_counts"] = np.asarray(
            adata[:, neg_probes].X.sum(axis=1)
        ).squeeze()

    adata.var["control"] = neg_probes

    # 1. MUST RUN THIS FIRST: Generates 'total_counts' and 'n_genes_by_counts'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["control"], percent_top=(10, 20, 50, 150), inplace=True
    )

    # 2. NOW WE CAN LOG (since total_counts exists)
    c_percent = (
        adata.obs["control_counts"].sum() / adata.obs["total_counts"].sum()
    ) * 100
    logger.info(f"Negative/Control probe count % : {c_percent:.4f}%")
    logger.info(
        f"Average number of transcripts per cell: {np.mean(adata.obs['total_counts']):.2f}"
    )
    logger.info(
        f"Average unique transcripts per cell: {np.mean(adata.obs['n_genes_by_counts']):.2f}"
    )
    logger.info(
        f"Max cell area: {np.max(adata.obs[cfg['area_col']])} | Min cell area: {np.min(adata.obs[cfg['area_col']])}"
    )

    # 3. Handle Nucleus logic
    if not cfg["has_dapi"]:
        if "nucleus_area" not in adata.obs.columns:
            adata.obs["nucleus_area"] = adata.obs[cfg["area_col"]]
    elif "Mean.DAPI" in adata.obs.columns:
        adata.obs["nucleus_area"] = adata.obs["Mean.DAPI"]

    return adata


def filter_and_normalize(adata, qc_config, cfg):
    """Applies thresholds and strictly preserves integer counts in .layers['counts']."""
    start_cells = adata.n_obs
    logger.info("--- Starting Cell Filtering ---")

    if (
        cfg["has_dapi"]
        and qc_config.min_dapi is not None
        and "Mean.DAPI" in adata.obs.columns
    ):
        pre = adata.n_obs
        adata = adata[adata.obs["Mean.DAPI"] >= qc_config.min_dapi].copy()
        logger.info(
            f"DAPI filter (>= {qc_config.min_dapi}): removed {pre - adata.n_obs} cells."
        )

    if qc_config.min_area is not None and qc_config.max_area is not None:
        pre = adata.n_obs
        adata = adata[
            (adata.obs[cfg["area_col"]] >= qc_config.min_area)
            & (adata.obs[cfg["area_col"]] <= qc_config.max_area)
        ].copy()
        logger.info(
            f"Area filter ({qc_config.min_area} - {qc_config.max_area}): removed {pre - adata.n_obs} cells."
        )

    if qc_config.min_counts is not None:
        pre = adata.n_obs
        adata = adata[adata.obs["total_counts"] >= qc_config.min_counts].copy()
        logger.info(
            f"Min transcripts filter (>= {qc_config.min_counts}): removed {pre - adata.n_obs} cells."
        )

    if qc_config.min_genes is not None:
        pre = adata.n_obs
        adata = adata[adata.obs["n_genes_by_counts"] >= qc_config.min_genes].copy()
        logger.info(
            f"Min unique genes filter (>= {qc_config.min_genes}): removed {pre - adata.n_obs} cells."
        )

    logger.info(
        f"Retained {(adata.n_obs / start_cells) * 100:.1f}% of cells ({adata.n_obs})."
    )

    if qc_config.min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=qc_config.min_cells)

    logger.info("Normalizing data...")
    adata.layers["counts"] = adata.X.copy()  # Safe harbor for PyDESeq2 integer counts

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, base=2)

    # RESTORED: .raw assignment moved to AFTER normalization so LFC math doesn't break
    adata.raw = adata

    adata.obs.index.name = None
    return adata


# ==========================================
# 2. Output & Plotting Functions
# ==========================================


def export_web_metrics(module_dir, adata, cfg, qc_config):
    """Exports raw metrics to JSON/CSV for Spatial-VisKit."""
    logger.info("Exporting raw QC metrics to CSV and JSON for web visualization...")
    try:
        thresholds = {
            "min_counts": qc_config.min_counts,
            "min_genes": qc_config.min_genes,
            "min_cells": qc_config.min_cells,
            "min_dapi": qc_config.min_dapi if cfg["has_dapi"] else None,
            "min_area": qc_config.min_area,
            "max_area": qc_config.max_area,
            "has_dapi": cfg["has_dapi"],
            "area_col": cfg["area_col"],
            "nucleus_col": cfg["nucleus_col"],
        }
        with open(module_dir / "qc_thresholds.json", "w") as f:
            json.dump(thresholds, f)

        df_qc = pd.DataFrame(
            {
                "total_counts": adata.obs["total_counts"],
                "n_genes_by_counts": adata.obs["n_genes_by_counts"],
                "area": adata.obs[cfg["area_col"]],
            }
        )

        if cfg["has_dapi"] and "Mean.DAPI" in adata.obs.columns:
            df_qc["nucleus_signal"] = adata.obs["Mean.DAPI"]
        else:
            df_qc["nucleus_signal"] = adata.obs["nucleus_area"] / adata.obs[
                cfg["area_col"]
            ].replace(0, np.nan)

        df_qc.to_csv(module_dir / "qc_metrics.csv", index=False)
        logger.info("Successfully exported QC metrics.")
    except Exception as e:
        logger.warning(f"Failed to export QC metrics for web: {e}")


def plot_metrics(module_dir, adata, cfg, qc_config):
    """Generates and saves summary histograms of key cell-level QC metrics."""
    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    def get_upper_limit(series):
        return max(np.percentile(series.dropna(), 99.5), 1)

    # 1. Total transcripts
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], color="blue")
    if qc_config.min_counts:
        axs[0].axvline(qc_config.min_counts, color="red", linestyle="--", linewidth=2)
    axs[0].set_xlim(0, get_upper_limit(adata.obs["total_counts"]))

    # 2. Unique transcripts
    axs[1].set_title("Unique genes per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], color="green")
    if qc_config.min_genes:
        axs[1].axvline(qc_config.min_genes, color="red", linestyle="--", linewidth=2)
    axs[1].set_xlim(0, get_upper_limit(adata.obs["n_genes_by_counts"]))

    # 3. Cell Area
    axs[2].set_title("Cell Area (Total)")
    sns.histplot(adata.obs[cfg["area_col"]], ax=axs[2], color="orange")
    if qc_config.min_area:
        axs[2].axvline(qc_config.min_area, color="red", linestyle="--", linewidth=2)
    if qc_config.max_area:
        axs[2].axvline(qc_config.max_area, color="red", linestyle="--", linewidth=2)
    axs[2].set_xlim(0, get_upper_limit(adata.obs[cfg["area_col"]]))

    # 4. Nucleus / DAPI plot
    if cfg["has_dapi"] and "Mean.DAPI" in adata.obs.columns:
        axs[3].set_title("Mean DAPI (Nucleus Signal)")
        sns.histplot(adata.obs["Mean.DAPI"], ax=axs[3], color="purple")
        if qc_config.min_dapi:
            axs[3].axvline(qc_config.min_dapi, color="red", linestyle="--", linewidth=2)
        axs[3].set_xlim(0, get_upper_limit(adata.obs["Mean.DAPI"]))
    else:
        axs[3].set_title("Nucleus ratio")
        nuc_ratio = adata.obs["nucleus_area"] / adata.obs[cfg["area_col"]].replace(
            0, np.nan
        )
        sns.histplot(nuc_ratio, ax=axs[3])
        axs[3].set_xlim(0, get_upper_limit(nuc_ratio))

    plt.tight_layout()
    out_file = module_dir / "cell_summary_histograms.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved readable plots to {out_file.absolute()}")


def plot_spatial_qc(module_dir, adata, cfg):
    """Generates and saves spatial scatter plots of structural cell metrics."""
    sc.settings.figdir = module_dir

    logger.info(f"Visualize {cfg['nucleus_col']} on tissue...")
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


# ==========================================
# 3. Pipeline Orchestrator
# ==========================================


def run_qc(data_type, module_dir, zarr_path, batch_key, sample_key, qc_config):
    """Main Orchestrator."""
    module_dir = Path(module_dir)
    module_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "area_col": "Area" if data_type == "CosMx" else "cell_area",
        "has_dapi": data_type == "CosMx",
        "control_pattern": "^NegPrb|^SystemControl"
        if data_type == "CosMx"
        else "control_probe|control_codeword",
        "spatial_key": "global" if data_type == "CosMx" else "spatial",
        "nucleus_col": "Mean.DAPI" if data_type == "CosMx" else "nucleus_area",
    }

    # 1. Load Original Zarr
    logger.info(f"Loading Zarr from {zarr_path}")
    sdata_orig = sd.read_zarr(zarr_path)
    adata = sdata_orig.tables["table"].copy()
    del sdata_orig
    gc.collect()

    # 2. Integrate Proseg Zarr (If provided)
    if qc_config.proseg_zarr_path is not None:
        logger.info(f"Integrating Proseg Zarr: {qc_config.proseg_zarr_path}")
        sdata_proseg = sd.read_zarr(qc_config.proseg_zarr_path)
        adata_proseg = sdata_proseg.tables["table"].copy()
        del sdata_proseg
        gc.collect()

        # Clean IDs
        proseg_ids = (
            adata_proseg.obs[qc_config.proseg_cell_id_col].astype(str)
            if qc_config.proseg_cell_id_col in adata_proseg.obs
            else adata_proseg.obs_names.astype(str)
        )
        adata_proseg.obs_names = proseg_ids.str.replace(
            r"\.0$", "", regex=True
        ).str.replace(r"^c_\d+_", "", regex=True)

        orig_ids = (
            adata.obs[qc_config.original_cell_id_col].astype(str)
            if qc_config.original_cell_id_col in adata.obs
            else adata.obs_names.astype(str)
        )
        adata.obs_names = orig_ids.str.replace(r"\.0$", "", regex=True).str.replace(
            r"^c_\d+_", "", regex=True
        )

        common_cells = adata_proseg.obs_names.intersection(adata.obs_names)
        logger.info(f"Overlap: {len(common_cells)} cells between Original and Proseg.")

        if len(common_cells) > 0:
            adata_proseg = adata_proseg[common_cells].copy()
            adata_orig_sub = adata[common_cells].copy()

            for col in adata_orig_sub.obs.columns:
                if col not in adata_proseg.obs.columns:
                    adata_proseg.obs[col] = adata_orig_sub.obs[col].copy()

            if "spatial" in adata_proseg.obsm:
                adata_proseg.obsm["global"] = adata_proseg.obsm["spatial"].copy()

            if "volume" in adata_proseg.obs.columns:
                adata_proseg.obs[cfg["area_col"]] = adata_proseg.obs["volume"]
            elif "area" in adata_proseg.obs.columns:
                adata_proseg.obs[cfg["area_col"]] = adata_proseg.obs["area"]

            adata = adata_proseg.copy()
            del adata_proseg, adata_orig_sub
            gc.collect()

    # 3. Keys and Axes
    slide_name = Path(zarr_path).stem
    for key in [batch_key, sample_key]:
        if key is not None and key not in adata.obs.columns:
            adata.obs[key] = slide_name

    if data_type == "CosMx" and cfg["spatial_key"] in adata.obsm:
        coords = adata.obsm[cfg["spatial_key"]].copy()
        coords[:, 1] = np.max(coords[:, 1]) - coords[:, 1]
        adata.obsm[cfg["spatial_key"]] = coords

    # 4. Apply Metadata CSV
    if qc_config.fov_metadata_path:
        meta_df = pd.read_csv(qc_config.fov_metadata_path)
        if qc_config.metadata_join_col in adata.obs.columns:
            logger.info(
                f"Concatenating metadata from CSV based on join column {qc_config.metadata_join_col}"
            )

            clean_col = f"{qc_config.metadata_join_col}_clean"

            adata.obs[clean_col] = (
                adata.obs[qc_config.metadata_join_col].astype(str).str.strip()
            )
            meta_df[clean_col] = (
                meta_df[qc_config.metadata_join_col].astype(str).str.strip()
            )

            meta_df.set_index(clean_col, inplace=True)

            for col in meta_df.columns:
                if col != qc_config.metadata_join_col:
                    adata.obs[col] = (
                        adata.obs[clean_col]
                        .map(meta_df[col])
                        .astype(str)
                        .astype("category")
                    )

            # Cleanup
            del adata.obs[clean_col]
        else:
            logger.warning(f"{qc_config.metadata_join_col} is not found in adata")
    else:
        logger.warning("Not adding any metadata..")

    # 5. Execute Pipeline Sequence
    adata = calculate_qc_metrics(adata, cfg)

    export_web_metrics(module_dir, adata, cfg, qc_config)
    plot_metrics(module_dir, adata, cfg, qc_config)
    plot_spatial_qc(module_dir, adata, cfg)

    adata = filter_and_normalize(adata, qc_config, cfg)

    # 6. Save
    out_path = module_dir / "adata.h5ad"
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
