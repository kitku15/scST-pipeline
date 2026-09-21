"""Quality control module."""

import gc
import json
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import spatialdata as sd
import squidpy as sq

matplotlib.use("Agg")


warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_qc(
    data_type,
    module_dir,
    zarr_path,
    min_counts,
    min_cells,
    min_genes,
    min_area,
    max_area,
    min_dapi,
    batch_key,
    sample_key,
    fov_metadata_path=None,
    proseg_zarr_path=None,
    proseg_cell_id_col="original_cell_id",
    original_cell_id_col="index",
):
    """
    Executes the quality control, filtering, and normalization pipeline for spatial data.
    """
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

    # 2. Load Original Data
    logger.info(f"Loading original Zarr from {zarr_path}")
    sdata_orig = sd.read_zarr(zarr_path)
    adata_orig = sdata_orig.tables["table"].copy()
    del sdata_orig
    gc.collect()

    # For using proseg data, we need to match the cell IDs between the original and proseg datasets.
    # this is only tested on CosMx data and may not work with others still needs to be tested
    if proseg_zarr_path is not None:
        logger.info(f"Proseg Zarr provided: {proseg_zarr_path}. Integrating...")

        sdata_proseg = sd.read_zarr(proseg_zarr_path)
        adata = sdata_proseg.tables["table"].copy()
        del sdata_proseg
        gc.collect()

        # 1. Setup Proseg IDs
        if proseg_cell_id_col in adata.obs.columns:
            proseg_ids = adata.obs[proseg_cell_id_col].astype(str)
        else:
            logger.warning(
                f"Proseg ID column '{proseg_cell_id_col}' not found. Using index."
            )
            proseg_ids = adata.obs_names.astype(str)

        # Clean: Strip '.0' and strip NanoString slide prefixes (e.g., 'c_1_', 'c_3_')
        adata.obs_names = proseg_ids.str.replace(r"\.0$", "", regex=True).str.replace(
            r"^c_\d+_", "", regex=True
        )
        adata.obs.index.name = None

        # 2. Setup Original IDs
        orig_index_backup = adata_orig.obs_names.copy()
        if (
            original_cell_id_col != "index"
            and original_cell_id_col in adata_orig.obs.columns
        ):
            orig_ids = adata_orig.obs[original_cell_id_col].astype(str)
        else:
            if original_cell_id_col != "index":
                logger.warning(
                    f"Original ID column '{original_cell_id_col}' not found. Falling back to index."
                )
            orig_ids = adata_orig.obs_names.astype(str)

        # Clean: Strip '.0' and strip NanoString slide prefixes (e.g., 'c_1_', 'c_3_')
        adata_orig.obs_names = orig_ids.str.replace(
            r"\.0$", "", regex=True
        ).str.replace(r"^c_\d+_", "", regex=True)
        adata_orig.obs.index.name = None

        # 3. Check Overlap
        common_cells = adata.obs_names.intersection(adata_orig.obs_names)
        logger.info(
            f"Original cells: {adata_orig.n_obs} | Proseg cells: {adata.n_obs} | Overlapping: {len(common_cells)}"
        )

        if len(common_cells) > 0:
            # Subset both to common cells
            adata = adata[common_cells].copy()
            adata_orig = adata_orig[common_cells].copy()

            # Transfer ALL metadata from adata_orig to our Proseg adata
            logger.info("Transferring metadata from original Zarr to Proseg data...")
            for col in adata_orig.obs.columns:
                if col not in adata.obs.columns:
                    adata.obs[col] = adata_orig.obs[col].copy()

            # Ensure spatial coordinates are where your pipeline expects them
            if "spatial" in adata.obsm:
                logger.info("Setting Proseg spatial coords to 'global' key...")
                adata.obsm["global"] = adata.obsm["spatial"].copy()

            # Standardize Area column
            if "volume" in adata.obs.columns:
                adata.obs[cfg["area_col"]] = adata.obs["volume"]
            elif "area" in adata.obs.columns:
                adata.obs[cfg["area_col"]] = adata.obs["area"]

        else:
            logger.warning(
                "Could not match Proseg cells to Original cells! Reverting to Original data."
            )
            adata_orig.obs_names = orig_index_backup
            adata_orig.obs.index.name = None
            adata = adata_orig.copy()

    else:
        adata = adata_orig.copy()

    del adata_orig
    gc.collect()

    # Map metadata from CSV if provided
    if fov_metadata_path is not None:
        logger.info(f"Loading FOV metadata from {fov_metadata_path}")
        meta_df = pd.read_csv(fov_metadata_path)

        # Standardize 'fov' column naming
        if "fov" not in adata.obs.columns and "FOV" in adata.obs.columns:
            adata.obs["fov"] = adata.obs["FOV"]

        if "fov" in adata.obs.columns:
            meta_df["fov_clean"] = pd.to_numeric(meta_df["fov"], errors="coerce")
            adata.obs["fov_clean"] = pd.to_numeric(adata.obs["fov"], errors="coerce")
            
            meta_df.set_index("fov_clean", inplace=True)

            for col in meta_df.columns:
                if col != "fov":  
                    adata.obs[col] = adata.obs["fov_clean"].map(meta_df[col]).astype(str).astype("category")
            
            if "DiseaseType" in adata.obs.columns:
                mapped_count = adata.obs["DiseaseType"].notna().sum()
                logger.info(f"VERIFICATION: {mapped_count} out of {adata.n_obs} cells received metadata!")

            del adata.obs["fov_clean"]
            logger.info(f"Successfully mapped metadata columns: {list(meta_df.columns)}")

    # Check if the requested batch/sample keys exist. If not, create them.
    for key in [batch_key, sample_key]:
        if key is not None and key not in adata.obs.columns:
            logger.info(
                f"Column '{key}' not found. Creating it and labeling all cells as 'Slide_1'."
            )
            adata.obs[key] = "Slide_1"

    # 3. Harmonize CosMx to Xenium logic
    if data_type == "CosMx":
        spatial_key = cfg["spatial_key"]
        if spatial_key in adata.obsm:
            coords = adata.obsm[spatial_key].copy()
            coords[:, 1] = np.max(coords[:, 1]) - coords[:, 1]
            adata.obsm[spatial_key] = coords
            logger.info(f"Permanently flipped Y-axis coordinates for {data_type}.")

        neg_probes = adata.var_names.str.contains("^NegPrb", case=False)
        sys_controls = adata.var_names.str.contains("^SystemControl", case=False)
        print(
            f"CosMx: Found {neg_probes.sum()} negative probes and {sys_controls.sum()} system controls in the data."
        )

        import scipy.sparse as sp

        if sp.issparse(adata.X):
            adata.obs["control_probe_counts"] = np.array(
                adata[:, neg_probes].X.sum(axis=1)
            ).flatten()
            adata.obs["control_codeword_counts"] = np.array(
                adata[:, sys_controls].X.sum(axis=1)
            ).flatten()
        else:
            adata.obs["control_probe_counts"] = adata[:, neg_probes].X.sum(axis=1)
            adata.obs["control_codeword_counts"] = adata[:, sys_controls].X.sum(axis=1)

        # Proxy for nucleus plot
        if "Mean.DAPI" not in adata.obs.columns:
            cfg["has_dapi"] = False
            cfg["nucleus_col"] = cfg["area_col"]
            adata.obs["nucleus_area"] = adata.obs[cfg["area_col"]]
        else:
            adata.obs["nucleus_area"] = adata.obs["Area"]

    # 4. Calculate QC metrics
    adata.var["control"] = adata.var_names.str.contains(
        cfg["control_pattern"], case=False, na=False
    )
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["control"], percent_top=(10, 20, 50, 150), inplace=True
    )

    # 5. Logging
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
    dapi_to_plot = min_dapi if cfg["has_dapi"] else None
    # We pass all thresholds to the plotting function so it can draw lines for them
    plot_metrics(
        module_dir, adata, cfg, min_counts, min_genes, min_area, max_area, dapi_to_plot
    )
    plot_spatial_qc(module_dir, adata, cfg)

    # Export QC metrics for Web UI
    logger.info("Exporting raw QC metrics to CSV and JSON for web visualization...")
    try:
        thresholds = {
            "min_counts": min_counts,
            "min_genes": min_genes,
            "min_cells": min_cells,
            "min_dapi": dapi_to_plot,
            "min_area": min_area,
            "max_area": max_area,
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

        if cfg["has_dapi"]:
            df_qc["nucleus_signal"] = adata.obs["Mean.DAPI"]
        else:
            df_qc["nucleus_signal"] = adata.obs["nucleus_area"] / adata.obs[
                cfg["area_col"]
            ].replace(0, np.nan)

        df_qc.to_csv(module_dir / "qc_metrics.csv", index=False)
        logger.info("Successfully exported QC metrics.")
    except Exception as e:
        logger.warning(f"Failed to export QC metrics for web: {e}")

    # Detailed Filtering and Normalization
    logger.info("--- Starting Cell Filtering ---")
    start_cells = adata.n_obs
    logger.info(f"Initial cell count: {start_cells}")

    # filter by dapi
    if cfg["has_dapi"] and min_dapi is not None:
        pre = adata.n_obs
        adata = adata[adata.obs["Mean.DAPI"] >= min_dapi].copy()
        logger.info(f"DAPI filter (>= {min_dapi}): removed {pre - adata.n_obs} cells.")

    # filter by area
    if min_area is not None and max_area is not None:
        pre = adata.n_obs
        adata = adata[
            (adata.obs[cfg["area_col"]] >= min_area)
            & (adata.obs[cfg["area_col"]] <= max_area)
        ].copy()
        logger.info(
            f"Area filter ({min_area} - {max_area}): removed {pre - adata.n_obs} cells."
        )

    # Filter by min counts
    if min_counts is not None:
        pre = adata.n_obs
        adata = adata[adata.obs["total_counts"] >= min_counts].copy()
        logger.info(
            f"Min transcripts filter (>= {min_counts}): removed {pre - adata.n_obs} cells."
        )

    # Filter by unique genes
    if min_genes is not None:
        pre = adata.n_obs
        adata = adata[adata.obs["n_genes_by_counts"] >= min_genes].copy()
        logger.info(
            f"Min unique genes filter (>= {min_genes}): removed {pre - adata.n_obs} cells."
        )

    logger.info("--- Cell Filtering Complete ---")
    logger.info(
        f"Final cell count: {adata.n_obs} (Retained {(adata.n_obs / start_cells) * 100:.1f}%)"
    )

    # Filter genes
    logger.info("--- Starting Gene Filtering ---")
    start_genes = adata.n_vars
    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
        logger.info(
            f"Min cells per gene filter (>= {min_cells}): removed {start_genes - adata.n_vars} genes."
        )
    logger.info(f"Final gene count: {adata.n_vars}")

    logger.info("Normalize data...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, base=2)

    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata

    # sc.pp.scale(adata, max_value=10)
    adata.obs.index.name = None

    # Save
    module_dir = Path(module_dir)
    out_path = module_dir / "adata.h5ad"
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
    logger.info("Quality control completed successfully.")


def plot_metrics(
    module_dir, adata, cfg, min_counts, min_genes, min_area, max_area, min_dapi
):
    """
    Generates and saves summary histograms of key cell-level QC metrics.
    Automatically zooms the x-axis to the 99.5th percentile to hide extreme outliers.
    """
    module_dir.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    # Helper function to crop outliers from the plot view
    def get_upper_limit(series):
        p99 = np.percentile(series.dropna(), 99.5)
        return max(p99, 1)

    # 1. Total transcripts
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0], color="blue")
    if min_counts is not None:
        axs[0].axvline(min_counts, color="red", linestyle="--", linewidth=2)
    axs[0].set_xlim(0, get_upper_limit(adata.obs["total_counts"]))

    # 2. Unique transcripts
    axs[1].set_title("Unique genes per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1], color="green")
    if min_genes is not None:
        axs[1].axvline(min_genes, color="red", linestyle="--", linewidth=2)
    axs[1].set_xlim(0, get_upper_limit(adata.obs["n_genes_by_counts"]))

    # 3. Cell Area
    axs[2].set_title("Cell Area (Total)")
    sns.histplot(adata.obs[cfg["area_col"]], ax=axs[2], color="orange")
    if min_area is not None:
        axs[2].axvline(min_area, color="red", linestyle="--", linewidth=2)
    if max_area is not None:
        axs[2].axvline(max_area, color="red", linestyle="--", linewidth=2)
    axs[2].set_xlim(0, get_upper_limit(adata.obs[cfg["area_col"]]))

    # 4. Nucleus / DAPI plot
    if cfg["has_dapi"]:
        axs[3].set_title("Mean DAPI (Nucleus Signal)")
        sns.histplot(adata.obs["Mean.DAPI"], ax=axs[3], color="purple")
        if min_dapi is not None:
            axs[3].axvline(min_dapi, color="red", linestyle="--", linewidth=2)
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
    """
    Generates and saves spatial scatter plots of structural cell metrics.
    """
    module_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = module_dir

    logger.info(f"Visualize {cfg['nucleus_col']!s} on tissue...")

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
