"""Prepares Files needed for Web Visualization Tool."""

import gc
import glob
import json
import os
import re
import shutil
import tarfile
import warnings
from logging import getLogger

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def export_cpdb_for_web_vis(
    cpdb_out_dir: str,
    adata_path: str,
    celltype_key: str,
    microenv_key: str,
    out_dir: str,
) -> None:
    """Exports CellPhoneDB significant means AND exact cell counts per microenvironment

    for an interactive D3 web visualization.

    Args:
        cpdb_out_dir: Directory containing CellPhoneDB output text files.
        adata_path: Path to the input .h5ad file to compute cell counts.
        celltype_key: Column name in adata.obs containing cell type labels.
        microenv_key: Column name in adata.obs containing microenvironment labels.
        out_dir: Target directory to save the JSON outputs.
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse and extract edges from CellPhoneDB results
    search_pattern = f"{cpdb_out_dir}/*significant_means*.txt"
    matched_files = glob.glob(search_pattern)

    if matched_files:
        target_file = matched_files[0]
        sig_means = pd.read_csv(target_file, sep="\t")
        logger.info(f"Successfully loaded edges from: {target_file}")
    else:
        logger.info(
            f"Error: No files matching '*significant_means*.txt' found in {cpdb_out_dir}"
        )
        return

    # CellPhoneDB columns representing interacting cell pairs are formatted as 'CellTypeA|CellTypeB'
    pair_cols = [c for c in sig_means.columns if "|" in c]
    edges = []

    for _, row in sig_means.iterrows():
        interaction = row["interacting_pair"]
        for pair in pair_cols:
            val = row[pair]
            # Only track interactions with statistically significant means (> 0)
            if pd.notna(val) and val > 0:
                source, target = pair.split("|")
                edges.append(
                    {
                        "source": source,
                        "target": target,
                        "interaction": interaction,
                        "value": float(val),
                    }
                )

    edges_out_path = os.path.join(out_dir, "cpdb_edges.json")
    with open(edges_out_path, "w") as f:
        json.dump(edges, f)
    logger.info(f"Saved {len(edges)} interactions to {edges_out_path}")

    # 2. Compute explicit cell counts per microenvironment
    logger.info(f"Loading AnnData to calculate cell counts: {adata_path}...")
    adata = sc.read_h5ad(adata_path)

    # Group by microenvironment and cell type to tally the population size
    counts = (
        adata.obs.groupby([microenv_key, celltype_key]).size().reset_index(name="count")
    )

    micro_map = {}
    for _, row in counts.iterrows():
        env = str(row[microenv_key])
        cell = str(row[celltype_key])
        cnt = int(row["count"])

        if cnt > 0:  # Omit empty pairs
            if env not in micro_map:
                micro_map[env] = {}
            micro_map[env][cell] = cnt

    envs_out_path = os.path.join(out_dir, "cpdb_microenvs.json")
    with open(envs_out_path, "w") as f:
        json.dump(micro_map, f)
    logger.info(f"Saved microenvironment cell counts to {envs_out_path}")


def prepare_vitessce_data(
    input_h5ad: str,
    output_zarr: str = "adata_vitessce.zarr",
    output_json: str = "spatial_metadata.json",
    coord_key: str = "global",
    slide_col: str = "slide_id",
    sample_col: str = "sample_id",
    microenv_key: str = None,
    mod6_dir: str = None,
    mod8b_dir: str = None,
) -> tuple[sc.AnnData, dict[str, list]]:
    """Processes an AnnData object to create masked spatial coordinates for hierarchical filtering

    (by Slide, Sample, and Microenvironment) inside Vitessce, saving the result to a Zarr store.

    Args:
        input_h5ad: Path to the input .h5ad file.
        output_zarr: Path to save the resulting unified .zarr store.
        output_json: Path to save the hierarchical structural dictionary.
        coord_key: The baseline spatial coordinate key found in adata.obsm.
        slide_col: The column in adata.obs containing slide IDs.
        sample_col: The column in adata.obs containing sample IDs.
        microenv_key: The column in adata.obs containing microenvironment metadata.

    Returns:
        A tuple containing the modified AnnData object and the hierarchy dictionary.
    """
    logger.info(f"Loading {input_h5ad}...")
    adata = sc.read_h5ad(input_h5ad)
    hierarchy = {}

    # --- NEW: TAG MUSPAN ROI CELLS ---
    adata.obs["muspan_region"] = "Outside ROI"

    if mod6_dir and os.path.exists(mod6_dir):
        # Look through all sample folders in Module 6
        for sample_name in os.listdir(mod6_dir):
            sample_path = os.path.join(mod6_dir, sample_name)
            if os.path.isdir(sample_path):
                morph_files = glob.glob(
                    f"{sample_path}/**/morphometrics.csv", recursive=True
                )
                if morph_files:
                    try:
                        df_morph = pd.read_csv(morph_files[0])
                        if "Cell_ID" in df_morph.columns:
                            valid_cells_set = set(
                                df_morph["Cell_ID"].dropna().astype(str)
                            )
                            matched_indices = []

                            # Pre-filter to only cells in THIS sample to speed up the loop
                            sample_mask = adata.obs[sample_col] == sample_name

                            for idx, row in adata.obs[sample_mask].iterrows():
                                # 1. CosMx Match: FOV_CellID (e.g. '66_1')
                                if "fov" in row and "cell_ID" in row:
                                    fov_str = str(row["fov"]).replace(".0", "")
                                    cid_str = str(row["cell_ID"]).replace(".0", "")
                                    cosmx_id = f"{fov_str}_{cid_str}"

                                    if cosmx_id in valid_cells_set:
                                        matched_indices.append(idx)
                                        continue

                                # 2. Xenium / Direct Match fallback
                                idx_str = str(idx)
                                if idx_str in valid_cells_set or (
                                    "cell_id" in row
                                    and str(row["cell_id"]) in valid_cells_set
                                ):
                                    matched_indices.append(idx)

                            if matched_indices:
                                adata.obs.loc[matched_indices, "muspan_region"] = (
                                    "In ROI"
                                )
                                logger.info(
                                    f"Tagged {len(matched_indices)} cells as 'In ROI' for sample {sample_name}"
                                )
                            else:
                                logger.warning(
                                    f"Matched 0 cells for MuSpAn ROI in sample {sample_name}."
                                )
                    except Exception as e:
                        logger.warning(
                            f"Could not parse ROI cells for {sample_name}: {e}"
                        )

    # Restore log-normalized, non-negative sparse counts for Vitessce
    if adata.raw is not None:
        logger.info("Restoring log-normalized sparse counts from .raw...")
        adata.X = adata.raw.X.copy()
        del adata.raw

    # Vitessce requires CSC format to query genes efficiently.
    if sp.issparse(adata.X):
        print(
            "Converting adata.X to CSC format for efficient Vitessce gene querying..."
        )
        adata.X = adata.X.tocsc()

    # --- 1. SLIDE & SAMPLE GEOMETRIC MASKING ---
    logger.info("Generating Slide and Sample hierarchical spatial coordinate masks...")
    for slide in adata.obs[slide_col].dropna().unique():
        slide_str = str(slide)
        hierarchy[slide_str] = []

        # Mask out anything that does not belong to this specific slide
        coords_slide = adata.obsm[coord_key].copy()
        coords_slide[adata.obs[slide_col] != slide] = np.nan
        adata.obsm[f"{coord_key}_{slide_str}"] = coords_slide

        # Isolate samples belonging to this slide
        slide_mask = adata.obs[slide_col] == slide
        samples_in_slide = adata.obs[slide_mask][sample_col].dropna().unique()

        for sample in samples_in_slide:
            sample_str = str(sample)
            hierarchy[slide_str].append(sample_str)

            # Mask out anything that does not belong to this specific sample
            coords_sample = adata.obsm[coord_key].copy()
            coords_sample[adata.obs[sample_col] != sample] = np.nan
            adata.obsm[f"{coord_key}_{sample_str}"] = coords_sample

    # --- 2. MICROENVIRONMENT GEOMETRIC MASKING ---
    if microenv_key:
        logger.info("Generating Microenvironment spatial coordinate masks...")
        for env in adata.obs[microenv_key].dropna().unique():
            env_str = str(env)

            # Mask out cells that do not live within this microenvironment
            coords_env = adata.obsm[coord_key].copy()
            coords_env[adata.obs[microenv_key] != env] = np.nan
            adata.obsm[f"spatial_microenv_{env_str}"] = coords_env

    # --- NEW: INJECT LIANA+ SINGLE-CELL SCORES FOR VITESSCE ---
    if mod8b_dir and os.path.exists(mod8b_dir):
        logger.info("Injecting LIANA+ Single-Cell CCC scores into Zarr...")
        lrdata_path = os.path.join(mod8b_dir, "lrdata.h5ad")
        nmf_path = os.path.join(mod8b_dir, "nmf_adata.h5ad")

        # Helper function to min-max scale arrays safely, removing NaNs
        def min_max_scale(arr):
            arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
            v_min, v_max = arr.min(), arr.max()
            if v_max > v_min:
                return (arr - v_min) / (v_max - v_min)
            return np.zeros_like(arr)

        # 1. Inject Top Ligand-Receptor Pairs
        if os.path.exists(lrdata_path):
            try:
                lrdata = sc.read_h5ad(lrdata_path)
                top_lrs = (
                    lrdata.var.sort_values("morans", ascending=False)
                    .head(15)
                    .index.tolist()
                )
                for pair in top_lrs:
                    val_arr = (
                        lrdata[:, pair].X.toarray().flatten()
                        if sp.issparse(lrdata.X)
                        else lrdata[:, pair].X.flatten()
                    )

                    # Use Pandas Series to guarantee cell barcodes match perfectly
                    scaled_vals = min_max_scale(val_arr)
                    adata.obs[f"LR_{pair}"] = pd.Series(
                        scaled_vals, index=lrdata.obs_names
                    )
                    adata.obs[f"LR_{pair}"] = (
                        adata.obs[f"LR_{pair}"].fillna(0.0).astype(float)
                    )

                logger.info(f"Added {len(top_lrs)} Ligand-Receptor tracks to adata.obs")
            except Exception as e:
                logger.warning(f"Failed to inject LR data: {e}")

        # 2. Inject NMF CCC Factors
        if os.path.exists(nmf_path):
            try:
                nmf_adata = sc.read_h5ad(nmf_path)
                for factor in nmf_adata.var_names:
                    val_arr = (
                        nmf_adata[:, factor].X.toarray().flatten()
                        if sp.issparse(nmf_adata.X)
                        else nmf_adata[:, factor].X.flatten()
                    )

                    # Use Pandas Series to guarantee cell barcodes match perfectly
                    scaled_vals = min_max_scale(val_arr)
                    adata.obs[f"CCC_{factor}"] = pd.Series(
                        scaled_vals, index=nmf_adata.obs_names
                    )
                    adata.obs[f"CCC_{factor}"] = (
                        adata.obs[f"CCC_{factor}"].fillna(0.0).astype(float)
                    )

                logger.info(
                    f"Added {len(nmf_adata.var_names)} CCC NMF Factors to adata.obs"
                )
            except Exception as e:
                logger.warning(f"Failed to inject NMF data: {e}")

    # --- 3. EXPORT UNIFIED ARTIFACTS ---
    logger.info("Stripping heavy metadata from AnnData before Zarr export...")
    # Delete massive neighbor graphs (often >50% of the file size)
    if hasattr(adata, "obsp"):
        del adata.obsp
    if hasattr(adata, "varp"):
        del adata.varp
    if hasattr(adata, "varm"):
        del adata.varm

    # Delete heavy embeddings (PCA, scVI) that Vitessce doesn't visualize
    for k in list(adata.obsm.keys()):
        if not (
            k == coord_key
            or k.startswith(coord_key + "_")
            or k.startswith("X_umap")
            or k.startswith("spatial_microenv_")
        ):
            del adata.obsm[k]

    logger.info(f"Saving consolidated Zarr store to {output_zarr}...")
    adata.write_zarr(output_zarr)

    logger.info(f"Saving hierarchy metadata JSON to {output_json}...")
    with open(output_json, "w") as f:
        json.dump(hierarchy, f)

    logger.info("Vitessce data preparation complete!")
    return adata, hierarchy


def export_spatial_stats_for_web(mod5_dir: str, mod6_dir: str, out_dir: str) -> None:
    """
    Crawls the SpatialStat and MuSpAn output directories to collect the JSON/CSV
    files and organizes them into the web visualization data folder.
    """
    target_stats_dir = os.path.join(out_dir, "spatial_stats")
    os.makedirs(target_stats_dir, exist_ok=True)

    logger.info(f"\n--- Collecting Spatial Statistics into {target_stats_dir} ---")

    # 1. Collect SquidPy Stats (from Module 5)
    if os.path.exists(mod5_dir):
        # Folders in mod5_dir represent the samples
        for sample_name in os.listdir(mod5_dir):
            sample_path = os.path.join(mod5_dir, sample_name)
            if os.path.isdir(sample_path):
                sample_out = os.path.join(target_stats_dir, sample_name)
                os.makedirs(sample_out, exist_ok=True)

                files_to_copy = [
                    f"centrality_scores_{sample_name}.json",
                    f"co_occurrence_{sample_name}.json",
                    f"nhood_enrichment_{sample_name}.json",
                    f"moranI_results_{sample_name}.csv",
                ]

                for f in files_to_copy:
                    src = os.path.join(sample_path, f)
                    if os.path.exists(src):
                        shutil.copy(src, os.path.join(sample_out, f))
                        logger.info(f"Copied {f}")

    # 2. Collect MuSpAn Stats (from Module 6)
    if os.path.exists(mod6_dir):
        # Folders in mod6_dir represent the samples
        for sample_name in os.listdir(mod6_dir):
            sample_path = os.path.join(mod6_dir, sample_name)
            if os.path.isdir(sample_path):
                sample_out = os.path.join(target_stats_dir, sample_name)
                os.makedirs(sample_out, exist_ok=True)

                # MuSpAn outputs are nested inside the selection folder, so we use glob
                pcf_files = glob.glob(
                    f"{sample_path}/**/cross_pcf_all.json", recursive=True
                )
                morph_files = glob.glob(
                    f"{sample_path}/**/morphometrics.csv", recursive=True
                )

                if pcf_files:
                    shutil.copy(
                        pcf_files[0], os.path.join(sample_out, "cross_pcf_all.json")
                    )
                    logger.info(f"Copied cross_pcf_all.json for {sample_name}")

                if morph_files:
                    shutil.copy(
                        morph_files[0], os.path.join(sample_out, "morphometrics.csv")
                    )
                    logger.info(f"Copied morphometrics.csv for {sample_name}")


def export_qc_for_web(mod1_dir: str, out_dir: str):
    target_qc_dir = os.path.join(out_dir, "qc")
    os.makedirs(target_qc_dir, exist_ok=True)
    logger.info(f"\n--- Processing QC metrics into {target_qc_dir} ---")

    qc_csvs = glob.glob(f"{mod1_dir}/**/qc_metrics.csv", recursive=True)
    qc_jsons = glob.glob(f"{mod1_dir}/**/qc_thresholds.json", recursive=True)

    histograms_by_slide = {}
    thresholds_by_slide = {}
    all_dfs = []

    # 1. Process Thresholds by Slide
    for json_path in qc_jsons:
        # Extract the slide name from the parent directory of the JSON
        slide_name = os.path.basename(os.path.dirname(json_path))
        if slide_name == os.path.basename(mod1_dir):
            slide_name = "Slide_1"  # Fallback if run without subfolders

        with open(json_path, "r") as f:
            thresholds_by_slide[slide_name] = json.load(f)

    # 2. Process Metrics & Histograms by Slide
    for csv_path in qc_csvs:
        slide_name = os.path.basename(os.path.dirname(csv_path))
        if slide_name == os.path.basename(mod1_dir):
            slide_name = "Slide_1"

        df = pd.read_csv(csv_path)
        all_dfs.append(df)

        histograms_by_slide[slide_name] = {}
        for col in ["total_counts", "n_genes_by_counts", "area", "nucleus_signal"]:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    counts, edges = np.histogram(data, bins=100)
                    histograms_by_slide[slide_name][col] = {
                        "counts": counts.tolist(),
                        "edges": np.round(edges, 2).tolist(),
                    }

    # 3. Create an "All" aggregate
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        histograms_by_slide["All"] = {}

        for col in ["total_counts", "n_genes_by_counts", "area", "nucleus_signal"]:
            if col in combined_df.columns:
                data = combined_df[col].dropna()
                if len(data) > 0:
                    counts, edges = np.histogram(data, bins=100)
                    histograms_by_slide["All"][col] = {
                        "counts": counts.tolist(),
                        "edges": np.round(edges, 2).tolist(),
                    }

        # For "All" thresholds, we'll just inherit the first slide's config
        if thresholds_by_slide:
            thresholds_by_slide["All"] = list(thresholds_by_slide.values())[0]

    # 4. Save the nested JSONs
    with open(os.path.join(target_qc_dir, "qc_histograms.json"), "w") as f:
        json.dump(histograms_by_slide, f)

    with open(os.path.join(target_qc_dir, "qc_thresholds.json"), "w") as f:
        json.dump(thresholds_by_slide, f)


def export_composition_for_web(
    adata_path: str, out_dir: str, slide_col: str, sample_col: str, anno_keywords: list
) -> None:
    """Exports precise cell counts for interactive pie charts"""
    logger.info(f"\n--- Generating cell composition JSON from {adata_path} ---")
    adata = sc.read_h5ad(adata_path)

    composition = {}

    # Extract dynamic annotation columns to compute
    cat_cols = [
        c
        for c in adata.obs.columns
        if any(k in c for k in anno_keywords)
        or c
        in [
            "fov",
            "DiseaseType",
            "TreatmentResponse",
            "spatial_microenvironment",
            "sample_id",
            "slide_id",
        ]
    ]

    def get_counts(df):
        res = {}
        for col in cat_cols:
            if col in df.columns:
                counts = df[col].value_counts().to_dict()
                res[col] = {str(k): int(v) for k, v in counts.items()}
        return res

    # 1. Global View
    composition["All_All"] = get_counts(adata.obs)

    # 2. Slice/Sample Views
    if slide_col in adata.obs.columns:
        for slide in adata.obs[slide_col].dropna().unique():
            slide_mask = adata.obs[slide_col] == slide
            composition[f"{slide}_All"] = get_counts(adata.obs[slide_mask])

            if sample_col in adata.obs.columns:
                samples = adata.obs[slide_mask][sample_col].dropna().unique()
                for sample in samples:
                    sample_mask = adata.obs[sample_col] == sample
                    composition[f"{slide}_{sample}"] = get_counts(
                        adata.obs[sample_mask]
                    )

    out_file = os.path.join(out_dir, "cell_composition.json")
    with open(out_file, "w") as f:
        json.dump(composition, f, indent=4)
    logger.info(f"Exported cell composition to {out_file}")


def export_segmentations_for_web(
    adata_path: str, data_type: str, settings: dict, out_dir: str
):
    """
    Parses original Flat Files (*polygons.csv) OR Proseg Zarrs to extract
    polygon vertices, applying the same Y-axis flips applied during QC.
    Safely handles multiple merged slides. Outputs sub-JSONs for dynamic filtering.
    """
    logger.info("\n--- Extracting Cell Segmentations ---")
    adata = sc.read_h5ad(adata_path)

    seg_all = {}
    seg_by_slide = {}
    seg_by_sample = {}
    seg_by_microenv = {}

    # Build metadata lookup to avoid repeated .get() on DataFrames
    obs_meta = {}
    for obs_name, row in adata.obs.iterrows():
        obs_meta[obs_name] = {
            "slide": str(row.get("slide_id", "All")),
            "sample": str(row.get("sample_id", "All")),
            "microenv": str(row.get("spatial_microenvironment", "All"))
            if "spatial_microenvironment" in row
            else None,
        }

    slide_ids = (
        adata.obs["slide_id"].dropna().unique()
        if "slide_id" in adata.obs.columns
        else [None]
    )

    for slide_id in slide_ids:
        # Resolve dataset paths from settings
        slide_settings = (
            settings.get("io", {}).get("raw_data", {}).get(slide_id, {})
            if slide_id
            else settings.get("io", {})
        )
        raw_dataset_dir = slide_settings.get("dataset_dir")

        # NOTE: Pass your Proseg Zarr path in your config under this key, e.g., settings["io"]["raw_data"]["Slide1"]["proseg_zarr_dir"] = "path/to/proseg.zarr"
        proseg_zarr_dir = slide_settings.get("proseg_zarr_dir")

        logger.info(f"Extracting segmentations for slide {slide_id}...")

        # =========================================================
        # PATH A: PROSEG ZARR PARSING
        # =========================================================
        if proseg_zarr_dir and os.path.exists(proseg_zarr_dir):
            try:
                import spatialdata

                logger.info(
                    f"Using Proseg Zarr for Slide {slide_id}: {proseg_zarr_dir}"
                )
                proseg_sdata = spatialdata.read_zarr(proseg_zarr_dir)

                proseg_gdf = proseg_sdata.shapes["cell_boundaries"]
                if hasattr(proseg_gdf, "compute"):
                    proseg_gdf = proseg_gdf.compute()

                proseg_obs = proseg_sdata.tables["table"].obs
                proseg_id_map = {}

                # Build mapping from (FOV, Cell_ID) -> Proseg Global Index
                if "fov" in proseg_obs.columns and "cell_ID" in proseg_obs.columns:
                    for idx, row in proseg_obs.iterrows():
                        f = str(row["fov"])
                        c = str(row["cell_ID"]).replace(".0", "")
                        proseg_id_map[(f, c)] = idx
                elif "original_cell_id" in proseg_obs.columns:
                    for idx, row in proseg_obs.iterrows():
                        parts = str(row["original_cell_id"]).split("_")
                        if len(parts) >= 2:
                            proseg_id_map[(parts[-2], parts[-1].replace(".0", ""))] = (
                                idx
                            )

                # Map AnnData obs_names to Proseg geometries
                for obs_name, row in adata.obs.iterrows():
                    if slide_id is not None and row.get("slide_id") != slide_id:
                        continue

                    if "fov" in row and "cell_ID" in row:
                        f = str(row["fov"])
                        c = str(row["cell_ID"]).replace(".0", "")
                        proseg_idx = proseg_id_map.get((f, c))

                        if proseg_idx is not None:
                            try:
                                geom = proseg_gdf.loc[int(proseg_idx)].geometry
                            except KeyError:
                                geom = proseg_gdf.loc[str(proseg_idx)].geometry

                            if geom.geom_type == "MultiPolygon":
                                geom = max(geom.geoms, key=lambda a: a.area)
                            if geom.geom_type != "Polygon":
                                continue

                            poly_coords = np.array(geom.exterior.coords)

                            # Note: Check if Proseg Y-axis needs flipping relative to your AnnData coords
                            # If it does, calculate max_y and do: poly_coords[:, 1] = max_y - poly_coords[:, 1]

                            coords = poly_coords[:, :2].tolist()  # Ensure [x, y] format

                            seg_all[obs_name] = coords
                            m = obs_meta.get(obs_name)
                            if m:
                                if m["slide"] != "All":
                                    seg_by_slide.setdefault(m["slide"], {})[
                                        obs_name
                                    ] = coords
                                if m["sample"] != "All":
                                    seg_by_sample.setdefault(m["sample"], {})[
                                        obs_name
                                    ] = coords
                                if m["microenv"]:
                                    seg_by_microenv.setdefault(m["microenv"], {})[
                                        obs_name
                                    ] = coords

            except Exception as e:
                logger.error(f"Failed to parse Proseg Zarr for Slide {slide_id}: {e}")

        # =========================================================
        # PATH B: COSMX FLAT CSV PARSING (Fallback)
        # =========================================================
        elif (
            data_type == "CosMx" and raw_dataset_dir and os.path.exists(raw_dataset_dir)
        ):
            logger.info(f"Using Flat CSV Polygons for Slide {slide_id}")
            poly_files = glob.glob(os.path.join(raw_dataset_dir, "*polygons*.csv"))
            if not poly_files:
                continue

            df_poly = pd.read_csv(poly_files[0])
            if "cellID" in df_poly.columns:
                df_poly.rename(columns={"cellID": "cell_ID"}, inplace=True)

            lookup = {}
            for obs_name, row in adata.obs.iterrows():
                if slide_id is None or row.get("slide_id") == slide_id:
                    if "fov" in row and "cell_ID" in row:
                        lookup[
                            (str(row["fov"]), str(row["cell_ID"]).replace(".0", ""))
                        ] = obs_name

            max_y = df_poly["y_global_px"].max()
            df_poly["y_global_px"] = max_y - df_poly["y_global_px"]  # Apply Y flip

            grouped = df_poly.groupby(["fov", "cell_ID"])
            for (fov, cell_id), group in grouped:
                obs_name = lookup.get((str(fov), str(cell_id).replace(".0", "")))
                if obs_name:
                    coords = group[["x_global_px", "y_global_px"]].values.tolist()
                    seg_all[obs_name] = coords

                    m = obs_meta.get(obs_name)
                    if m:
                        if m["slide"] != "All":
                            seg_by_slide.setdefault(m["slide"], {})[obs_name] = coords
                        if m["sample"] != "All":
                            seg_by_sample.setdefault(m["sample"], {})[obs_name] = coords
                        if m["microenv"]:
                            seg_by_microenv.setdefault(m["microenv"], {})[obs_name] = (
                                coords
                            )

        # =========================================================
        # PATH C: XENIUM PARSING
        # =========================================================
        elif (
            data_type == "Xenium"
            and raw_dataset_dir
            and os.path.exists(raw_dataset_dir)
        ):
            bound_csv = os.path.join(raw_dataset_dir, "cell_boundaries.csv.gz")
            bound_pq = os.path.join(raw_dataset_dir, "cell_boundaries.parquet")

            df_poly = None
            if os.path.exists(bound_pq):
                df_poly = pd.read_parquet(bound_pq)
            elif os.path.exists(bound_csv):
                df_poly = pd.read_csv(bound_csv)

            if df_poly is not None:
                lookup = {}
                for obs_name, row in adata.obs.iterrows():
                    if slide_id is None or row.get("slide_id") == slide_id:
                        lookup[str(row.get("cell_id", obs_name))] = obs_name

                for cell_id, group in df_poly.groupby("cell_id"):
                    obs_name = lookup.get(str(cell_id))
                    if obs_name:
                        coords = group[["vertex_x", "vertex_y"]].values.tolist()
                        seg_all[obs_name] = coords

                        m = obs_meta.get(obs_name)
                        if m:
                            if m["slide"] != "All":
                                seg_by_slide.setdefault(m["slide"], {})[obs_name] = (
                                    coords
                                )
                            if m["sample"] != "All":
                                seg_by_sample.setdefault(m["sample"], {})[obs_name] = (
                                    coords
                                )
                            if m["microenv"]:
                                seg_by_microenv.setdefault(m["microenv"], {})[
                                    obs_name
                                ] = coords

    # =========================================================
    # EXPORT EVERYTHING AS COMPRESSED JSON SUBSETS
    # =========================================================
    os.makedirs(os.path.join(out_dir, "segmentations"), exist_ok=True)

    def save_seg(data_dict, filename):
        rounded_data = {k: np.round(v, 1).tolist() for k, v in data_dict.items()}
        with open(os.path.join(out_dir, "segmentations", filename), "w") as f:
            json.dump(rounded_data, f, separators=(",", ":"))

    save_seg(seg_all, "segmentations.json")

    for k, v in seg_by_slide.items():
        save_seg(v, f"segmentations_{k}.json")
    for k, v in seg_by_sample.items():
        save_seg(v, f"segmentations_{k}.json")
    for k, v in seg_by_microenv.items():
        save_seg(v, f"segmentations_microenv_{k}.json")

    logger.info(
        f"Successfully exported {len(seg_all)} cell boundaries across all subsets."
    )


def export_gene_list(adata_path, out_dir):
    adata = sc.read_h5ad(adata_path)

    genes = sorted([str(x) for x in adata.var_names])

    with open(os.path.join(out_dir, "genes.json"), "w") as f:
        json.dump(genes, f)


def export_de_analysis_for_web(mod3_dir: str, out_dir: str):
    target_de_dir = os.path.join(out_dir, "de_analysis")
    os.makedirs(target_de_dir, exist_ok=True)
    logger.info(f"\n--- Collecting and Optimizing DE Analysis into {target_de_dir} ---")
    de_metadata = {}

    top_de_files = glob.glob(f"{mod3_dir}/**/top_DEgenes_*.csv", recursive=True)

    for top_file in top_de_files:
        filename = os.path.basename(top_file)
        annotation_col = filename.replace("top_DEgenes_", "").replace(".csv", "")
        shutil.copy(top_file, os.path.join(target_de_dir, filename))

        parent_dir = os.path.dirname(top_file)
        cluster_files = glob.glob(
            os.path.join(parent_dir, "DEgenes", "cluster_*_data.csv")
        )

        clusters = []
        for c_file in cluster_files:
            c_filename = os.path.basename(c_file)
            cluster_name = c_filename.replace("cluster_", "").replace("_data.csv", "")
            clusters.append(cluster_name)

            # --- OPTIMIZATION: Convert bulky CSV to lightweight JSON array ---
            df = pd.read_csv(c_file)
            # Filter non-significant fluff to save space
            df = df[(df["pvals_adj"] < 0.1) | (abs(df["logfoldchanges"]) > 0.5)]

            clean_data = {
                "names": df["names"].tolist(),
                "logfc": np.round(df["logfoldchanges"].fillna(0), 3).tolist(),
                "pvals": df["pvals_adj"].fillna(1.0).tolist(),
            }
            json_filename = f"{annotation_col}_cluster_{cluster_name}.json"
            with open(os.path.join(target_de_dir, json_filename), "w") as f:
                json.dump(clean_data, f)

        de_metadata[annotation_col] = sorted(clusters)

    with open(os.path.join(target_de_dir, "de_metadata.json"), "w") as f:
        json.dump(de_metadata, f)


def sanitize_filename(name):
    """Ensure gene names don't break file paths (e.g., if a gene has a slash)."""
    return re.sub(r'[\\/*?:"<>|]', "_", str(name))


def export_genes_and_obs_for_web(
    adata_path: str,
    out_dir: str,
    extra_cols: list = None,
    anno_keywords: list = None,
    standard_meta: list = None,
):
    """
    1. Exports one tiny sparse JSON file per gene.
    2. Exports cell clusters arrays for the DE Analysis Violins.
    3. Exports optimized spatial locations for Multiplexing.
    """
    logger.info("\n--- Exporting Highly Optimized JSON Data ---")
    adata = sc.read_h5ad(adata_path)

    # 1. Export Sparse Genes in CHUNKS
    gene_dir = os.path.join(out_dir, "genes")
    os.makedirs(gene_dir, exist_ok=True)

    if adata.raw is not None:
        sparse_matrix = adata.raw.X
    elif "counts" in adata.layers:
        sparse_matrix = adata.layers["counts"]
    else:
        sparse_matrix = adata.X

    X = (
        sparse_matrix.tocsc()
        if sp.issparse(sparse_matrix)
        else sp.csc_matrix(sparse_matrix)
    )

    gene_list = []
    gene_to_chunk_map = {}
    chunk_size = 200

    genes = list(adata.var_names)

    for chunk_idx in range(0, len(genes), chunk_size):
        chunk_genes = genes[chunk_idx : chunk_idx + chunk_size]
        chunk_data = {}

        for gene in chunk_genes:
            i = genes.index(gene)
            safe_gene = sanitize_filename(gene)

            # Save mapping info
            gene_list.append({"original": str(gene), "safe": safe_gene})
            gene_to_chunk_map[safe_gene] = f"chunk_{chunk_idx // chunk_size}"

            col_data = X.getcol(i)
            if col_data.nnz > 0:
                chunk_data[safe_gene] = {
                    "i": col_data.indices.tolist(),
                    "v": np.round(col_data.data, 3).tolist(),
                }
            else:
                chunk_data[safe_gene] = {"i": [], "v": []}

        # Export the chunk
        with open(
            os.path.join(gene_dir, f"chunk_{chunk_idx // chunk_size}.json"), "w"
        ) as f:
            json.dump(chunk_data, f, separators=(",", ":"))

    # Save the gene list mapping
    with open(os.path.join(out_dir, "genes_list.json"), "w") as f:
        json.dump(gene_list, f)

    # Save the chunk index so the frontend knows where to look for a gene
    with open(os.path.join(out_dir, "gene_chunk_index.json"), "w") as f:
        json.dump(gene_to_chunk_map, f)

    # 2. Export Cell Clusters (Used by DE Analysis Violins)
    if anno_keywords is None:
        anno_keywords = ["leiden", "CellTypist", "sctype"]

    cat_cols = [c for c in adata.obs.columns if any(k in c for k in anno_keywords)]

    for m in standard_meta:
        if m in adata.obs.columns and m not in cat_cols:
            cat_cols.append(m)

    if extra_cols:
        for col in extra_cols:
            if col in adata.obs.columns and col not in cat_cols:
                cat_cols.append(col)

    obs_data = {}
    for col in cat_cols:
        obs_data[col] = adata.obs[col].astype(str).tolist()

    with open(os.path.join(out_dir, "cell_clusters.json"), "w") as f:
        json.dump(obs_data, f)

    # 3. Export Optimized Locations (Used by Multiplexing)
    spatial_key = "global" if "global" in adata.obsm else "spatial"

    loc_data = {
        "id": adata.obs_names.tolist(),
        "x": np.round(adata.obsm[spatial_key][:, 0], 2).tolist(),
        "y": np.round(adata.obsm[spatial_key][:, 1], 2).tolist(),
        "slide": adata.obs["slide_id"].tolist()
        if "slide_id" in adata.obs
        else ["All"] * adata.n_obs,
        "sample": adata.obs["sample_id"].tolist()
        if "sample_id" in adata.obs
        else ["All"] * adata.n_obs,
    }
    with open(os.path.join(out_dir, "locations_optimized.json"), "w") as f:
        json.dump(loc_data, f)

    logger.info(f"Exported {len(gene_list)} sparse gene files and optimized metadata.")


def generate_sankeys(zarr_path, out_dir, anno_keywords: list):
    logger.info("Loading Zarr store...")
    # Read directly from the Zarr directory instead of h5ad
    adata = ad.read_zarr(zarr_path)

    # Create the 'sankeys' folder in the current directory
    out_path = f"{out_dir}/sankeys"
    os.makedirs(out_path, exist_ok=True)

    # Find all columns that have annotation data (leiden or CellTypist)
    cols = [c for c in adata.obs.columns if any(k in c for k in anno_keywords)]
    logger.info(f"Found {len(cols)} annotation columns. Generating combinations...")

    # Loop through combinations and create JSONs
    count = 0
    for col_a in cols:
        for col_b in cols:
            if col_a == col_b:
                continue

            df = adata.obs[[col_a, col_b]].dropna().copy()

            # Prefix nodes so we don't merge identical cluster numbers
            df["source_node"] = col_a + "_C" + df[col_a].astype(str)
            df["target_node"] = col_b + "_C" + df[col_b].astype(str)

            # Count the flows
            flows = (
                df.groupby(["source_node", "target_node"])
                .size()
                .reset_index(name="value")
            )
            flows = flows[flows["value"] > 0]

            # Map to nodes and links
            unique_nodes = list(
                pd.unique(flows[["source_node", "target_node"]].values.ravel("K"))
            )
            node_map = {name: i for i, name in enumerate(unique_nodes)}

            nodes = [{"name": name} for name in unique_nodes]
            links = [
                {
                    "source": node_map[row["source_node"]],
                    "target": node_map[row["target_node"]],
                    "value": int(row["value"]),
                }
                for _, row in flows.iterrows()
            ]

            # Save to the sankeys folder
            filename = f"{out_path}/{col_a}_vs_{col_b}.json"
            with open(filename, "w") as f:
                json.dump({"nodes": nodes, "links": links}, f)

            count += 1

    logger.info(f"Done! Created {count} Sankey JSON files inside the /sankeys folder.")


def export_conditions_de_for_web(
    targeted_de_dir: str, out_dir: str, celltype_col: str, treatment_col: str
):
    """
    Parses targeted Pairwise DE outputs (Cell Type -> Comparisons) and formats
    them into highly optimized JSON/CSV files for the React web tool.
    """
    target_de_dir = os.path.join(out_dir, "conditions_de_analysis")
    os.makedirs(target_de_dir, exist_ok=True)

    # We will nest the comparisons inside this dictionary alongside the configuration
    de_metadata = {
        "config": {"celltype_col": celltype_col, "treatment_col": treatment_col},
        "comparisons": {},
    }

    logger.info(
        f"\n--- Collecting and Optimizing Conditions DE Analysis into {target_de_dir} ---"
    )

    # Path where targeted_pairwise_DE saved its results
    base_de_path = os.path.join(targeted_de_dir)

    if not os.path.exists(base_de_path):
        logger.info(f"Error: Could not find {base_de_path}")
        return

    # Loop through each Cell Type folder
    for celltype_folder in os.listdir(base_de_path):
        ct_path = os.path.join(base_de_path, celltype_folder)
        if not os.path.isdir(ct_path):
            continue

        logger.info(f"Processing cell type: {celltype_folder}")
        comparisons = []
        summary_rows = []

        # Find all the "_all_genes.csv" files which represent our comparisons
        comparison_files = glob.glob(os.path.join(ct_path, "*_all_genes.csv"))

        for c_file in comparison_files:
            filename = os.path.basename(c_file)
            comparison_name = filename.replace("_all_genes.csv", "")
            comparisons.append(comparison_name)

            # 1. Read the full DE data
            df = pd.read_csv(c_file)

            # Ensure p-values don't become absolute 0 (breaks log math in JS)
            df["pvals_adj"] = df["pvals_adj"].fillna(1.0)
            df.loc[df["pvals_adj"] < 1e-300, "pvals_adj"] = 1e-300

            # --- 2. OPTIMIZATION: Convert to lightweight JSON for Volcano Plot ---
            # Filter non-significant fluff to save space
            plot_df = df[(df["pvals_adj"] < 0.1) | (abs(df["logfoldchanges"]) > 0.5)]

            clean_data = {
                "names": plot_df["names"].tolist(),
                "logfc": np.round(plot_df["logfoldchanges"].fillna(0), 3).tolist(),
                "pvals": plot_df["pvals_adj"].tolist(),
            }

            json_filename = f"{celltype_folder}_comparison_{comparison_name}.json"
            with open(os.path.join(target_de_dir, json_filename), "w") as f:
                json.dump(clean_data, f)

            # --- 3. Extract Top Genes for the React Summary Table ---
            # Sort by significance and fold change
            sig_df = df[df["pvals_adj"] < 0.05]

            # Positive LogFC = Upregulated in the Test Group
            upregulated = sig_df[sig_df["logfoldchanges"] > 0].sort_values(
                by="logfoldchanges", ascending=False
            )
            top_up = upregulated["names"].head(5).tolist()

            # Negative LogFC = Upregulated in the Reference Group (Downregulated in Test)
            downregulated = sig_df[sig_df["logfoldchanges"] < 0].sort_values(
                by="logfoldchanges", ascending=True
            )
            top_down = downregulated["names"].head(5).tolist()

            # Append to our summary list
            summary_rows.append(
                {
                    "Comparison": comparison_name,
                    "Top Upregulated": top_up,
                    "Top Downregulated": top_down,
                }
            )

        # 4. Save the Summary Table for this Cell Type
        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_csv_path = os.path.join(
                target_de_dir, f"summary_{celltype_folder}.csv"
            )
            summary_df.to_csv(summary_csv_path, index=False)

        de_metadata["comparisons"][celltype_folder] = sorted(comparisons)

    # 5. Save the hierarchy metadata so React knows what dropdowns to generate
    with open(os.path.join(target_de_dir, "conditions_de_metadata.json"), "w") as f:
        json.dump(de_metadata, f)

    logger.info(f"Exported DE conditions for {len(de_metadata.keys())} cell types.")


def export_causal_for_web(mod8c_dir: str, out_dir: str):
    """Exports LIANA+ Condition-Specific CCC & Causal Network outputs for Web."""
    if not mod8c_dir or not os.path.exists(mod8c_dir):
        return

    target_dir = os.path.join(out_dir, "causal_ccc")
    os.makedirs(target_dir, exist_ok=True)
    logger.info(f"\n--- Exporting Causal CCC into {target_dir} ---")

    # Find all generated networks recursively across all comparison subfolders
    network_files = glob.glob(
        os.path.join(mod8c_dir, "**", "causal_net_*.csv"), recursive=True
    )

    metadata_map = {}

    for net_path in network_files:
        filename = os.path.basename(net_path)
        # Extract the comparison and cell types from the filename
        match = re.search(r"causal_net_(.*?)_(.*?)_to_(.*?)\.csv", filename)
        if not match:
            continue

        comp_name, source_ct, target_ct = match.groups()

        causal_data = {
            "lr_interactions": [],
            "tf_activities": [],
            "network": {"nodes": [], "edges": []},
        }

        # Determine the parent subfolder where this specific network lives
        comp_folder = os.path.dirname(net_path)

        # 1. Ligand-Receptor interactions
        lr_path = os.path.join(comp_folder, f"liana_lr_{comp_name}.csv")
        real_source = source_ct
        real_target = target_ct

        if os.path.exists(lr_path):
            df_lr = pd.read_csv(lr_path)

            # Find the actual original cell type names by sanitizing to match
            def sanitize_name(n):
                return re.sub(r"[^\w\s-]", "", str(n)).replace(" ", "_")

            df_lr["safe_source"] = df_lr["source"].apply(sanitize_name)
            df_lr["safe_target"] = df_lr["target"].apply(sanitize_name)

            pair_df = df_lr[
                (df_lr["safe_source"] == source_ct)
                & (df_lr["safe_target"] == target_ct)
            ]
            if not pair_df.empty:
                real_source = pair_df.iloc[0]["source"]
                real_target = pair_df.iloc[0]["target"]

                pair_df = pair_df.sort_values(
                    "interaction_stat", ascending=False, key=abs
                ).head(50)
                for _, row in pair_df.iterrows():
                    causal_data["lr_interactions"].append(
                        {
                            "source": str(row["source"]),
                            "target": str(row["target"]),
                            "ligand": str(row["ligand"]),
                            "receptor": str(row["receptor"]),
                            "stat": float(row.get("interaction_stat", 0)),
                            "pval": float(row.get("padj", row.get("pvalue", 1))),
                        }
                    )

        # 2. Transcription Factor Activities
        tf_path = os.path.join(comp_folder, f"tf_estimates_{comp_name}.csv")
        if os.path.exists(tf_path):
            df_tf = pd.read_csv(tf_path, index_col=0)
            if real_target in df_tf.index:
                row = df_tf.loc[real_target]
                top_tfs = (
                    row[row.abs() > 0].sort_values(ascending=False, key=abs).head(30)
                )
                for tf, stat in top_tfs.items():
                    causal_data["tf_activities"].append(
                        {
                            "cell_type": str(real_target),
                            "tf": str(tf),
                            "stat": float(stat),
                        }
                    )

        # Step 3. Causal Network Edges
        df_net = pd.read_csv(net_path)

        node_dict = {}

        # Helper to map the CSV types directly
        def map_node_type(raw_type):
            if raw_type == "input":
                return "Receptor"
            if raw_type == "output":
                return "TF"
            return "Kinase/Protein"

        for _, row in df_net.iterrows():
            if pd.isna(row["source"]) or pd.isna(row["target"]):
                continue

            src = str(row["source"])
            tgt = str(row["target"])

            src_type = map_node_type(str(row.get("source_type", "unmeasured")))
            tgt_type = map_node_type(str(row.get("target_type", "unmeasured")))

            # Prioritize assigning Receptor/TF if a node acts as multiple types
            if src not in node_dict or src_type != "Kinase/Protein":
                node_dict[src] = src_type
            if tgt not in node_dict or tgt_type != "Kinase/Protein":
                node_dict[tgt] = tgt_type

            causal_data["network"]["edges"].append(
                {
                    "source": src,
                    "target": tgt,
                    # Safe fallback in case the column names slightly differ
                    "weight": float(row.get("target_weight", row.get("weight", 1.0))),
                    "sign": int(float(row.get("edge_type", 1))),
                }
            )

        for n, t in node_dict.items():
            causal_data["network"]["nodes"].append({"id": n, "type": t})

        # Save specific JSON for this combo to the React Data Directory
        out_json_name = f"causal_data_{comp_name}_{source_ct}_to_{target_ct}.json"
        with open(os.path.join(target_dir, out_json_name), "w") as f:
            json.dump(causal_data, f, separators=(",", ":"))

        # Append to metadata so React dropdowns know this combination exists
        display_comp = comp_name.replace("_vs_", " vs ")
        if display_comp not in metadata_map:
            metadata_map[display_comp] = []
        metadata_map[display_comp].append(
            {
                "source": str(real_source),
                "target": str(real_target),
                "file": out_json_name,
            }
        )

    # Save metadata mapping file
    with open(os.path.join(target_dir, "causal_metadata.json"), "w") as f:
        json.dump(metadata_map, f, indent=4)

    logger.info(f"Exported Causal CCC metadata mapping {len(metadata_map)} conditions.")


def run_webvisprep(
    module_dir,
    input_adata_path,
    batch_key,
    sample_key,
    celltype_key,
    microenv_key,
    module_1_dir,
    module_3_dir,
    module_5_dir,
    module_6_dir,
    module_7_dir,
    module_8_dir,
    module_8b_dir,
    module_8c_dir,
    module_9_dir,
    analysis_name,
    spatial_key,
    data_type,
    settings,
    DEAnalysis,
    anno_keywords,
    standard_meta,
):
    cpdb_adata_path = f"{module_8_dir}/adata.h5ad"

    # --- NEW: Sync Microenvironments to the Main AnnData ---
    # We must copy the microenvironment column into the main Mod 7 adata
    # so Vitessce, Segmentations, and cell composition functions have access to it!
    try:
        logger.info(
            f"Checking if main adata {input_adata_path} has microenvironment metadata..."
        )
        main_adata = sc.read_h5ad(input_adata_path, backed="r")

        if microenv_key and microenv_key not in main_adata.obs.columns:
            if os.path.exists(cpdb_adata_path):
                logger.info(f"Transferring {microenv_key} from CPDB to main AnnData...")

                # Load fully into memory so we can save changes
                main_adata = main_adata.to_memory()
                cpdb_adata = sc.read_h5ad(cpdb_adata_path, backed="r")

                if microenv_key in cpdb_adata.obs:
                    main_adata.obs[microenv_key] = cpdb_adata.obs[microenv_key].copy()
                    main_adata.write_h5ad(input_adata_path)
                    logger.info("Successfully synced microenvironments to main adata.")
        del main_adata
        gc.collect()
    except Exception as e:
        logger.warning(f"Could not sync microenvironments: {e}")

    # Step 1: Export D3 JSON Files for CellPhoneDB Network Analytics
    if os.path.exists(f"{module_8_dir}/cpdb_out"):
        export_cpdb_for_web_vis(
            cpdb_out_dir=f"{module_8_dir}/cpdb_out",
            adata_path=cpdb_adata_path,
            celltype_key=celltype_key,
            microenv_key=microenv_key,
            out_dir=module_dir,
        )

    logger.info("\n" + "=" * 50 + "\n")

    # Step 2: Prepare the unified Vitessce Spatial Zarr and Metadata Structure
    # Note: Using the single combined function prevents opening/saving Zarr multiple times!
    adata_processed, hierarchy_dict = prepare_vitessce_data(
        input_h5ad=input_adata_path,
        output_zarr=f"{module_dir}/adata_{analysis_name}.zarr",
        output_json=f"{module_dir}/spatial_metadata_{analysis_name}.json",
        coord_key=spatial_key,
        slide_col=batch_key,
        sample_col=sample_key,  # sample_id
        microenv_key=microenv_key,
        mod6_dir=module_6_dir,
        mod8b_dir=module_8b_dir,
    )

    # do for the transcription activity scores
    tf_input_adata_path = f"{module_7_dir}/tf_activity_scores.h5ad"
    if os.path.exists(tf_input_adata_path):
        _, _ = prepare_vitessce_data(
            input_h5ad=tf_input_adata_path,
            output_zarr=f"{module_dir}/adata_{analysis_name}_tf.zarr",
            output_json=f"{module_dir}/spatial_metadata_{analysis_name}.json",
            coord_key=spatial_key,
            slide_col=batch_key,
            sample_col=sample_key,
        )
    else:
        logger.warning("Module 7 data not found. Skipping TF Zarr export.")

    # Step 3: Export Spatial Statistics JSONs and CSVs
    export_spatial_stats_for_web(
        mod5_dir=module_5_dir, mod6_dir=module_6_dir, out_dir=module_dir
    )

    # Step 4: Export Quality Control Data
    export_qc_for_web(mod1_dir=module_1_dir, out_dir=module_dir)

    # Step 5: Export Cell compositions for interactuve explorer
    export_composition_for_web(
        adata_path=input_adata_path,
        out_dir=module_dir,
        slide_col=batch_key,
        sample_col=sample_key,
        anno_keywords=anno_keywords,
    )

    # Step 6: Export Segmentations (Polygons) directly from flat files
    export_segmentations_for_web(
        adata_path=input_adata_path,
        data_type=data_type,
        settings=settings,
        out_dir=module_dir,
    )

    # Step 7 & 8: Export Highly Optimized JSONs for React
    logger.info("\n--- Exporting Optimized Web Data ---")
    extra_export_cols = []
    if DEAnalysis:
        celltype_col = settings["modules"]["DEAnalysis"].get("celltype_col")
        treatment_col = settings["modules"]["DEAnalysis"].get("treatment_col")
        extra_export_cols = [celltype_col, treatment_col]

    export_genes_and_obs_for_web(
        input_adata_path,
        module_dir,
        extra_cols=extra_export_cols,
        anno_keywords=anno_keywords,
        standard_meta=standard_meta,
    )

    # Step 9: Generate Sankeys
    logger.info("\n--- Exporting Json for Sankey plot ---")
    generate_sankeys(
        zarr_path=f"{module_dir}/adata_{analysis_name}.zarr",
        out_dir=module_dir,
        anno_keywords=anno_keywords,
    )

    logger.info("\n--- Exporting Cell Type DE Analysis ---")
    export_de_analysis_for_web(mod3_dir=module_3_dir, out_dir=module_dir)

    # Step 10: Export Highly Optimized JSONs for React SPECIFIC FOR TREATMENT TYPE
    if DEAnalysis:
        export_conditions_de_for_web(
            targeted_de_dir=module_9_dir,
            out_dir=module_dir,
            celltype_col=celltype_col,
            treatment_col=treatment_col,
        )

        # Step 11: Export LIANA Causal Network
        export_causal_for_web(mod8c_dir=module_8c_dir, out_dir=module_dir)

    # Step 12: Archive the output directory using native Python tarfile
    logger.info(f"\n--- Archiving {module_dir} to {module_dir}.tar ---")
    try:
        # Resolve absolute paths to prevent issues with relative paths in container environments
        abs_module_dir = os.path.abspath(module_dir)
        tar_filename = f"{abs_module_dir.rstrip('/')}.tar"

        with tarfile.open(tar_filename, "w") as tar:
            tar.add(abs_module_dir, arcname=os.path.basename(abs_module_dir))

        logger.info(f"Successfully created archive: {tar_filename}")
    except Exception as e:
        logger.warning(f"Failed to create tar archive in Apptainer: {e}")
