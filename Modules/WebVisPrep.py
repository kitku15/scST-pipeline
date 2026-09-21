"""Module 10: Prepares Zarr and Auxiliary Data for the FastAPI Backend."""

import os
import gc
import json
import glob
import shutil
import tarfile
import warnings
from logging import getLogger
from pathlib import Path
import anndata as ad
import re

import scanpy as sc
import numpy as np
import pandas as pd
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
    os.makedirs(out_dir, exist_ok=True)
    search_pattern = f"{cpdb_out_dir}/*significant_means*.txt"
    matched_files = glob.glob(search_pattern)

    if matched_files:
        target_file = matched_files[0]
        sig_means = pd.read_csv(target_file, sep="\t")
        logger.info(f"Successfully loaded edges from: {target_file}")
    else:
        logger.info(f"Error: No files matching '*significant_means*.txt' found in {cpdb_out_dir}")
        return

    pair_cols = [c for c in sig_means.columns if "|" in c]
    edges = []

    for _, row in sig_means.iterrows():
        interaction = row["interacting_pair"]
        for pair in pair_cols:
            val = row[pair]
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

    logger.info(f"Loading AnnData to calculate cell counts: {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    counts = adata.obs.groupby([microenv_key, celltype_key]).size().reset_index(name="count")

    micro_map = {}
    for _, row in counts.iterrows():
        env = str(row[microenv_key])
        cell = str(row[celltype_key])
        cnt = int(row["count"])

        if cnt > 0:
            if env not in micro_map:
                micro_map[env] = {}
            micro_map[env][cell] = cnt

    envs_out_path = os.path.join(out_dir, "cpdb_microenvs.json")
    with open(envs_out_path, "w") as f:
        json.dump(micro_map, f)


def export_spatial_stats_for_web(mod5_dir: str, mod6_dir: str, out_dir: str) -> None:
    target_stats_dir = os.path.join(out_dir, "spatial_stats")
    os.makedirs(target_stats_dir, exist_ok=True)

    if os.path.exists(mod5_dir):
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

    if os.path.exists(mod6_dir):
        for sample_name in os.listdir(mod6_dir):
            sample_path = os.path.join(mod6_dir, sample_name)
            if os.path.isdir(sample_path):
                sample_out = os.path.join(target_stats_dir, sample_name)
                os.makedirs(sample_out, exist_ok=True)

                pcf_files = glob.glob(f"{sample_path}/**/cross_pcf_all.json", recursive=True)
                morph_files = glob.glob(f"{sample_path}/**/morphometrics.csv", recursive=True)

                if pcf_files:
                    shutil.copy(pcf_files[0], os.path.join(sample_out, "cross_pcf_all.json"))
                if morph_files:
                    shutil.copy(morph_files[0], os.path.join(sample_out, "morphometrics.csv"))


def export_qc_for_web(mod1_dir: str, out_dir: str):
    target_qc_dir = os.path.join(out_dir, "qc")
    os.makedirs(target_qc_dir, exist_ok=True)

    qc_csvs = glob.glob(f"{mod1_dir}/**/qc_metrics.csv", recursive=True)
    qc_jsons = glob.glob(f"{mod1_dir}/**/qc_thresholds.json", recursive=True)

    histograms_by_slide = {}
    thresholds_by_slide = {}
    all_dfs = []

    for json_path in qc_jsons:
        slide_name = os.path.basename(os.path.dirname(json_path))
        if slide_name == os.path.basename(mod1_dir):
            slide_name = "Slide_1"
        with open(json_path, "r") as f:
            thresholds_by_slide[slide_name] = json.load(f)

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
        if thresholds_by_slide:
            thresholds_by_slide["All"] = list(thresholds_by_slide.values())[0]

    with open(os.path.join(target_qc_dir, "qc_histograms.json"), "w") as f:
        json.dump(histograms_by_slide, f)
    with open(os.path.join(target_qc_dir, "qc_thresholds.json"), "w") as f:
        json.dump(thresholds_by_slide, f)


def export_segmentations_for_web(adata_path: str, data_type: str, settings: dict, out_dir: str):
    logger.info("\n--- Extracting Cell Segmentations ---")
    adata = sc.read_h5ad(adata_path)
    spatial_key = settings.get("project", {}).get("spatial_key", "global" if data_type == "CosMx" else "spatial")

    seg_all = {}
    seg_by_slide = {}
    seg_by_sample = {}
    seg_by_microenv = {}

    obs_meta = {}
    for obs_name, row in adata.obs.iterrows():
        obs_meta[obs_name] = {
            "slide": str(row.get("slide_id", "All")),
            "sample": str(row.get("sample_id", "All")),
            "microenv": str(row.get("spatial_microenvironment", "All")) if "spatial_microenvironment" in row else None,
        }

    slide_ids = adata.obs["slide_id"].dropna().unique() if "slide_id" in adata.obs.columns else [None]

    for slide_id in slide_ids:
        slide_settings = settings.get("io", {}).get("raw_data", {}).get(slide_id, {}) if slide_id else settings.get("io", {})
        raw_dataset_dir = slide_settings.get("dataset_dir")
        proseg_zarr_dir = slide_settings.get("proseg_zarr_dir")

        slide_polys = {}

        # PATH A: PROSEG ZARR PARSING
        if proseg_zarr_dir and os.path.exists(proseg_zarr_dir):
            try:
                import spatialdata
                logger.info(f"Slide {slide_id}: Extracting polygons from Proseg Zarr...")
                proseg_sdata = spatialdata.read_zarr(proseg_zarr_dir)
                proseg_gdf = proseg_sdata.shapes["cell_boundaries"]
                if hasattr(proseg_gdf, "compute"):
                    proseg_gdf = proseg_gdf.compute()

                proseg_obs = proseg_sdata.tables["table"].obs
                proseg_id_map = {}

                if "fov" in proseg_obs.columns and "cell_ID" in proseg_obs.columns:
                    for idx, row in proseg_obs.iterrows():
                        proseg_id_map[(str(row["fov"]), str(row["cell_ID"]).replace(".0", ""))] = idx
                elif "original_cell_id" in proseg_obs.columns:
                    for idx, row in proseg_obs.iterrows():
                        parts = str(row["original_cell_id"]).split("_")
                        if len(parts) >= 2:
                            proseg_id_map[(parts[-2], parts[-1].replace(".0", ""))] = idx

                for obs_name, row in adata.obs.iterrows():
                    if slide_id is not None and row.get("slide_id") != slide_id: continue
                    if "fov" in row and "cell_ID" in row:
                        f, c = str(row["fov"]), str(row["cell_ID"]).replace(".0", "")
                        proseg_idx = proseg_id_map.get((f, c))

                        if proseg_idx is not None:
                            try:
                                geom = proseg_gdf.loc[int(proseg_idx)].geometry
                            except KeyError:
                                geom = proseg_gdf.loc[str(proseg_idx)].geometry

                            if geom.geom_type == "MultiPolygon":
                                geom = max(geom.geoms, key=lambda a: a.area)
                            if geom.geom_type != "Polygon": continue

                            poly_coords = np.array(geom.exterior.coords)
                            coords = poly_coords[:, :2].tolist()
                            slide_polys[obs_name] = coords

            except Exception as e:
                logger.error(f"Failed to parse Proseg Zarr: {e}")

        # PATH B: COSMX FLAT CSV PARSING (VECTORIZED + 0.12 PIXEL SCALING)
        elif data_type == "CosMx" and raw_dataset_dir and os.path.exists(raw_dataset_dir):
            poly_files = glob.glob(os.path.join(raw_dataset_dir, "*polygons*.csv"))
            if poly_files:
                logger.info(f"Slide {slide_id}: Extracting, scaling (0.12 µm/px), and aligning polygons...")
                df_poly = pd.read_csv(poly_files[0])
                if "cellID" in df_poly.columns:
                    df_poly.rename(columns={"cellID": "cell_ID"}, inplace=True)
                
                # Format IDs to match
                df_poly["fov"] = df_poly["fov"].astype(str)
                df_poly["cell_ID"] = df_poly["cell_ID"].astype(str).str.replace(".0", "", regex=False)

                # 1. Create a dataframe of the EXACT points from adata
                adata_obs = adata.obs.copy()
                if slide_id is not None and "slide_id" in adata_obs.columns:
                    adata_obs = adata_obs[adata_obs["slide_id"] == slide_id]

                # THE FIX: Get the absolute global integer positions of these specific cells
                global_indices = adata.obs_names.get_indexer(adata_obs.index)
                
                df_adata = pd.DataFrame({
                    "obs_name": adata_obs.index,
                    "fov": adata_obs["fov"].astype(str),
                    "cell_ID": adata_obs["cell_ID"].astype(str).str.replace(".0", "", regex=False),
                    "pt_x": adata.obsm[spatial_key][global_indices, 0],
                    "pt_y": adata.obsm[spatial_key][global_indices, 1]
                })

                # 2. Calculate the original centers of the raw polygons
                df_centroids = df_poly.groupby(["fov", "cell_ID"]).agg(
                    poly_x=("x_global_px", "mean"),
                    poly_y=("y_global_px", "mean")
                ).reset_index()

                # 3. Merge adata points with polygon centroids
                df_mapping = pd.merge(df_centroids, df_adata, on=["fov", "cell_ID"], how="inner")

                # 4. Merge mapping back to the full polygon vertices
                df_poly_mapped = pd.merge(
                    df_poly, 
                    df_mapping[["fov", "cell_ID", "obs_name", "pt_x", "pt_y", "poly_x", "poly_y"]], 
                    on=["fov", "cell_ID"], 
                    how="inner"
                )

                # 5. Apply the hardware physical constant
                PIXEL_SIZE = 0.12

                # 6. Vectorized Math: Center the polygon, scale by 0.12, flip Y, and snap to the point
                df_poly_mapped["final_x"] = df_poly_mapped["pt_x"] + ((df_poly_mapped["x_global_px"] - df_poly_mapped["poly_x"]) * PIXEL_SIZE)
                df_poly_mapped["final_y"] = df_poly_mapped["pt_y"] - ((df_poly_mapped["y_global_px"] - df_poly_mapped["poly_y"]) * PIXEL_SIZE)

                # 7. Group the final coordinates back into lists
                for obs_name, group in df_poly_mapped.groupby("obs_name"):
                    slide_polys[obs_name] = group[["final_x", "final_y"]].values.tolist()

        # PATH C: XENIUM PARSING
        elif data_type == "Xenium" and raw_dataset_dir and os.path.exists(raw_dataset_dir):
            bound_csv = os.path.join(raw_dataset_dir, "cell_boundaries.csv.gz")
            bound_pq = os.path.join(raw_dataset_dir, "cell_boundaries.parquet")

            df_poly = None
            if os.path.exists(bound_pq): df_poly = pd.read_parquet(bound_pq)
            elif os.path.exists(bound_csv): df_poly = pd.read_csv(bound_csv)

            if df_poly is not None:
                lookup = {}
                for obs_name, row in adata.obs.iterrows():
                    if slide_id is None or row.get("slide_id") == slide_id:
                        lookup[str(row.get("cell_id", obs_name))] = obs_name

                for cell_id, group in df_poly.groupby("cell_id"):
                    obs_name = lookup.get(str(cell_id))
                    if obs_name:
                        coords = group[["vertex_x", "vertex_y"]].values.tolist()
                        slide_polys[obs_name] = coords

        # Append local slide_polys to the global dictionaries
        for obs_name, coords in slide_polys.items():
            seg_all[obs_name] = coords
            m = obs_meta.get(obs_name)
            if m:
                if m["slide"] != "All": seg_by_slide.setdefault(m["slide"], {})[obs_name] = coords
                if m["sample"] != "All": seg_by_sample.setdefault(m["sample"], {})[obs_name] = coords
                if m["microenv"]: seg_by_microenv.setdefault(m["microenv"], {})[obs_name] = coords

    os.makedirs(os.path.join(out_dir, "segmentations"), exist_ok=True)

    def save_seg(data_dict, filename):
        rounded_data = {k: np.round(v, 1).tolist() for k, v in data_dict.items()}
        with open(os.path.join(out_dir, "segmentations", filename), "w") as f:
            json.dump(rounded_data, f, separators=(",", ":"))

    save_seg(seg_all, "segmentations.json")
    for k, v in seg_by_slide.items(): save_seg(v, f"segmentations_{k}.json")
    for k, v in seg_by_sample.items(): save_seg(v, f"segmentations_{k}.json")
    for k, v in seg_by_microenv.items(): save_seg(v, f"segmentations_microenv_{k}.json")

    logger.info(f"Successfully exported {len(seg_all)} cell boundaries across all subsets.")


def export_de_analysis_for_web(mod3_dir: str, out_dir: str):
    target_de_dir = os.path.join(out_dir, "de_analysis")
    os.makedirs(target_de_dir, exist_ok=True)
    de_metadata = {}

    top_de_files = glob.glob(f"{mod3_dir}/**/top_DEgenes_*.csv", recursive=True)

    for top_file in top_de_files:
        filename = os.path.basename(top_file)
        annotation_col = filename.replace("top_DEgenes_", "").replace(".csv", "")
        shutil.copy(top_file, os.path.join(target_de_dir, filename))

        parent_dir = os.path.dirname(top_file)
        cluster_files = glob.glob(os.path.join(parent_dir, "DEgenes", "cluster_*_data.csv"))

        clusters = []
        for c_file in cluster_files:
            c_filename = os.path.basename(c_file)
            cluster_name = c_filename.replace("cluster_", "").replace("_data.csv", "")
            clusters.append(cluster_name)

            df = pd.read_csv(c_file)

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


def export_conditions_de_for_web(targeted_de_dir: str, out_dir: str, celltype_col: str, treatment_col: str):
    target_de_dir = os.path.join(out_dir, "conditions_de_analysis")
    os.makedirs(target_de_dir, exist_ok=True)

    de_metadata = {
        "config": {"celltype_col": celltype_col, "treatment_col": treatment_col},
        "comparisons": {},
    }

    base_de_path = os.path.join(targeted_de_dir)
    if not os.path.exists(base_de_path): return

    for celltype_folder in os.listdir(base_de_path):
        ct_path = os.path.join(base_de_path, celltype_folder)
        if not os.path.isdir(ct_path): continue

        comparisons = []
        summary_rows = []

        comparison_files = glob.glob(os.path.join(ct_path, "*_all_genes.csv"))

        for c_file in comparison_files:
            filename = os.path.basename(c_file)
            comparison_name = filename.replace("_all_genes.csv", "")
            comparisons.append(comparison_name)

            df = pd.read_csv(c_file)
            df["pvals_adj"] = df["pvals_adj"].fillna(1.0)
            df.loc[df["pvals_adj"] < 1e-300, "pvals_adj"] = 1e-300

            plot_df = df[(df["pvals_adj"] < 0.1) | (abs(df["logfoldchanges"]) > 0.5)]

            clean_data = {
                "names": plot_df["names"].tolist(),
                "logfc": np.round(plot_df["logfoldchanges"].fillna(0), 3).tolist(),
                "pvals": plot_df["pvals_adj"].tolist(),
            }

            json_filename = f"{celltype_folder}_comparison_{comparison_name}.json"
            with open(os.path.join(target_de_dir, json_filename), "w") as f:
                json.dump(clean_data, f)

            sig_df = df[df["pvals_adj"] < 0.05]
            upregulated = sig_df[sig_df["logfoldchanges"] > 0].sort_values(by="logfoldchanges", ascending=False)
            top_up = upregulated["names"].head(5).tolist()

            downregulated = sig_df[sig_df["logfoldchanges"] < 0].sort_values(by="logfoldchanges", ascending=True)
            top_down = downregulated["names"].head(5).tolist()

            summary_rows.append({"Comparison": comparison_name, "Top Upregulated": top_up, "Top Downregulated": top_down})

        if summary_rows:
            pd.DataFrame(summary_rows).to_csv(os.path.join(target_de_dir, f"summary_{celltype_folder}.csv"), index=False)

        de_metadata["comparisons"][celltype_folder] = sorted(comparisons)

    with open(os.path.join(target_de_dir, "conditions_de_metadata.json"), "w") as f:
        json.dump(de_metadata, f)


def export_causal_for_web(mod8c_dir: str, out_dir: str):
    if not mod8c_dir or not os.path.exists(mod8c_dir): return
    target_dir = os.path.join(out_dir, "causal_ccc")
    os.makedirs(target_dir, exist_ok=True)

    network_files = glob.glob(os.path.join(mod8c_dir, "**", "causal_net_*.csv"), recursive=True)
    metadata_map = {}

    for net_path in network_files:
        filename = os.path.basename(net_path)
        match = re.search(r"causal_net_(.*?)_(.*?)_to_(.*?)\.csv", filename)
        if not match: continue

        comp_name, source_ct, target_ct = match.groups()
        causal_data = {"lr_interactions": [], "tf_activities": [], "network": {"nodes": [], "edges": []}}
        comp_folder = os.path.dirname(net_path)

        lr_path = os.path.join(comp_folder, f"liana_lr_{comp_name}.csv")
        real_source, real_target = source_ct, target_ct

        if os.path.exists(lr_path):
            df_lr = pd.read_csv(lr_path)
            df_lr["safe_source"] = df_lr["source"].apply(lambda n: re.sub(r"[^\w\s-]", "", str(n)).replace(" ", "_"))
            df_lr["safe_target"] = df_lr["target"].apply(lambda n: re.sub(r"[^\w\s-]", "", str(n)).replace(" ", "_"))

            pair_df = df_lr[(df_lr["safe_source"] == source_ct) & (df_lr["safe_target"] == target_ct)]
            if not pair_df.empty:
                real_source = pair_df.iloc[0]["source"]
                real_target = pair_df.iloc[0]["target"]
                pair_df = pair_df.sort_values("interaction_stat", ascending=False, key=abs).head(50)
                for _, row in pair_df.iterrows():
                    causal_data["lr_interactions"].append({
                        "source": str(row["source"]), "target": str(row["target"]), "ligand": str(row["ligand"]),
                        "receptor": str(row["receptor"]), "stat": float(row.get("interaction_stat", 0)),
                        "pval": float(row.get("padj", row.get("pvalue", 1)))
                    })

        tf_path = os.path.join(comp_folder, f"tf_estimates_{comp_name}.csv")
        if os.path.exists(tf_path):
            df_tf = pd.read_csv(tf_path, index_col=0)
            if real_target in df_tf.index:
                row = df_tf.loc[real_target]
                for tf, stat in row[row.abs() > 0].sort_values(ascending=False, key=abs).head(30).items():
                    causal_data["tf_activities"].append({"cell_type": str(real_target), "tf": str(tf), "stat": float(stat)})

        df_net = pd.read_csv(net_path)
        node_dict = {}

        def map_node_type(raw_type):
            if raw_type == "input": return "Receptor"
            if raw_type == "output": return "TF"
            return "Kinase/Protein"

        for _, row in df_net.iterrows():
            if pd.isna(row["source"]) or pd.isna(row["target"]): continue
            src, tgt = str(row["source"]), str(row["target"])
            src_type, tgt_type = map_node_type(str(row.get("source_type", "unmeasured"))), map_node_type(str(row.get("target_type", "unmeasured")))

            if src not in node_dict or src_type != "Kinase/Protein": node_dict[src] = src_type
            if tgt not in node_dict or tgt_type != "Kinase/Protein": node_dict[tgt] = tgt_type

            causal_data["network"]["edges"].append({
                "source": src, "target": tgt, "weight": float(row.get("target_weight", row.get("weight", 1.0))),
                "sign": int(float(row.get("edge_type", 1)))
            })

        for n, t in node_dict.items(): causal_data["network"]["nodes"].append({"id": n, "type": t})

        out_json_name = f"causal_data_{comp_name}_{source_ct}_to_{target_ct}.json"
        with open(os.path.join(target_dir, out_json_name), "w") as f:
            json.dump(causal_data, f, separators=(",", ":"))

        display_comp = comp_name.replace("_vs_", " vs ")
        if display_comp not in metadata_map: metadata_map[display_comp] = []
        metadata_map[display_comp].append({"source": str(real_source), "target": str(real_target), "file": out_json_name})

    with open(os.path.join(target_dir, "causal_metadata.json"), "w") as f:
        json.dump(metadata_map, f, indent=4)


def min_max_scale(arr):
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    v_min, v_max = arr.min(), arr.max()
    if v_max > v_min:
        return (arr - v_min) / (v_max - v_min)
    return np.zeros_like(arr)


def prepare_zarr_for_fastapi(input_h5ad, output_zarr, spatial_key, mod8b_dir):
    logger.info(f"Loading {input_h5ad}...")
    adata = sc.read_h5ad(input_h5ad)

    # delete the duplicate slide_ID
    if "slide_ID" in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=["slide_ID"])

    if mod8b_dir and os.path.exists(mod8b_dir):
        logger.info("Injecting LIANA+ Single-Cell CCC scores into Zarr...")
        lrdata_path = os.path.join(mod8b_dir, "lrdata.h5ad")
        nmf_path = os.path.join(mod8b_dir, "nmf_adata.h5ad")

        if os.path.exists(lrdata_path):
            try:
                lrdata = sc.read_h5ad(lrdata_path)
                top_lrs = lrdata.var.sort_values("morans", ascending=False).head(50).index.tolist()
                for pair in top_lrs:
                    val_arr = lrdata[:, pair].X.toarray().flatten() if sp.issparse(lrdata.X) else lrdata[:, pair].X.flatten()
                    adata.obs[f"LR_{pair}"] = pd.Series(min_max_scale(val_arr), index=lrdata.obs_names).fillna(0.0).astype(float)
            except Exception as e:
                logger.warning(f"Failed to inject LR data: {e}")

        if os.path.exists(nmf_path):
            try:
                nmf_adata = sc.read_h5ad(nmf_path)
                for factor in nmf_adata.var_names:
                    val_arr = nmf_adata[:, factor].X.toarray().flatten() if sp.issparse(nmf_adata.X) else nmf_adata[:, factor].X.flatten()
                    adata.obs[f"CCC_{factor}"] = pd.Series(min_max_scale(val_arr), index=nmf_adata.obs_names).fillna(0.0).astype(float)
            except Exception as e:
                logger.warning(f"Failed to inject NMF data: {e}")

    if adata.raw is not None:
        adata.X = adata.raw.X.copy()
        del adata.raw
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsc()

    if hasattr(adata, "obsp"): del adata.obsp
    if hasattr(adata, "varp"): del adata.varp
    if hasattr(adata, "varm"): del adata.varm

    keys_to_keep = [spatial_key] + [k for k in adata.obsm.keys() if k.startswith(spatial_key + "_") or k.startswith("X_umap") or k.startswith("spatial_microenv_")]
    for k in list(adata.obsm.keys()):
        if k not in keys_to_keep:
            del adata.obsm[k]

    logger.info(f"Saving optimized Zarr to {output_zarr}...")
    adata.write_zarr(output_zarr)
    return adata


def run_web_backend_prep(
    module_dir, input_adata_path, batch_key, sample_key, celltype_key, microenv_key,
    module_1_dir, module_3_dir, module_5_dir, module_6_dir, module_7_dir, 
    module_8_dir, module_8b_dir, module_8c_dir, module_9_dir, 
    analysis_name, spatial_key, data_type, settings, DEAnalysis, anno_keywords
):
    os.makedirs(module_dir, exist_ok=True)
    aux_dir = os.path.join(module_dir, "aux_data")
    os.makedirs(aux_dir, exist_ok=True)

    zarr_filename = f"adata_{analysis_name}_web.zarr"
    zarr_path = os.path.join(module_dir, zarr_filename)
    
    # 1. Zarr Creation
    prepare_zarr_for_fastapi(input_adata_path, zarr_path, spatial_key, module_8b_dir)

    tf_input_adata_path = f"{module_7_dir}/tf_activity_scores.h5ad"
    if os.path.exists(tf_input_adata_path):
        tf_zarr_path = os.path.join(module_dir, f"adata_{analysis_name}_tf_web.zarr")
        prepare_zarr_for_fastapi(tf_input_adata_path, tf_zarr_path, spatial_key, None)

    # 2. Collect Auxiliary Data
    logger.info("\n--- Collecting Pre-Computed Analytics for Backend ---")
    if os.path.exists(f"{module_8_dir}/cpdb_out"):
        export_cpdb_for_web_vis(f"{module_8_dir}/cpdb_out", f"{module_8_dir}/adata.h5ad", celltype_key, microenv_key, aux_dir)

    export_spatial_stats_for_web(module_5_dir, module_6_dir, aux_dir)
    export_qc_for_web(module_1_dir, aux_dir)
    export_segmentations_for_web(input_adata_path, data_type, settings, aux_dir)
    export_de_analysis_for_web(module_3_dir, aux_dir)

    if DEAnalysis:
        treatment_col = settings["modules"]["DEAnalysis"].get("treatment_col")
        export_conditions_de_for_web(module_9_dir, aux_dir, celltype_key, treatment_col)
        export_causal_for_web(module_8c_dir, aux_dir)

    # 3. Write Backend Configuration
    logger.info("\n--- Writing Backend Configuration ---")
    # We load the adata purely to check if the batch/sample columns exist
    temp_adata = sc.read_h5ad(input_adata_path, backed='r')
    backend_config = {
        "zarr_filename": zarr_filename,
        "spatial_key": spatial_key,
        "slide_col": batch_key if batch_key in temp_adata.obs.columns else None,
        "sample_col": sample_key if sample_key in temp_adata.obs.columns else None
    }
    with open(os.path.join(module_dir, "dataset_config.json"), "w") as f:
        json.dump(backend_config, f, indent=4)

    # 4. Archive
    logger.info(f"\n--- Archiving {module_dir} to {module_dir}.tar ---")
    abs_module_dir = os.path.abspath(module_dir)
    with tarfile.open(f"{abs_module_dir.rstrip('/')}.tar", "w") as tar:
        tar.add(abs_module_dir, arcname=os.path.basename(abs_module_dir))