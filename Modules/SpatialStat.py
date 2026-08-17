"""Spatial statistics module."""

import json
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_spatial_statistics(
    module_dir,
    input_adata_path,
    sample_key,
    cluster_name,
    condition_key=None,
    reference_condition=None,
    skip_compute=False,
):
    """Run spatial statistics."""

    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading data...")
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    if adata.obs[cluster_name].dtype.name != "category":
        adata.obs[cluster_name] = adata.obs[cluster_name].astype("category")

    valid_clusters = list(adata.obs[cluster_name].cat.categories)
    num_clusters = len(valid_clusters)
    logger.info(f"Detected {num_clusters} unique clusters. Adjusting plot sizes...")

    cent_width = max(16.0, num_clusters * 1.5)
    cent_figsize = (cent_width, max(5.0, num_clusters * 0.3))
    co_size = max(10.0, num_clusters * 4)
    co_figsize = (co_size, 10.0)
    nhood_size = max(8.0, num_clusters * 0.6)
    nhood_figsize = (nhood_size, nhood_size)

    do_aggregation = condition_key is not None and condition_key in adata.obs.columns
    if do_aggregation:
        agg_nhood = {}
        agg_centrality = []
        agg_moran = []
    else:
        logger.warning(
            "No valid condition_key provided. Skipping condition-level aggregations."
        )

    # --- PER-SAMPLE LOOP ---
    for sample in adata.obs[sample_key].unique():
        logger.info(f"--- Processing Spatial Statistics for Sample: {sample} ---")

        adata_sample = adata[adata.obs[sample_key] == sample].copy()
        sample_dir = module_dir / sample
        sample_dir.mkdir(exist_ok=True)

        if do_aggregation:
            cond = adata_sample.obs[condition_key].iloc[0]
            if cond not in agg_nhood:
                agg_nhood[cond] = []

        # =========================================================
        # FAST RESCUE MODE (Load from disk, skip heavy compute)
        # =========================================================
        if skip_compute:
            logger.info(f"Rescue Mode Active: Loading saved JSONs/CSVs for {sample}...")

            # Load Neighborhood Matrix and ALIGN SHAPES
            nhood_json = sample_dir / f"nhood_enrichment_{sample}.json"
            if nhood_json.exists() and do_aggregation:
                with open(nhood_json, "r") as f:
                    data = json.load(f)

                local_clusters = data["clusters"]
                local_matrix = np.array(data["zscores"])

                # Align to global shape
                aligned_matrix = np.full((num_clusters, num_clusters), np.nan)
                for i, c1 in enumerate(local_clusters):
                    for j, c2 in enumerate(local_clusters):
                        if c1 in valid_clusters and c2 in valid_clusters:
                            gi, gj = valid_clusters.index(c1), valid_clusters.index(c2)
                            aligned_matrix[gi, gj] = local_matrix[i, j]

                agg_nhood[cond].append(aligned_matrix)

            # Load Centrality
            cent_json = sample_dir / f"centrality_scores_{sample}.json"
            if cent_json.exists() and do_aggregation:
                cent_df = pd.read_json(cent_json, orient="index")
                cent_df["Sample"] = sample
                cent_df["Condition"] = cond
                cent_df["Cluster"] = cent_df.index
                agg_centrality.append(cent_df)

            # Load Moran's I
            moran_csv = sample_dir / f"moranI_results_{sample}.csv"
            if moran_csv.exists() and do_aggregation:
                moran_df = pd.read_csv(moran_csv, index_col=0)
                moran_df["Sample"] = sample
                moran_df["Condition"] = cond
                moran_df["Gene"] = moran_df.index
                agg_moran.append(moran_df)

            continue  # Skip the rest of the loop and go to the next sample

        # =========================================================
        # NORMAL MODE (Compute everything)
        # =========================================================
        logger.info("Building spatial neighborhood graph...")
        sq.gr.spatial_neighbors(adata_sample, coord_type="generic", delaunay=True)

        # CENTRALITY
        logger.info("Computing and plotting centrality scores...")
        sq.gr.centrality_scores(adata_sample, cluster_key=cluster_name)
        try:
            cent_df = adata_sample.uns[f"{cluster_name}_centrality_scores"]
            if do_aggregation:
                temp_cent = cent_df.copy()
                temp_cent["Sample"] = sample
                temp_cent["Condition"] = cond
                temp_cent["Cluster"] = temp_cent.index
                agg_centrality.append(temp_cent)

            cent_json_path = sample_dir / f"centrality_scores_{sample}.json"
            cent_df.replace([np.inf, -np.inf, np.nan], None).to_json(
                cent_json_path, orient="index"
            )
        except Exception as e:
            logger.warning(f"Failed to export Centrality JSON: {e}")

        with plt.rc_context({"figure.facecolor": "white", "axes.facecolor": "white"}):
            sq.pl.centrality_scores(
                adata_sample, cluster_key=cluster_name, figsize=cent_figsize
            )
            plt.savefig(
                sample_dir / f"centrality_scores_{sample}.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

        # CO-OCCURRENCE
        logger.info("Computing co-occurrence probability...")
        adata_subsample = sc.pp.subsample(adata_sample, fraction=0.5, copy=True)
        sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
        sq.gr.co_occurrence(
            adata_subsample, cluster_key=cluster_name, n_jobs=16, backend="loky"
        )
        valid_clusters_sub = list(adata_subsample.obs[cluster_name].cat.categories)

        try:
            co_occ_data = adata_subsample.uns[f"{cluster_name}_co_occurrence"]
            distances = co_occ_data["interval"].tolist()
            occ_matrix = co_occ_data["occ"]

            co_occ_export = {
                "clusters": valid_clusters_sub,
                "distances": distances,
                "probabilities": {},
            }
            for i, c1 in enumerate(valid_clusters_sub):
                for j, c2 in enumerate(valid_clusters_sub):
                    probs = np.nan_to_num(occ_matrix[i, j, :], nan=0.0).tolist()
                    co_occ_export["probabilities"][f"{c1}|{c2}"] = probs

            with open(sample_dir / f"co_occurrence_{sample}.json", "w") as f:
                json.dump(co_occ_export, f)
        except Exception as e:
            logger.warning(f"Failed to export Co-occurrence JSON: {e}")

        with plt.rc_context({"figure.facecolor": "white", "axes.facecolor": "white"}):
            sq.pl.co_occurrence(
                adata_subsample,
                cluster_key=cluster_name,
                clusters=valid_clusters_sub,
                figsize=co_figsize,
            )
            plt.savefig(
                sample_dir / f"co_occurrence_{sample}.png", dpi=300, bbox_inches="tight"
            )
            plt.close()

        # NEIGHBORHOOD ENRICHMENT
        logger.info("Performing neighborhood enrichment analysis...")
        sq.gr.nhood_enrichment(adata_sample, cluster_key=cluster_name)

        try:
            nhood_zscore = adata_sample.uns[f"{cluster_name}_nhood_enrichment"][
                "zscore"
            ]
            clusters = list(adata_sample.obs[cluster_name].cat.categories)

            if do_aggregation:
                aligned_matrix = np.full((num_clusters, num_clusters), np.nan)
                for i, c1 in enumerate(clusters):
                    for j, c2 in enumerate(clusters):
                        if c1 in valid_clusters and c2 in valid_clusters:
                            gi, gj = valid_clusters.index(c1), valid_clusters.index(c2)
                            aligned_matrix[gi, gj] = nhood_zscore[i, j]
                agg_nhood[cond].append(aligned_matrix)

            nhood_export = {
                "clusters": clusters,
                "zscores": np.nan_to_num(nhood_zscore, nan=0.0).tolist(),
            }
            with open(sample_dir / f"nhood_enrichment_{sample}.json", "w") as f:
                json.dump(nhood_export, f)
        except Exception as e:
            logger.warning(f"Failed to export Neighborhood Enrichment JSON: {e}")

        with plt.rc_context({"figure.facecolor": "white", "axes.facecolor": "white"}):
            fig, ax = plt.subplots(1, 1, figsize=nhood_figsize, facecolor="white")
            sq.pl.nhood_enrichment(
                adata_sample,
                cluster_key=cluster_name,
                title=f"Neighborhood enrichment ({sample})",
                ax=ax,
            )
            plt.savefig(
                sample_dir / f"nhood_enrichment_{sample}.png",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

        # MORAN'S I
        logger.info("Calculating Moran's I...")
        sq.gr.spatial_autocorr(
            adata_subsample, mode="moran", use_raw=True, n_perms=100, n_jobs=1
        )
        moran_df = adata_subsample.uns["moranI"]

        if do_aggregation:
            temp_moran = moran_df.copy()
            temp_moran["Sample"] = sample
            temp_moran["Condition"] = cond
            temp_moran["Gene"] = temp_moran.index
            agg_moran.append(temp_moran)

        moran_df.to_csv(sample_dir / f"moranI_results_{sample}.csv", index=True)
        logger.info(f"Finished sample {sample}.")

    # --- POST-LOOP CONDITION-LEVEL AGGREGATION & COMPARISONS ---
    if do_aggregation:
        logger.info("=== Running Condition-Level Aggregations ===")
        agg_dir = module_dir / "Aggregated_Results"
        agg_dir.mkdir(exist_ok=True)

        mean_nhoods = {}
        for condition, matrices in agg_nhood.items():
            if matrices:
                stacked = np.stack(matrices, axis=0)  # Now guaranteed to work!
                mean_matrix = np.nanmean(stacked, axis=0)
                mean_nhoods[condition] = mean_matrix

                plt.figure(figsize=nhood_figsize)
                sns.heatmap(
                    mean_matrix,
                    xticklabels=valid_clusters,
                    yticklabels=valid_clusters,
                    cmap="viridis",
                    center=0,
                    cbar_kws={"label": "Mean Z-score"},
                )
                plt.title(f"Consensus Neighborhood Enrichment: {condition}")
                plt.tight_layout()
                plt.savefig(
                    agg_dir / f"Consensus_Nhood_{condition}.png",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()

        if reference_condition and reference_condition in mean_nhoods:
            for condition, mean_matrix in mean_nhoods.items():
                if condition != reference_condition:
                    diff_matrix = mean_matrix - mean_nhoods[reference_condition]
                    plt.figure(figsize=nhood_figsize)
                    vmax = np.nanmax(np.abs(diff_matrix))
                    sns.heatmap(
                        diff_matrix,
                        xticklabels=valid_clusters,
                        yticklabels=valid_clusters,
                        cmap="coolwarm",
                        center=0,
                        vmin=-vmax,
                        vmax=vmax,
                        cbar_kws={
                            "label": f"Δ Z-score ({condition} - {reference_condition})"
                        },
                    )
                    plt.title(
                        f"Differential Neighborhood: {condition} vs {reference_condition}"
                    )
                    plt.tight_layout()
                    plt.savefig(
                        agg_dir
                        / f"Diff_Nhood_{condition}_vs_{reference_condition}.png",
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()

        if agg_centrality:
            df_cent = pd.concat(agg_centrality)
            if "degree_centrality" in df_cent.columns:
                plt.figure(figsize=(max(12, num_clusters), 6))
                sns.boxplot(
                    data=df_cent,
                    x="Cluster",
                    y="degree_centrality",
                    hue="Condition",
                    palette="Set2",
                )
                plt.xticks(rotation=45, ha="right")
                plt.title("Degree Centrality by Condition")
                plt.tight_layout()
                plt.savefig(
                    agg_dir / "Boxplot_Degree_Centrality.png",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()
                df_cent.to_csv(agg_dir / "Aggregated_Centrality.csv", index=False)

        if agg_moran:
            df_moran = pd.concat(agg_moran)
            top_genes = (
                df_moran.groupby("Gene")["I"]
                .mean()
                .sort_values(ascending=False)
                .head(20)
                .index
            )
            df_moran_top = df_moran[df_moran["Gene"].isin(top_genes)]

            plt.figure(figsize=(14, 6))
            sns.boxplot(
                data=df_moran_top, x="Gene", y="I", hue="Condition", palette="Set2"
            )
            plt.xticks(rotation=45, ha="right")
            plt.title("Moran's I (Spatial Clustering) of Top 20 Genes by Condition")
            plt.tight_layout()
            plt.savefig(agg_dir / "Boxplot_Morans_I.png", dpi=300, bbox_inches="tight")
            plt.close()
            df_moran.to_csv(agg_dir / "Aggregated_MoransI.csv", index=False)

    out_path = module_dir / input_adata_path.name
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
    logger.info("Spatial statistics module completed successfully.")
