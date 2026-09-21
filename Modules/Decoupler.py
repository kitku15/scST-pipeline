"""Decoupler module."""

import json
import warnings
from logging import getLogger
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.cluster.hierarchy as sch
import squidpy as sq
from config import settings

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def make_tffile(module_dir, score, celltype_key, active_tfs_file_name, organism):
    # 1. Run statistical tests to find TFs active in each microenvironment
    # (This uses the 'score' AnnData object we made earlier, which holds the ULM estimates)
    marker_tfs = dc.tl.rankby_group(
        adata=score,
        groupby=celltype_key,
        reference="rest",
        method="t-test_overestim_var",
    )

    # 2. Filter for TFs that are strictly "Active"
    # We only want TFs that are positively enriched (stat > 0) and statistically significant (padj < 0.05)
    active_tfs_filtered = marker_tfs[
        (marker_tfs["stat"] > 0) & (marker_tfs["padj"] < 0.05)
    ].copy()

    # 3. Format the DataFrame for CellphoneDB
    # Rename the columns to 'cluster' and 'TF'
    active_tfs_filtered = active_tfs_filtered.rename(
        columns={"group": "cell_type", "name": "TF"}
    )

    # Keep only the two required columns
    final_active_tf_df = active_tfs_filtered[["cell_type", "TF"]].copy()

    # make cell types in the TF file to match the format used by CellPhoneDB.py
    final_active_tf_df["cell_type"] = (
        final_active_tf_df["cell_type"]
        .astype(str)
        .str.replace(r"[^a-zA-Z0-9 ]", "_", regex=True)
    )

    if organism == "mouse":
        final_active_tf_df["TF"] = final_active_tf_df["TF"].astype(str).str.upper()

    # 4. Save to a tab-separated text file
    save_path = module_dir / active_tfs_file_name
    final_active_tf_df.to_csv(save_path, sep="\t", index=False)

    logger.info(
        f"Exported {len(final_active_tf_df)} active TF-cluster pairs to {active_tfs_file_name}."
    )


def tf_enrichment(
    module_dir,
    web_dir,
    input_adata_path,
    sample_key,
    celltype_key,
    organism,
    grn,
    dorothea_levels=["A", "B"],
    active_tfs_file_name=None,
):
    """
    Transcription factor enrichment analysis using user's choice of gene regulatory network (GRN)
    could choose from:
        - collectri (CollecTRI)
        - dorothea (DoRothEA)
        - (CellOracle) not sure how to do this yet
        - (pySCENIC) not sure how to do this yet
        - (SCENIC+)' not sure how to do this yet

    user needs to set which organism their data is from
    could choose from:
        - human
        - mouse
        - rat
    """

    module_dir.mkdir(exist_ok=True)
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    # setting spatial key
    data_type = settings["project"]["data_type"]
    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    def grn_selection():
        # setting grn
        grn_list = ["collectri", "dorothea"]

        if grn not in grn_list:
            # user selected model thats unavailable
            message = (
                f"Invalid model selection: {grn}\n"
                "Please select one of the following:\n- " + "\n- ".join(grn_list)
            )

            logger.error(message)
            raise ValueError(message)
        elif grn == "collectri":
            collectri = dc.op.collectri(organism=organism)
            logger.info("Using collectri network")
            return collectri
        elif grn == "dorothea":
            dorothea = dc.op.dorothea(organism=organism, levels=dorothea_levels)
            logger.info("Using dorothea network")
            return dorothea

    network = grn_selection()

    # 2. Match gene casing
    if organism == "human":
        adata.var_names = adata.var_names.str.upper()
    elif organism == "mouse":
        adata.var_names = adata.var_names.str.capitalize()

    # 3. RUN ULM
    result = dc.mt.ulm(data=adata, net=network, raw=True, verbose=True)

    # If Decoupler dropped empty cells, it returns a new AnnData object. Catch it!
    if isinstance(result, sc.AnnData):
        adata = result
        logger.info(
            "Decoupler dropped empty cells and returned a repaired AnnData object."
        )

    # extract scores as new anndata object
    score = dc.pp.get_obsm(adata=adata, key="score_ulm")

    # 1. PURGE UNUSED CATEGORIES & CACHED DENDROGRAMS
    if hasattr(score.obs[celltype_key], "cat"):
        score.obs[celltype_key] = score.obs[celltype_key].cat.remove_unused_categories()

    dendro_key = f"dendrogram_{celltype_key}"
    if dendro_key in score.uns:
        del score.uns[dendro_key]

    # Extract unique TFs from the 'source' column
    tf_list = network["source"].unique().tolist()
    logger.info(f"Total number of transcription factors: {len(tf_list)}")

    # Filter out rare cell types (< 2 cells) so the t-test doesn't divide by zero
    val_counts = score.obs[celltype_key].value_counts()
    valid_groups = val_counts[val_counts >= 2].index.tolist()
    
    if len(valid_groups) == 0:
        logger.error(f"No groups with >= 2 cells found in {celltype_key}. Cannot run TF enrichment.")
        return
        
    score_filtered = score[score.obs[celltype_key].isin(valid_groups)].copy()
    score_filtered.obs[celltype_key] = score_filtered.obs[celltype_key].cat.remove_unused_categories()

    # identifying marker TFs for each spatial microenvironment
    df = dc.tl.rankby_group(
        adata=score_filtered,
        groupby=celltype_key,
        reference="rest",
        method="t-test_overestim_var",
    )

    # Get all active categories present in score.obs[celltype_key]
    if hasattr(score.obs[celltype_key], "cat"):
        all_categories = list(score.obs[celltype_key].cat.categories)
    else:
        all_categories = list(score.obs[celltype_key].unique())

    # Find top 3 TF markers per microenvironment while guaranteeing ALL categories exist in source_markers
    n_markers = 3
    source_markers = {}

    for cat in all_categories:
        df_cat = df[df["group"] == cat]
        top_tfs = []
        if not df_cat.empty:
            df_pos = df_cat[df_cat["stat"] > 0.0]
            if not df_pos.empty:
                top_tfs = (
                    df_pos.sort_values("stat", ascending=False)["name"]
                    .unique()[:n_markers]
                    .tolist()
                )
            else:
                top_tfs = (
                    df_cat.sort_values("stat", ascending=False)["name"]
                    .unique()[:n_markers]
                    .tolist()
                )

        # 2. FALLBACK: If a category has no specific TFs, assign top overall TFs so Scanpy never crashes
        if not top_tfs:
            top_tfs = (
                df.sort_values("stat", ascending=False)["name"]
                .unique()[:n_markers]
                .tolist()
            )

        source_markers[cat] = top_tfs

    # check how many of the marker TFs I actually have in my RNA panel
    all_markers = []
    for tfs in source_markers.values():
        all_markers.extend(tfs)
    present_tfs = [tf for tf in all_markers if tf in adata.var_names]
    logger.info(
        f"Out of {len(all_markers)} transcription factors checked, you have {len(present_tfs)} in your RNA panel.\n"
    )

    logger.info("Generating TF heatmap...")

    sc.pl.matrixplot(
        adata=score,
        var_names=source_markers,
        groupby=celltype_key,
        dendrogram=True,
        standard_scale="var",
        colorbar_title="Z-scaled scores",
        cmap="RdBu_r",
        title="Top 3 TFs per microenvironment",
        figsize=(10, 10),
        show=False,
    )

    save_path = module_dir / "top3tfs_heatmap.png"

    # Grab the current figure generated by Scanpy and save it
    plt.savefig(
        save_path,
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close()

    logger.info(f"Saved tf heatmap plot to {save_path}")

    def export_json_heatmap():
        # 1. Flatten top TFs into a single list
        all_heatmap_tfs = []
        for tfs in source_markers.values():
            for tf in tfs:
                if tf not in all_heatmap_tfs:
                    all_heatmap_tfs.append(tf)

        # 2. Get mean scores for all TFs per cell type
        mean_scores = score[:, all_heatmap_tfs].to_df()
        mean_scores[celltype_key] = score.obs[celltype_key].values
        mean_scores_grouped = mean_scores.groupby(celltype_key).mean()

        # 3. standard_scale="var": Scale each column (TF) individually from 0 to 1
        denom = mean_scores_grouped.max() - mean_scores_grouped.min()
        denom[denom == 0] = 1.0  # Prevent division by zero
        scaled_df = (mean_scores_grouped - mean_scores_grouped.min()) / denom

        # 4. dendrogram=True: Perform Hierarchical Clustering ONLY on the cell types (Rows)
        row_linkage = sch.linkage(scaled_df, method="ward")
        row_order = sch.leaves_list(row_linkage)

        # Get the newly sorted cell types
        ordered_celltypes = scaled_df.index[row_order].tolist()

        # 5. var_names=source_markers: Reorder the X-axis to match the Y-axis clustered order!
        ordered_tfs = []
        for ct in ordered_celltypes:
            if ct in source_markers:
                for tf in source_markers[ct]:
                    if tf not in ordered_tfs:
                        ordered_tfs.append(tf)

        # 6. Apply both ordered lists to the dataframe
        final_df = scaled_df.loc[ordered_celltypes, ordered_tfs]

        # 7. Export for React Plotly Heatmap
        heatmap_data = {
            "x": final_df.columns.tolist(),  # TFs ordered sequentially by their clustered cell type
            "y": final_df.index.tolist(),  # Clustered Cell Types
            "z": final_df.values.tolist(),  # Scaled scores matrix
        }

        web_dir.mkdir(exist_ok=True)
        aux_dir = web_dir / "aux_data"
        aux_dir.mkdir(exist_ok=True)
        with open(aux_dir / "tf_heatmap_data.json", "w") as f:
            json.dump(heatmap_data, f, indent=4)

    # export to json for webtool
    export_json_heatmap()

    # visualize top marker TFs for each microenvironment on a spatial plot
    top_tfs = []
    top_tfs_microenvs = []

    # get only the top TF for each microenv
    for microenv, tfs in source_markers.items():
        if tfs:
            top_tfs.append(tfs[0])
            top_tfs_microenvs.append(microenv)

    logger.info("Visualizing top tfs per microenvironment on tissue...")

    # LOOP OVER SAMPLES TO PREVENT SQUIDPY GRID CRASH
    for sample in score.obs[sample_key].unique():
        score_sample = score[score.obs[sample_key] == sample].copy()

        # We generate custom titles for this specific sample
        titles = [
            f"{t} (Microenv {m}) - {sample}" for t, m in zip(top_tfs, top_tfs_microenvs)
        ]

        sq.pl.spatial_scatter(
            score_sample,
            color=top_tfs,
            spatial_key=spatial_key,
            shape=None,
            size=1,
            legend_loc="right margin",
            cmap="RdBu_r",
            figsize=(5, 5),
            dpi=300,
            title=titles,
        )

        save_path = module_dir / f"microenvironment_tfs_{sample}.png"
        plt.savefig(save_path, dpi=300, facecolor="white", bbox_inches="tight")
        plt.close()

    logger.info(f"Saved microenvironment TF plots to {module_dir}")

    # make TF activity file for cellphonedb
    make_tffile(module_dir, score, celltype_key, active_tfs_file_name, organism)

    # generate pseudobulk adata for DE analysis
    logger.info("Generating pseudobulk AnnData...")
    try:
        # 1. Create the raw pseudobulk
        psbulk = dc.pp.pseudobulk(
            adata,
            sample_col=sample_key,
            groups_col=celltype_key,
            layer="counts",
            mode="sum",
        )

        # 2. BULLETPROOF: Manually calculate the exact total counts per pseudobulk sample
        psbulk.obs["manual_total_counts"] = np.array(psbulk.X.sum(axis=1)).flatten()

        # 3. BULLETPROOF: Manually calculate the exact number of cells using the original adata
        cell_counts_dict = (
            adata.obs.groupby([sample_key, celltype_key]).size().to_dict()
        )

        psbulk_n_cells = []
        for _, row in psbulk.obs.iterrows():
            # Get the exact sample and cell type for this pseudobulk row
            s_val = row[sample_key]
            c_val = row[celltype_key]
            psbulk_n_cells.append(cell_counts_dict.get((s_val, c_val), 0))

        psbulk.obs["manual_n_cells"] = psbulk_n_cells

        # 4. Apply the thresholds safely
        valid_mask = (psbulk.obs["manual_n_cells"] >= 10) & (
            psbulk.obs["manual_total_counts"] >= 1000
        )
        psbulk = psbulk[valid_mask].copy()

        # 5. Filter unexpressed genes
        try:
            dc.pp.filter_by_prop(psbulk, min_prop=0.2, min_smpls=2)
        except AttributeError:
            # Fallback if Decoupler version is too old for filter_by_prop
            sc.pp.filter_genes(psbulk, min_cells=2)

        psbulk_path = module_dir / "pseudobulk.h5ad"
        psbulk.write_h5ad(psbulk_path)
        logger.info(
            f"Pseudobulk data saved to {psbulk_path} ({psbulk.n_obs} valid pseudobulk samples remaining)"
        )

    except Exception as e:
        logger.error(f"Failed to generate pseudobulk data: {e}")

    # save adata
    out_path = module_dir / input_adata_path.name
    adata.write_h5ad(out_path)
    score.write_h5ad(module_dir / "tf_activity_scores.h5ad")
    logger.info(f"Data saved to {out_path}")
