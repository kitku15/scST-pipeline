import re
import warnings
from logging import getLogger
from pathlib import Path

import corneto as cn
import decoupler as dc
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import omnipath as op
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def sanitize_name(name: str) -> str:
    """Helper to ensure cell type and condition names are safe for filenames."""
    return re.sub(r"[^\w\s-]", "", str(name)).replace(" ", "_")


def run_pseudobulk_dea_with_shrinkage(
    pdata: sc.AnnData,
    groupby: str,
    condition_key: str,
    control_group: str,
    stim_group: str,
    min_count_expr: int = 5,
    min_total_count_expr: int = 10,
    quiet: bool = True,
) -> pd.DataFrame:
    """Runs PyDESeq2 DEA with LFC shrinkage on clean, picklable pseudobulk data."""
    logger.info("Running PyDESeq2 with LFC shrinkage...")
    pdata.obs[condition_key] = pdata.obs[condition_key].astype(str)

    dea_results = {}
    unique_groups = pdata.obs[groupby].unique()

    for cell_group in unique_groups:
        ctdata = pdata[pdata.obs[groupby] == cell_group].copy()
        observed_conditions = ctdata.obs[condition_key].unique()

        if (
            control_group not in observed_conditions
            or stim_group not in observed_conditions
        ):
            logger.warning(
                f"Skipping group '{cell_group}': missing control or stim conditions."
            )
            continue

        n_ctrl = (ctdata.obs[condition_key] == control_group).sum()
        n_stim = (ctdata.obs[condition_key] == stim_group).sum()
        if n_ctrl < 2 or n_stim < 2:
            logger.warning(
                f"Skipping group '{cell_group}': insufficient replicates (Control: {n_ctrl}, Stim: {n_stim})."
            )
            continue

        genes = dc.pp.filter_by_expr(
            ctdata,
            group=condition_key,
            min_count=min_count_expr,
            min_total_count=min_total_count_expr,
            inplace=False,
        )

        if genes is None or len(genes) == 0 or not np.any(genes):
            continue

        ctdata = ctdata[:, genes].copy()

        # Ensure categories are set properly
        ctdata.obs[condition_key] = pd.Categorical(
            ctdata.obs[condition_key],
            categories=[control_group]
            + [c for c in observed_conditions if c != control_group],
            ordered=False,
        )

        # Reconstruct a clean AnnData object
        # Strip away any nested sparse view matrices or unpicklable metadata
        counts = (
            ctdata.X.toarray() if hasattr(ctdata.X, "toarray") else np.asarray(ctdata.X)
        )
        counts = np.round(counts).astype(int)

        clean_ctdata = sc.AnnData(
            X=counts,
            obs=pd.DataFrame(ctdata.obs.copy()),
            var=pd.DataFrame(index=ctdata.var_names),
        )

        try:
            dds = DeseqDataSet(
                adata=clean_ctdata,
                design_factors=condition_key,
                refit_cooks=True,
                n_cpus=1,
                quiet=quiet,
            )
            dds.deseq2()

            stat_res = DeseqStats(
                dds,
                contrast=[condition_key, stim_group, control_group],
                n_cpus=1,
                quiet=quiet,
            )
            stat_res.summary()

            coeff_name = None
            for c in dds.varm["LFC"].columns:
                if stim_group in c:
                    coeff_name = c
                    break

            if coeff_name is not None:
                # lfc_shrink inherits the single-threaded context
                stat_res.lfc_shrink(coeff=coeff_name)

            dea_results[cell_group] = stat_res.results_df
            logger.info(f"PyDESeq2 and LFC shrinkage completed for {cell_group}")

        except Exception as e:
            logger.error(f"PyDESeq2 run failed for {cell_group}: {e}")

    if not dea_results:
        raise ValueError(
            "No cell groups had sufficient data to complete PyDESeq2 analysis."
        )

    dea_df = pd.concat(dea_results)
    dea_df = (
        dea_df.reset_index()
        .rename(columns={"level_0": groupby, "level_1": "index"})
        .set_index("index")
    )
    return dea_df


def run_condition_ccc_pipeline(
    module_dir: Path,
    input_adata_path: Path,
    pseudobulk_adata_path: Path,
    settings: dict,
):
    """Executes condition-specific CCC, TF activity inference, and causal modeling."""
    module_dir.mkdir(parents=True, exist_ok=True)

    # Load configuration parameters
    # sample_key = settings.get("sample_key")
    groupby = settings.get("celltype_col")
    condition_key = settings.get("treatment_col")
    comparisons = settings.get("comparisons", [])
    cell_type_pairs = settings.get("cell_type_pairs", [])
    organism = settings.get("organism", "human")
    solver = settings.get("solver", "scipy")

    if not comparisons:
        logger.warning("Comparisons missing from config. Skipping pipeline entirely.")
        return

    if not cell_type_pairs:
        logger.info(
            "No cell_type_pairs provided. Will run DEA, LIANA, and TF inference, but skip causal network inference."
        )

    logger.info("Loading preprocessed single-cell data and pseudobulk data...")
    adata = sc.read_h5ad(input_adata_path)
    pdata = sc.read_h5ad(pseudobulk_adata_path)

    # 1. Outer Loop: Process one Condition Comparison at a time sequentially
    for stim_group, control_group in comparisons:
        comp_name = f"{stim_group}_vs_{control_group}"
        safe_comp_name = sanitize_name(comp_name)

        comp_dir = module_dir / safe_comp_name
        comp_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            f"\n{'=' * 50}\nProcessing Condition: {comp_name}\nOutputs routing to: {comp_dir.name}/\n{'=' * 50}"
        )

        try:
            dea_df = run_pseudobulk_dea_with_shrinkage(
                pdata=pdata,
                groupby=groupby,
                condition_key=condition_key,
                control_group=control_group,
                stim_group=stim_group,
                quiet=True,
            )
            dea_df.to_csv(comp_dir / f"dea_{safe_comp_name}.csv")
        except ValueError as e:
            logger.warning(f"Skipping comparison {comp_name}: {e}")
            continue

        logger.info(f"Mapping condition-specific LR interactions for {comp_name}...")
        adata_sub = adata[adata.obs[condition_key] == stim_group].copy()

        if "counts" in adata_sub.layers:
            adata_sub.X = adata_sub.layers["counts"].copy()
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)

        stat_keys = ["stat", "pvalue", "padj"]
        lr_res = li.multi.df_to_lr(
            adata_sub,
            dea_df=dea_df,
            resource_name="consensus" if organism == "human" else "mouseconsensus",
            expr_prop=settings.get("expr_prop", 0.1),
            groupby=groupby,
            stat_keys=stat_keys,
            use_raw=False,
            verbose=False,
        )
        lr_res.to_csv(comp_dir / f"liana_lr_{safe_comp_name}.csv", index=False)

        logger.info("Generating diagnostic plots...")
        _save_plots(comp_dir, lr_res)

        logger.info(f"Inferring TF activities for {comp_name}...")
        net = dc.op.collectri(organism=organism, remove_complexes=False, verbose=False)
        dea_wide = (
            dea_df[[groupby, "stat"]]
            .reset_index(names="genes")
            .pivot(index=groupby, columns="genes", values="stat")
            .fillna(0)
        )
        tf_estimates, _ = dc.mt.ulm(data=dea_wide, net=net)
        tf_estimates.to_csv(comp_dir / f"tf_estimates_{safe_comp_name}.csv")

        if not cell_type_pairs:
            logger.info(
                f"Skipping causal network inference for {comp_name} (no cell_type_pairs provided)."
            )
            continue

        if cn is None:
            logger.warning(
                "Corneto is not installed. Skipping causal network inference."
            )
            continue

        ppis = op.interactions.OmniPath().get(genesymbols=True, organism=organism)
        ppis["mor"] = ppis["is_stimulation"].astype(int) - ppis["is_inhibition"].astype(
            int
        )
        filtered_ppis = ppis[
            (ppis["mor"] != 0)
            & (ppis["curation_effort"] >= settings.get("curation_effort", 5))
            & ppis["consensus_direction"]
        ]

        if filtered_ppis.empty:
            logger.warning("No OmniPath interactions met the curation threshold.")
            continue

        input_pkn = filtered_ppis[["source_genesymbol", "mor", "target_genesymbol"]]
        input_pkn.columns = ["source", "mor", "target"]

        # 2. Inner Loop: Process each cell type pair SEQUENTIALLY
        for source_ct, target_ct in cell_type_pairs:
            pair_name = f"{source_ct}_to_{target_ct}"
            safe_pair_name = sanitize_name(pair_name)
            logger.info(f"  --> Inferring Causal Network: {source_ct} -> {target_ct}")

            lr_stats = lr_res[
                lr_res["source"].isin([source_ct]) & lr_res["target"].isin([target_ct])
            ].copy()
            lr_stats = lr_stats.sort_values(
                "interaction_stat", ascending=False, key=abs
            )

            # Drop NaNs safely before dictionary conversion
            lr_stats_clean = lr_stats.dropna(subset=["interaction_stat"])
            lr_dict = lr_stats_clean.set_index("receptor")["interaction_stat"].to_dict()
            input_scores = {
                k: v
                for i, (k, v) in enumerate(
                    sorted(lr_dict.items(), key=lambda item: abs(item[1]), reverse=True)
                )
                if i < settings.get("num_top_receptors", 10)
            }

            # Extract TFs, drop NaNs, convert to dict
            if target_ct in tf_estimates.index:
                tf_dict = tf_estimates.loc[target_ct].dropna().to_dict()
            else:
                tf_dict = {}

            output_scores = {
                k: v
                for i, (k, v) in enumerate(
                    sorted(tf_dict.items(), key=lambda item: abs(item[1]), reverse=True)
                )
                if i < settings.get("num_top_tfs", 5)
            }

            if not input_scores or not output_scores:
                logger.info(
                    f"      [Skipped] Insufficient LR or TF scores for {pair_name}"
                )
                continue

            try:
                prior_graph = li.mt.build_prior_network(
                    input_pkn, input_scores, output_scores, verbose=False
                )

                # Pre-calculate cell-type specific node weights
                temp_ct = adata[adata.obs[groupby] == target_ct].copy()
                if temp_ct.n_obs > 0:
                    matrix = (
                        temp_ct.layers["counts"]
                        if "counts" in temp_ct.layers
                        else temp_ct.X
                    )
                    non_zeros = (
                        matrix.getnnz(axis=0)
                        if hasattr(matrix, "getnnz")
                        else np.count_nonzero(matrix, axis=0)
                    )
                    node_weights = dict(
                        zip(temp_ct.var_names, non_zeros / temp_ct.n_obs)
                    )
                else:
                    logger.warning(
                        f"      [Warning] 0 cells found for {target_ct}. Defaulting weights to 0."
                    )
                    node_weights = {gene: 0.0 for gene in adata.var_names}

                out_csv = comp_dir / f"causal_net_{safe_comp_name}_{safe_pair_name}.csv"

                causal_net_res, _ = li.mt.find_causalnet(
                    prior_graph,
                    input_scores,
                    output_scores,
                    node_weights,
                    node_cutoff=settings.get("node_cutoff", 0.1),
                    max_penalty=1.0,
                    min_penalty=0.01,
                    edge_penalty=0.1,
                    verbose=False,
                    max_runs=settings.get("max_runs", 50),
                    stable_runs=settings.get("stable_runs", 10),
                    solver=solver,
                )
                causal_net_res = causal_net_res.drop_duplicates(
                    subset=["source", "target"]
                )
                causal_net_res.to_csv(out_csv, index=False)
                logger.info(f"      [Success] Saved to {out_csv.name}")

            except Exception as e:
                logger.error(
                    f"      [Failed] Network inference failed for {pair_name}: {e}"
                )
                continue


def _save_plots(output_dir: Path, lr_res: pd.DataFrame):
    """Generates standard diagnostic plots for condition-specific ligand-receptor interactions."""
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        lr_res["interaction_stat"].hist(bins=50, ax=ax)
        ax.set_title("Distribution of Interaction Statistics")
        ax.set_xlabel("Interaction Stat (Wald)")
        ax.set_ylabel("Count")
        plt.tight_layout()
        plt.savefig(output_dir / "interaction_stat_histogram.png", dpi=300)
        plt.close()
    except Exception as e:
        logger.warning(f"Could not generate interaction histogram: {e}")

    try:
        tile_gg = li.pl.tileplot(
            liana_res=lr_res,
            fill="expr",
            label="padj",
            label_fun=lambda x: "*" if x < 0.05 else np.nan,
            top_n=15,
            orderby="interaction_stat",
            orderby_ascending=False,
            orderby_absolute=False,
            source_title="Ligand",
            target_title="Receptor",
        )
        tile_gg.save(str(output_dir / "liana_tileplot.png"), dpi=300, width=8, height=8)
    except Exception as e:
        logger.warning(f"Could not generate LIANA tileplot: {e}")
