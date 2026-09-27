"""cellphonedb module."""

import gc
import os
import warnings
from logging import getLogger
from pathlib import Path
from typing import List, Optional, Union

import ktplotspy as kpy
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scipy.io as sio
import squidpy as sq
from anndata import AnnData
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from cellphonedb.utils import db_utils
from config import settings
from plotnine import theme

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def find_leiden_resolution(
    adata: AnnData,
    target_clusters: int,
    microenv_key: str,
    microenv_nghbr_key: str,
    min_res: float = 0.01,
    max_res: float = 3.0,
    max_iters: int = 30,
) -> float:
    """Performs a binary search to find the Leiden resolution that yields a specific number of clusters."""
    logger.info(f"Searching for resolution to yield {target_clusters} clusters...")

    res = 1.0
    for i in range(max_iters):
        res = (min_res + max_res) / 2
        sc.tl.leiden(
            adata,
            resolution=res,
            neighbors_key=microenv_nghbr_key,
            key_added=microenv_key,
        )
        current_clusters = adata.obs[microenv_key].nunique()

        logger.info(
            f"Iteration {i + 1}: Resolution = {res:.4f} -> Found {current_clusters} clusters"
        )

        if current_clusters == target_clusters:
            logger.info(
                f"Success! Target {target_clusters} clusters reached at resolution {res:.4f}."
            )
            return res
        elif current_clusters < target_clusters:
            min_res = res
        else:
            max_res = res

    logger.info(
        f"Max iterations reached. Settled on {current_clusters} clusters at resolution {res:.4f}."
    )
    return res


def microenvironment_split(
    adata: AnnData,
    resolution: float,
    spatial_key: str,
    sample_key: str,
    module_dir: Path,
) -> str:
    """Defines spatial microenvironments for CellPhonedb and Decoupler analysis."""
    microenv_key = "spatial_microenvironment"
    microenv_nghbr_key = "spatial_neighbors"

    if "niche_composition" in adata.obsm:
        logger.info(
            "Found 'niche_composition'! Clustering cells based purely on spatial neighborhood niches..."
        )
        sc.pp.neighbors(
            adata,
            use_rep="niche_composition",
            n_neighbors=30,
            key_added=microenv_nghbr_key,
        )
    else:
        logger.info(
            "'niche_composition' not found. Building neighbourhood graph based on physical space..."
        )
        sq.gr.spatial_neighbors(
            adata,
            coord_type="generic",
            n_neighs=30,
            spatial_key=spatial_key,
            library_key=sample_key,
        )

    logger.info(
        f"Running leiden clustering at resolution {resolution} to define microenvironments..."
    )
    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key=microenv_nghbr_key,
        key_added=microenv_key,
    )
    logger.info(
        f"Generated {adata.obs[microenv_key].nunique()} spatial microenvironments."
    )

    adata.uns.pop(f"{microenv_key}_colors", None)
    logger.info("Visualizing microenvironments on tissue...")

    for sample in adata.obs[sample_key].unique():
        adata_sample = adata[adata.obs[sample_key] == sample]
        fig, ax = plt.subplots(figsize=(10, 10), facecolor="white")
        sq.pl.spatial_scatter(
            adata_sample,
            color=microenv_key,
            spatial_key=spatial_key,
            shape=None,
            size=1,
            ax=ax,
            fig=fig,
        )

        plots_dir = module_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        save_path = plots_dir / f"microenvironment_res{resolution}_{sample}.png"
        fig.savefig(save_path, dpi=150, facecolor="white", bbox_inches="tight")
        plt.close(fig)

        del adata_sample
        gc.collect()

    logger.info(f"Saved microenvironments plots to {module_dir}")
    return microenv_key


def _apply_cellsign_mask(pvals: pd.DataFrame, cellsign_file: Path) -> pd.DataFrame:
    """Maps CellSign binary TF logic onto the CellphoneDB P-value matrix."""
    logger.info(
        "use_cellsign=True: Mapping CellSign active TFs onto full p-values matrix..."
    )
    cellsign_df = pd.read_csv(cellsign_file, sep="\t")
    pair_cols = [c for c in pvals.columns if "|" in c]

    # Reset ALL p-values to 1.0 (Completely insignificant / Hidden)
    for c in pair_cols:
        pvals[c] = 1.0

    pvals.set_index("id_cp_interaction", drop=False, inplace=True)
    cellsign_df = cellsign_df.drop_duplicates(subset=["id_cp_interaction"]).set_index(
        "id_cp_interaction", drop=False
    )

    common_rows = [i for i in cellsign_df.index if i in pvals.index]
    common_pairs = [c for c in cellsign_df.columns if c in pair_cols]

    for row_id in common_rows:
        for col_name in common_pairs:
            if pd.to_numeric(cellsign_df.loc[row_id, col_name], errors="coerce") == 1:
                pvals.loc[row_id, col_name] = 0.0

    pvals.reset_index(drop=True, inplace=True)
    return pvals


def _plot_dotplot(
    adata,
    celltype,
    means,
    pvals,
    celltype_key,
    target_genes,
    dot_sig_only,
    dot_title,
    env_str,
    save_dir,
):
    """Generates the ktplots-py Dotplot."""
    try:
        p = kpy.plot_cpdb(
            adata=adata,
            cell_type1=celltype,
            cell_type2=".",
            means=means,
            pvals=pvals,
            celltype_key=celltype_key,
            genes=target_genes,
            title=dot_title,
            keep_id_cp_interaction=True,
            figsize=(19, 8),
            default_style=False,
            keep_significant_only=dot_sig_only,
            alpha=0.05,
        ) + theme(legend_position="right")

        outfile = save_dir / f"{celltype}_dotplot_{env_str}.png"
        p.save(filename=str(outfile), dpi=150, limitsize=False)
        logger.info(f"Dot plot saved to {outfile}")
        del p
    except Exception as e:
        logger.error(f"Failed to generate Dotplot for {celltype}: {e}")


def _plot_chordplot(
    adata,
    celltype,
    means,
    pvals,
    decon,
    celltype_key,
    target_genes,
    chord_sig_only,
    env_str,
    save_dir,
):
    """Generates the ktplots-py Chord diagram."""
    try:
        fig = plt.figure(figsize=(12, 10))
        kpy.plot_cpdb_chord(
            adata=adata,
            means=means,
            pvals=pvals,
            deconvoluted=decon,
            celltype_key=celltype_key,
            interaction=target_genes,
            keep_significant_only=chord_sig_only,
            link_kwargs={"direction": 1, "allow_twist": True, "r1": 95, "r2": 90},
            sector_text_kwargs={
                "color": "black",
                "size": 12,
                "r": 105,
                "adjust_rotation": True,
            },
            link_offset=1,
            legend_save_path=str(
                save_dir / f"{celltype}_chordplot_legend_{env_str}.png"
            ),
        )
        fig = plt.gcf()
        fig.set_size_inches(12, 10)
        fig.savefig(str(save_dir / f"{celltype}_chordplot_{env_str}.png"), dpi=150)
        plt.close(fig)
        logger.info(
            f"Chord Diagram saved to {save_dir / f'{celltype}_chordplot_{env_str}.png'}"
        )
    except Exception as e:
        logger.error(f"Failed to generate Chord plot for {celltype}: {e}")


def plot_cellphonedb(
    module_dir: Union[str, Path],
    cpdb_counts_path: Union[str, Path],
    celltype_key: Optional[str] = None,
    celltypes: Optional[List[str]] = None,
    gene_family: Optional[Union[str, List[str]]] = None,
    target_microenv: Optional[Union[List[str], str]] = None,
    target_sample: Optional[Union[List[str], str]] = None,
    sample_key: Optional[str] = None,
    microenv_key: str = "spatial_microenvironment",
    target_genes: Optional[List[str]] = None,
    use_cellsign: bool = False,
):
    """cellphonedb plotting using ktplots-py"""
    module_dir = Path(module_dir)
    cpbd_out_dir = module_dir / "cpdb_out"
    adata = sc.read_h5ad(cpdb_counts_path)

    if target_genes:
        valid_genes = [g for g in target_genes if g in adata.var_names]
        adata = adata[:, valid_genes].copy()
        gc.collect()

    means_file = next(cpbd_out_dir.glob("*_means_*.txt"))
    pvals_file = next(cpbd_out_dir.glob("*_pvalues_*.txt"))
    decon_file = next(
        f
        for f in cpbd_out_dir.glob("*_deconvoluted_*.txt")
        if "percents" not in f.name and "CellSign" not in f.name
    )

    means = pd.read_csv(means_file, sep="\t")
    pvals = pd.read_csv(pvals_file, sep="\t")
    decon = pd.read_csv(decon_file, sep="\t")

    cellsign_files = list(cpbd_out_dir.glob("*CellSign_active_interactions_[0-9]*.txt"))

    if use_cellsign and cellsign_files:
        pvals = _apply_cellsign_mask(pvals, cellsign_files[0])
        dot_sig_only, chord_sig_only = False, True
    else:
        if use_cellsign:
            logger.warning(
                "use_cellsign=True but no CellSign file found! Falling back to standard plots."
            )
        logger.info("Using standard results (use_cellsign=False).")
        dot_sig_only, chord_sig_only = False, True

    def clean_df(df: pd.DataFrame, fill_val: float) -> pd.DataFrame:
        for c in [col for col in df.columns if "|" in col]:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(fill_val)
        return df

    means = clean_df(means, 1e-5)
    pvals = clean_df(pvals, 1.0)

    logger.info(f"Loaded Means from: {means_file.name}")
    logger.info(f"Loaded P-values from: {pvals_file.name}")
    logger.info(f"Loaded Deconvoluted from: {decon_file.name}")

    # Heatmaps
    hm_dir = module_dir / "heatmaps"
    hm_dir.mkdir(exist_ok=True)
    logger.info("Generating symmetrical heatmap...")
    g1 = kpy.plot_cpdb_heatmap(pvals=pvals, figsize=(10, 10), title="Symmetrical")
    g1.savefig(
        str(hm_dir / "cpdb_heatmap_symmetrical.png"), dpi=150, bbox_inches="tight"
    )

    logger.info("Generating asymmetrical heatmap...")
    g2 = kpy.plot_cpdb_heatmap(
        pvals=pvals, figsize=(10, 10), title="Asymmetrical", symmetrical=False
    )
    g2.savefig(
        str(hm_dir / "cpdb_heatmap_asymmetrical.png"), dpi=150, bbox_inches="tight"
    )
    plt.close("all")

    # Target Logic
    if celltypes:
        for celltype in celltypes:
            if (
                (target_sample or target_microenv)
                and microenv_key in adata.obs
                and celltype_key
            ):
                if target_sample and sample_key and sample_key in adata.obs:
                    target_sample = (
                        [target_sample]
                        if isinstance(target_sample, str)
                        else target_sample
                    )
                    env_cells = adata.obs[adata.obs[sample_key].isin(target_sample)]
                    logger.info(
                        f"Auto-detected microenvironments for {target_sample}: {env_cells[microenv_key].unique().tolist()}"
                    )
                    env_str = "_".join([str(s) for s in target_sample])
                    plot_title = f"{celltype} interactions (Sample: {', '.join([str(s) for s in target_sample])})"
                else:
                    target_microenv = (
                        [target_microenv]
                        if isinstance(target_microenv, str)
                        else target_microenv
                    )
                    env_cells = adata.obs[adata.obs[microenv_key].isin(target_microenv)]
                    env_str_list = [str(m) for m in target_microenv]
                    env_str = (
                        "_".join(env_str_list[:5]) + "_etc"
                        if len(env_str_list) > 5
                        else "_".join(env_str_list)
                    )
                    plot_title = f"{celltype} interactions (Microenvs: {', '.join(env_str_list)})"

                ct_counts = env_cells[celltype_key].value_counts()
                neighbor_celltypes = [
                    ct
                    for ct in ct_counts[ct_counts >= 50].index.tolist()
                    if ct != celltype
                ]

                if not neighbor_celltypes:
                    logger.warning(
                        f"No other cell types found in selected region for {celltype}. Skipping."
                    )
                    continue

                valid_pairs = (
                    [f"{celltype}|{celltype}"]
                    + [f"{celltype}|{n}" for n in neighbor_celltypes]
                    + [f"{n}|{celltype}" for n in neighbor_celltypes]
                )
                keep_cols = [
                    c for c in pvals.columns if ("|" not in c) or (c in valid_pairs)
                ]

                plot_means, plot_pvals = (
                    means[keep_cols].copy(),
                    pvals[keep_cols].copy(),
                )
            else:
                plot_means, plot_pvals = means.copy(), pvals.copy()
                plot_title, env_str = f"{celltype} interactions (Global)", "Global"

            dot_dir = module_dir / "dot_plots"
            dot_dir.mkdir(exist_ok=True)
            chord_dir = module_dir / "chord_plots"
            chord_dir.mkdir(exist_ok=True)

            _plot_dotplot(
                adata,
                celltype,
                plot_means,
                plot_pvals,
                celltype_key,
                target_genes,
                dot_sig_only,
                plot_title,
                env_str,
                dot_dir,
            )
            _plot_chordplot(
                adata,
                celltype,
                plot_means,
                plot_pvals,
                decon,
                celltype_key,
                target_genes,
                chord_sig_only,
                env_str,
                chord_dir,
            )
            gc.collect()

    if gene_family:
        try:
            p = kpy.plot_cpdb(
                adata=adata,
                cell_type1=".",
                cell_type2=".",
                means=means,
                pvals=pvals,
                celltype_key=celltype_key,
                gene_family=gene_family,
                highlight_size=1,
                figsize=(20, 8),
                keep_id_cp_interaction=True,
                keep_significant_only=False,
                default_style=False,
                alpha=1.01,
            )
            plots_dir = module_dir / "plots"
            plots_dir.mkdir(exist_ok=True)
            p.save(
                filename=str(plots_dir / f"{gene_family}_interactions.png"),
                dpi=150,
                limitsize=False,
            )
        except Exception as e:
            logger.error(f"Failed to generate Gene Family plot: {e}")


def run_cellphonedb(
    module_dir: Path,
    input_adata_path: Path,
    sample_key: str,
    chosen_cluster: str,
    cpdb_version: str,
    active_tfs_file_path: Optional[Path],
    human: bool = True,
    microenv_resolution: float = 0.5,
) -> Path:
    module_dir = Path(module_dir)
    cpdb_counts_path = module_dir / "adata.h5ad"
    cpdb_out = module_dir / "cpdb_out"

    if (
        cpdb_out.exists()
        and list(cpdb_out.glob("*_means_*.txt"))
        and cpdb_counts_path.exists()
    ):
        logger.info(
            f"CellPhoneDB results already exist in '{cpdb_out}'. Skipping preprocessing and execution."
        )
        return cpdb_counts_path

    os.makedirs(module_dir, exist_ok=True)
    adata = sc.read_h5ad(input_adata_path)

    adata.obs[chosen_cluster] = (
        adata.obs[chosen_cluster]
        .astype(str)
        .str.replace(r"[^a-zA-Z0-9 ]", "_", regex=True)
    )

    spatial_key = "global" if settings["project"]["data_type"] == "CosMx" else "spatial"
    microenv_key = microenvironment_split(
        adata, microenv_resolution, spatial_key, sample_key, module_dir
    )

    logger.info("Generating microenvironment mapping file...")
    microenv_mapping = adata.obs[[chosen_cluster, microenv_key]].drop_duplicates()
    microenv_mapping.columns = ["cell_type", "microenvironment"]
    microenvs_file_path = module_dir / "cpdb_microenvironments.tsv"
    microenv_mapping.to_csv(microenvs_file_path, sep="\t", index=False)

    cpdb_target_dir = os.path.join("data/cellphonedb", cpdb_version)
    db_utils.download_database(cpdb_target_dir, cpdb_version)

    X = adata.layers["counts"]
    mtx_path = module_dir / "counts.mtx"
    sio.mmwrite(mtx_path, X.T)
    del X
    gc.collect()

    if not human:
        adata.var_names = adata.var_names.str.upper()
        adata.var_names_make_unique()

    adata.X = adata.layers["counts"]
    del adata.layers["counts"]
    gc.collect()

    sc.pp.normalize_total(adata, target_sum=10000)
    adata.obsp.clear()
    adata.obsm.clear()
    adata.uns.clear()
    for layer in list(adata.layers.keys()):
        del adata.layers[layer]
    gc.collect()
    adata.write(cpdb_counts_path)

    meta = pd.DataFrame(
        {"Cell": adata.obs.index, "cell_type": adata.obs[chosen_cluster]}
    )
    cpdb_meta_path = module_dir / "cpdb_meta.tsv"
    meta.to_csv(cpdb_meta_path, sep="\t", index=False)
    logger.info("Data preparation for cellphonedb complete!")

    cpdb_path = os.path.join(cpdb_target_dir, "cellphonedb.zip")
    cpdb_out.mkdir(exist_ok=True)

    logger.info("Running CellPhoneDB...")
    _ = cpdb_statistical_analysis_method.call(
        cpdb_file_path=cpdb_path,
        meta_file_path=cpdb_meta_path,
        counts_file_path=cpdb_counts_path,
        counts_data="hgnc_symbol",
        microenvs_file_path=microenvs_file_path,
        active_tfs_file_path=str(active_tfs_file_path)
        if active_tfs_file_path
        else None,
        output_path=str(cpdb_out),
    )
    logger.info(f"Finished! Results are saved in the '{cpdb_out}' directory.")

    return cpdb_counts_path
