"""cellphonedb module."""

import gc
import os
import warnings
from logging import getLogger
from pathlib import Path

import ktplotspy as kpy
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scipy.io as sio
import squidpy as sq
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from cellphonedb.utils import db_utils
from config import settings
from plotnine import theme

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def find_leiden_resolution(
    adata,
    target_clusters,
    microenv_key,
    microenv_nghbr_key,
    min_res=0.01,
    max_res=3.0,
    max_iters=30,
):
    """
    Performs a binary search to find the Leiden resolution that yields a specific number of clusters.
    """
    logger.info(f"Searching for resolution to yield {target_clusters} clusters...")

    for i in range(max_iters):
        # Calculate the midpoint resolution
        res = (min_res + max_res) / 2

        # Run Leiden with your specific parameters
        sc.tl.leiden(
            adata,
            resolution=res,
            neighbors_key=microenv_nghbr_key,
            key_added=microenv_key,
        )

        # Count the number of unique clusters generated
        current_clusters = adata.obs[microenv_key].nunique()

        logger.info(
            f"Iteration {i + 1}: Resolution = {res:.4f} -> Found {current_clusters} clusters"
        )

        # Check if we hit the target
        if current_clusters == target_clusters:
            logger.info(
                f"Success! Target {target_clusters} clusters reached at resolution {res:.4f}."
            )
            return res

        # Adjust search window
        elif current_clusters < target_clusters:
            # Too few clusters -> we need a higher resolution
            min_res = res
        else:
            # Too many clusters -> we need a lower resolution
            max_res = res

    logger.info(
        f"Max iterations reached. Settled on {current_clusters} clusters at resolution {res:.4f}."
    )


def microenvironment_split(adata, resolution, spatial_key, sample_key, module_dir):
    """
    User defines a leiden resolution to define spatial microenvironments for CellPhonedb and Decoupler
    analysis. Also produces a plot visualizing the split.
    """
    microenv_key = "spatial_microenvironment"

    # FIX: Use standard 'spatial_neighbors' key so Squidpy and Scanpy communicate perfectly
    microenv_nghbr_key = "spatial_neighbors"

    # Use scVIVA's niche composition to define pure spatial tissue regions
    if "niche_composition" in adata.obsm:
        logger.info(
            "Found 'niche_composition'! Clustering cells based purely on their spatial neighborhood to define tissue niches..."
        )

        sc.pp.neighbors(
            adata,
            use_rep="niche_composition",
            n_neighbors=30,
            key_added=microenv_nghbr_key,
        )
    else:
        # Fallback approach: Calculate nearest neighbors based strictly on physical space
        logger.info(
            "'niche_composition' not found. Building neighbourhood graph based on physical space (per sample)..."
        )
        sq.gr.spatial_neighbors(
            adata,
            coord_type="generic",
            n_neighs=30,
            spatial_key=spatial_key,
            library_key=sample_key,
            # Removed key_added here so Squidpy saves correctly to its default location
        )

    # run leiden clustering at a fixed resolution
    logger.info(
        f"Running leiden clustering at resolution {resolution} to define microenvironments..."
    )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key=microenv_nghbr_key,
        key_added=microenv_key,
    )

    num_clusters = adata.obs[microenv_key].nunique()
    logger.info(f"Generated {num_clusters} spatial microenvironments.")

    # remove color col from previous so the plotting doesnt crash
    adata.uns.pop(f"{microenv_key}_colors", None)

    logger.info("Visualizing microenvironments on tissue...")

    # loop over samples
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
        save_path = module_dir / f"microenvironment_res{resolution}_{sample}.png"
        fig.savefig(
            save_path,
            dpi=150,
            facecolor="white",
            bbox_inches="tight",
        )
        plt.close(fig)
        del adata_sample
        import gc

        gc.collect()

    logger.info(f"Saved microenvironments plots to {module_dir}")

    return microenv_key


def plot_cellphonedb(
    module_dir: str | Path,
    cpdb_counts_path: str | Path,
    celltype_key: str | None = None,
    celltypes: list[str] | None = None,
    gene_family: str | list[str] | None = None,
    target_microenv: list[str] | str | None = None,
    target_sample: list[str] | str | None = None,
    sample_key: str | None = None,
    microenv_key: str = "spatial_microenvironment",
    target_genes: list[str] | None = None,
    use_cellsign: bool = False,
):
    """
    cellphonedb plotting using ktplots-py
    """

    adata = sc.read_h5ad(cpdb_counts_path)
    cpbd_out_dir = Path(f"{module_dir}/cpdb_out")

    # If target genes are provided, drop all other genes from RAM immediately
    if target_genes:
        valid_genes = [g for g in target_genes if g in adata.var_names]
        adata = adata[:, valid_genes].copy()
        gc.collect()

    # 1. ALWAYS load the standard files to get the correct matrix shape, colors, and decon data
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
        logger.info(
            "use_cellsign=True: Mapping CellSign active TFs onto full p-values matrix..."
        )

        # Load the binary CellSign matrix
        cellsign_df = pd.read_csv(cellsign_files[0], sep="\t")

        # Identify the columns that represent cluster pairs
        pair_cols = [c for c in pvals.columns if "|" in c]

        # 1. Reset ALL p-values to 1.0 (Completely insignificant / Hidden)
        for c in pair_cols:
            pvals[c] = 1.0

        # 2. Safely map the active interactions onto the full matrix
        pvals.set_index("id_cp_interaction", drop=False, inplace=True)

        # Drop duplicates in CellSign just in case
        cellsign_df = cellsign_df.drop_duplicates(subset=["id_cp_interaction"])
        cellsign_df.set_index("id_cp_interaction", drop=False, inplace=True)

        common_rows = [i for i in cellsign_df.index if i in pvals.index]
        common_pairs = [c for c in cellsign_df.columns if c in pair_cols]

        # 3. Loop through ONLY the active rows and turn the specific dots "ON" (0.0)
        for row_id in common_rows:
            for col_name in common_pairs:
                val = cellsign_df.loc[row_id, col_name]
                if pd.to_numeric(val, errors="coerce") == 1:
                    pvals.loc[row_id, col_name] = 0.0

        pvals.reset_index(drop=True, inplace=True)
        del cellsign_df
        gc.collect()

        # 4. Set our plotting rules
        dot_sig_only = (
            False  # Keep the dot plot grid, but inactive ones will be tiny dots
        )
        chord_sig_only = True  # Strictly hide inactive ribbons in the chord plot

    else:
        if use_cellsign and not cellsign_files:
            logger.warning(
                "use_cellsign=True but no CellSign file found! Falling back to standard plots."
            )
        else:
            logger.info("Using standard results (use_cellsign=False).")

        dot_sig_only = False
        chord_sig_only = True

    # Clean dataframes
    def clean_df(df, fill_val):
        pair_cols = [c for c in df.columns if "|" in c]
        for c in pair_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(fill_val)
        return df

    means = clean_df(means, 1e-5)
    pvals = clean_df(pvals, 1.0)

    logger.info(f"Loaded Means from: {means_file.name}")
    logger.info(f"Loaded P-values from: {pvals_file.name}")
    logger.info(f"Loaded Deconvoluted from: {decon_file.name}")

    # --- Heatmaps ---
    logger.info("Generating symmetrical heatmap...")
    g1 = kpy.plot_cpdb_heatmap(pvals=pvals, figsize=(10, 10), title="Symmetrical")
    g1.savefig(
        str(module_dir / "cpdb_heatmap_symmetrical.png"), dpi=150, bbox_inches="tight"
    )
    plt.close("all")

    logger.info("Generating asymmetrical heatmap...")
    g2 = kpy.plot_cpdb_heatmap(
        pvals=pvals, figsize=(10, 10), title="Asymmetrical", symmetrical=False
    )
    g2.savefig(
        str(module_dir / "cpdb_heatmap_asymmetrical.png"), dpi=150, bbox_inches="tight"
    )
    plt.close("all")

    # ========================================================
    # Dot plot and chord diagram per celltype
    # ========================================================

    # --- Dot plot and chord diagram ---
    if celltypes:
        for celltype in celltypes:
            # Use target_sample OR target_microenv to subset the plotting interactions
            if (
                (target_sample or target_microenv)
                and microenv_key in adata.obs
                and celltype_key
            ):
                # --- NEW LOGIC: If a sample is selected ---
                if target_sample and sample_key and sample_key in adata.obs:
                    if isinstance(target_sample, str):
                        target_sample = [target_sample]

                    # 1. Grab only the cells in the targeted sample(s)
                    env_cells = adata.obs[adata.obs[sample_key].isin(target_sample)]

                    # 2. Automatically find out what microenvironments make up this sample (for logging)
                    auto_microenvs = env_cells[microenv_key].unique().tolist()
                    logger.info(
                        f"Auto-detected microenvironments for {target_sample}: {auto_microenvs}"
                    )

                    # 3. Naming variables safely based on sample name instead of 30+ microenvironments
                    env_str = "_".join([str(s) for s in target_sample])
                    plot_title = f"{celltype} interactions (Sample: {', '.join([str(s) for s in target_sample])})"

                # --- FALLBACK LOGIC: If manual microenvironments are set ---
                else:
                    if target_microenv and isinstance(target_microenv, str):
                        target_microenv = [target_microenv]

                    env_cells = adata.obs[adata.obs[microenv_key].isin(target_microenv)]

                    # Prevent huge filenames if many microenvs are selected
                    env_str_list = [str(m) for m in target_microenv]
                    if len(env_str_list) > 5:
                        env_str = "_".join(env_str_list[:5]) + "_etc"
                    else:
                        env_str = "_".join(env_str_list)

                    plot_title = f"{celltype} interactions (Microenvs: {', '.join(env_str_list)})"

                # Calculate neighbor cell types based on the isolated footprint
                ct_counts = env_cells[celltype_key].value_counts()
                neighbor_celltypes = ct_counts[ct_counts >= 50].index.tolist()

                if celltype in neighbor_celltypes:
                    neighbor_celltypes.remove(celltype)

                if not neighbor_celltypes:
                    logger.warning(
                        f"No other cell types found in selected region for {celltype}. Skipping."
                    )
                    continue

                logger.info(
                    f"Filtering to neighbors for plotting: {neighbor_celltypes}"
                )

                valid_pairs = [f"{celltype}|{celltype}"]
                for n in neighbor_celltypes:
                    valid_pairs.append(f"{celltype}|{n}")
                    valid_pairs.append(f"{n}|{celltype}")

                keep_cols = [
                    c for c in pvals.columns if ("|" not in c) or (c in valid_pairs)
                ]

                plot_means = means[keep_cols].copy()
                plot_pvals = pvals[keep_cols].copy()

            else:
                plot_means = means.copy()
                plot_pvals = pvals.copy()
                plot_title = f"{celltype} interactions (Global)"
                env_str = "Global"

            dot_title = plot_title

            # ================================================================
            # DOT PLOT
            # ================================================================
            try:
                p = kpy.plot_cpdb(
                    adata=adata,
                    cell_type1=celltype,
                    cell_type2=".",
                    means=plot_means,
                    pvals=plot_pvals,
                    celltype_key=celltype_key,
                    genes=target_genes,
                    title=dot_title,
                    keep_id_cp_interaction=True,
                    figsize=(19, 8),
                    default_style=False,
                    keep_significant_only=dot_sig_only,
                    alpha=0.05,
                ) + theme(legend_position="right")

                save_dir = module_dir / celltype
                save_dir.mkdir(exist_ok=True)

                outfile = save_dir / f"dotplot_{env_str}.png"

                p.save(
                    filename=str(outfile),
                    dpi=150,
                    limitsize=False,
                )
                del p
                gc.collect()

                logger.info(
                    f"Dot plot saved to {save_dir / f'dotplot_{env_str}.png'!s}"
                )

            except Exception as e:
                logger.error(f"Failed to generate Dotplot for {celltype}: {e}")

            # ========================================================
            # CHORD PLOT
            # ========================================================
            try:
                fig = plt.figure(figsize=(12, 10))
                kpy.plot_cpdb_chord(
                    adata=adata,
                    means=plot_means,
                    pvals=plot_pvals,
                    deconvoluted=decon,
                    celltype_key=celltype_key,
                    interaction=target_genes,
                    keep_significant_only=chord_sig_only,
                    link_kwargs={
                        "direction": 1,
                        "allow_twist": True,
                        "r1": 95,
                        "r2": 90,
                    },
                    sector_text_kwargs={
                        "color": "black",
                        "size": 12,
                        "r": 105,
                        "adjust_rotation": True,
                    },
                    link_offset=1,
                    legend_save_path=str(save_dir / f"chordplot_legend_{env_str}.png"),
                )

                fig = plt.gcf()
                fig.set_size_inches(12, 10)

                save_dir = module_dir / celltype
                save_dir.mkdir(exist_ok=True)

                fig.savefig(str(save_dir / f"chordplot_{env_str}.png"), dpi=150)
                plt.close(fig)
                del fig
                gc.collect()

                logger.info(
                    f"Chord Diagram saved to {save_dir / f'chordplot_{env_str}.png'!s}"
                )
            except Exception as e:
                logger.error(f"Failed to generate Chord plot for {celltype}: {e}")

            if "plot_means" in locals():
                del plot_means
            if "plot_pvals" in locals():
                del plot_pvals
            gc.collect()

    # ========================================================
    # GENE FAMILY PLOT
    # ========================================================
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
                keep_significant_only=False,  # Changed to False to prevent subsetting out rows
                default_style=False,
                alpha=1.01,
            )
            p.save(
                filename=str(module_dir / f"{gene_family}_interactions.png"),
                dpi=150,
                limitsize=False,
            )
        except Exception as e:
            logger.error(f"Failed to generate Gene Family plot: {e}")


def run_cellphonedb(
    module_dir,
    input_adata_path,
    sample_key,
    chosen_cluster,
    cpdb_version,
    active_tfs_file_path,
    human=True,
    microenv_resolution=0.5,
):
    """
    chosen_cluster could be clusters from leiden or celltypes from CellTypist / ScType
    """
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
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    adata.obs[chosen_cluster] = (
        adata.obs[chosen_cluster]
        .astype(str)
        .str.replace(r"[^a-zA-Z0-9 ]", "_", regex=True)
    )

    # setting spatial key
    data_type = settings["project"]["data_type"]
    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    microenv_key = microenvironment_split(
        adata, microenv_resolution, spatial_key, sample_key, module_dir
    )
    logger.info("Generating microenvironment mapping file...")
    # Get unique combinations of cell types and the microenvironments they live in
    microenv_mapping = adata.obs[[chosen_cluster, microenv_key]].drop_duplicates()
    microenv_mapping.columns = ["cell_type", "microenvironment"]
    microenvs_file_path = module_dir / "cpdb_microenvironments.tsv"
    microenv_mapping.to_csv(microenvs_file_path, sep="\t", index=False)

    # download database
    cpdb_target_dir = os.path.join("data/cellphonedb", cpdb_version)
    db_utils.download_database(cpdb_target_dir, cpdb_version)

    # create counts mtx file
    X = adata.layers["counts"]
    counts = X.T
    mtx_path = module_dir / "counts.mtx"
    sio.mmwrite(mtx_path, counts)

    del X, counts
    gc.collect()

    # if data not from human, assumes its from mouse
    # and convert mouse gene symbols to uppercase for human HGNC symbol
    if not human:
        adata.var_names = adata.var_names.str.upper()
        adata.var_names_make_unique()

    # normalize raw counts (this is what cellphonedb wants)
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]  # Free up the duplicate layer memory
    gc.collect()
    sc.pp.normalize_total(adata, target_sum=10000)
    adata.obsp.clear()
    adata.obsm.clear()
    adata.uns.clear()
    for layer in list(adata.layers.keys()):
        del adata.layers[layer]
    gc.collect()
    adata.write(cpdb_counts_path)

    # create metadata file which is a 2-column text file: 'Cell' and 'cell_type'
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
