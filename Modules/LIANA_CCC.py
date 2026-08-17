"""LIANA+ Single-Cell Spatial Cell-Cell Communication module."""

import itertools
import warnings
from logging import getLogger
from pathlib import Path

import liana as li
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_liana_pipeline(
    module_dir: Path,
    input_adata_path: Path,
    tf_adata_path: Path,
    spatial_key: str,
    sample_key: str,
    settings: dict,
):
    """Main runner for LIANA+ Spatial CCC."""
    module_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load Settings
    bandwidth = settings.get("bandwidth", 400)
    cutoff = settings.get("cutoff", 0.1)
    nz_prop = settings.get("nz_prop", 0.05)
    n_nmf_components = settings.get("n_nmf_components", 5)
    resource_name = settings.get(
        "resource_name", "consensus"
    )  # Use 'mouseconsensus' for mouse

    # 2. Load and Preprocess Main Spatial Data
    logger.info(f"Loading spatial data from {input_adata_path}...")
    adata = sc.read_h5ad(input_adata_path)

    # 2.5 Restore Raw Counts (Avoid taking log1p of negative scaled values)
    if "counts" in adata.layers:
        logger.info("Restoring raw counts from adata.layers['counts']...")
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        logger.info("Restoring data from adata.raw...")
        adata.X = adata.raw.X.copy()

    # Clean up any potential NaNs just to be bulletproof

    if sp.issparse(adata.X):
        adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Now it is safe to normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 3. Build Spatial Graph
    # logger.info(f"Building spatial graph (bandwidth={bandwidth}, cutoff={cutoff})...")
    # li.ut.spatial_neighbors(
    #     adata,
    #     bandwidth=bandwidth,
    #     cutoff=cutoff,
    #     kernel='gaussian',
    #     set_diag=True,
    #     spatial_key=spatial_key
    # )

    # 3. Build Spatial Graph (Sample-Aware)
    logger.info(
        f"Building spatial graph per sample (bandwidth={bandwidth}, cutoff={cutoff})..."
    )

    row_list, col_list, data_list = [], [], []

    # Iterate over each unique sample to avoid cross-sample "ghost" connections
    for sample in adata.obs[sample_key].unique():
        # 1. Find exact integer indices in the main adata for this sample
        mask = (adata.obs[sample_key] == sample).values
        idx = np.where(mask)[0]

        # 2. Create a temporary subset for the specific sample
        adata_sub = adata[mask].copy()

        # 3. Run LIANA's spatial neighbors safely on just this sample
        li.ut.spatial_neighbors(
            adata_sub,
            bandwidth=bandwidth,
            cutoff=cutoff,
            kernel="gaussian",
            set_diag=True,
            spatial_key=spatial_key,
        )

        # 4. Extract the subset connectivities and convert to coordinate (COO) format
        sub_conn = adata_sub.obsp["spatial_connectivities"].tocoo()

        # 5. Map the row/col indices from the subset back to the global adata indices
        mapped_rows = idx[sub_conn.row]
        mapped_cols = idx[sub_conn.col]

        row_list.append(mapped_rows)
        col_list.append(mapped_cols)
        data_list.append(sub_conn.data)

    # Combine all the localized matrices
    all_rows = np.concatenate(row_list)
    all_cols = np.concatenate(col_list)
    all_data = np.concatenate(data_list)

    # Build the global sparse matrix for the whole AnnData
    global_connectivities = sp.csr_matrix(
        (all_data, (all_rows, all_cols)), shape=(adata.n_obs, adata.n_obs)
    )

    # Store back into the main object
    adata.obsp["spatial_connectivities"] = global_connectivities
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "params": {
            "bandwidth": bandwidth,
            "cutoff": cutoff,
            "kernel": "gaussian",
            "set_diag": True,
        },
    }

    # 4. Run Bivariate Ligand-Receptor Analysis
    logger.info(f"Running Bivariate LR Analysis (resource={resource_name})...")
    lrdata = li.mt.bivariate(
        adata,
        resource_name=resource_name,
        local_name="cosine",
        global_name="morans",
        n_perms=100,
        mask_negatives=False,
        add_categories=True,
        nz_prop=nz_prop,
        use_raw=False,
        verbose=False,
    )

    lrdata_path = module_dir / "lrdata.h5ad"
    lrdata.write_h5ad(lrdata_path)
    logger.info(f"Saved LR data to {lrdata_path}")

    # Plot top 4 LR pairs by global Moran's R (Spatial autocorrelation)
    top_lrs = lrdata.var.sort_values("morans", ascending=False).head(4).index.tolist()
    if top_lrs:
        logger.info(f"Plotting top LR pairs: {top_lrs}")
        target_sample = adata.obs[sample_key].unique()[0]  # Grabs the first slide
        sq.pl.spatial_scatter(
            lrdata,
            spatial_key=spatial_key,
            library_key=sample_key,
            color=top_lrs,
            library_id=target_sample,
            shape=None,
            size=1,
            wspace=0.4,
            dpi=300,
            figsize=(10, 6),
        )
        plt.savefig(
            module_dir / "top_LR_spatial_patterns.png", dpi=300, bbox_inches="tight"
        )
        plt.close()

    # 5. Extract CCC Spatial Patterns via NMF
    logger.info(
        f"Extracting spatial CCC patterns via NMF (components={n_nmf_components})..."
    )
    li.multi.nmf(
        lrdata,
        n_components=n_nmf_components,
        inplace=True,
        random_state=0,
        max_iter=200,
        verbose=False,
    )

    lr_loadings = li.ut.get_variable_loadings(lrdata, varm_key="NMF_H").set_index(
        "index"
    )
    nmf_adata = sc.AnnData(
        X=lrdata.obsm["NMF_W"],
        obs=lrdata.obs,
        var=pd.DataFrame(index=lr_loadings.columns.astype(str)),
        uns=lrdata.uns,
        obsm=lrdata.obsm,
    )

    nmf_path = module_dir / "nmf_adata.h5ad"
    nmf_adata.write_h5ad(nmf_path)
    logger.info(f"Saved NMF patterns to {nmf_path}")

    # Plot NMF Factors
    target_sample = adata.obs[sample_key].unique()[0]
    sq.pl.spatial_scatter(
        nmf_adata,
        spatial_key=spatial_key,
        library_key=sample_key,
        color=[*nmf_adata.var.index],
        library_id=target_sample,
        shape=None,
        size=1,
        wspace=0.4,
        dpi=300,
        figsize=(10, 6),
    )
    plt.savefig(module_dir / "NMF_spatial_factors.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 6. Generalized Bivariate (Cross-talk between TFs and LR pairs)
    if tf_adata_path.exists():
        logger.info(
            "Module 7 TF data found. Running Generalized Bivariate (TF vs LR)..."
        )
        tf_adata = sc.read_h5ad(tf_adata_path)

        # Calculate Coefficient of Variation to find highly variable TFs
        tf_means = tf_adata.X.mean(axis=0)
        tf_stds = tf_adata.X.std(axis=0)
        # Avoid division by zero
        tf_means[tf_means == 0] = 1e-12
        cv = np.abs(tf_stds / tf_means)

        if isinstance(cv, np.matrix):
            cv = cv.A1  # flatten if matrix

        tf_adata.var["cv"] = cv

        # Select Top 10 TFs and Top 10 LRs to cross-correlate (prevents combinatorial explosion)
        top_tfs = (
            tf_adata.var.sort_values("cv", ascending=False).head(10).index.tolist()
        )
        interactions = list(
            itertools.product(top_tfs, top_lrs[:10])
        )  # Compare Top TFs to Top LRs

        # Find common cells across adata, tf_adata, and lrdata
        common_obs = adata.obs_names.intersection(tf_adata.obs_names).intersection(
            lrdata.obs_names
        )
        logger.info(
            f"Aligning adata, tf_adata, and lrdata to {len(common_obs)} common cells..."
        )

        # Subset all objects to the exact same valid cells
        tf_adata = tf_adata[common_obs].copy()
        lrdata_sub = lrdata[common_obs].copy()
        adata_sub = adata[common_obs].copy()

        # Build MuData object with perfectly aligned cells
        mdata = md.MuData({"tf": tf_adata, "lr": lrdata_sub})

        # Share spatial connectivities from the subsetted adata
        mdata.obsp["spatial_connectivities"] = adata_sub.obsp[
            "spatial_connectivities"
        ].copy()

        logger.info(
            f"Calculating spatial cross-correlation for {len(interactions)} TF-LR combinations..."
        )
        bdata = li.mt.bivariate(
            mdata,
            x_mod="tf",
            y_mod="lr",
            x_transform=sc.pp.scale,
            y_transform=sc.pp.scale,
            local_name="cosine",
            interactions=interactions,
            mask_negatives=True,
            add_categories=True,
            x_use_raw=False,
            y_use_raw=False,
            xy_sep="<->",
            x_name="tf",
            y_name="lr",
        )

        bdata_path = module_dir / "bdata_tf_lr.h5ad"
        bdata.write_h5ad(bdata_path)
        logger.info(f"Saved TF-LR cross-correlation to {bdata_path}")

        # Ensure sample metadata and spatial coordinates exist in bdata
        bdata.obs[sample_key] = adata_sub.obs[sample_key].copy()

        if spatial_key in adata_sub.obsm:
            bdata.obsm[spatial_key] = adata_sub.obsm[spatial_key].copy()
        elif "tf:centroid_x" in bdata.obs and "tf:centroid_y" in bdata.obs:
            bdata.obsm[spatial_key] = bdata.obs[
                ["tf:centroid_x", "tf:centroid_y"]
            ].values.astype(float)

        # Plot Top TF-LR correlations
        top_tf_lrs = (
            bdata.var.sort_values("mean", ascending=False).head(6).index.tolist()
        )
        if top_tf_lrs:
            target_sample = adata.obs[sample_key].unique()[0]

            # Determine correct library key for squidpy plotting
            lib_key = (
                f"tf:{sample_key}" if f"tf:{sample_key}" in bdata.obs else sample_key
            )

            sq.pl.spatial_scatter(
                bdata,
                spatial_key=spatial_key,
                library_key=lib_key,
                color=top_tf_lrs,
                library_id=target_sample,
                shape=None,
                size=1,
                wspace=0.4,
                dpi=300,
                figsize=(10, 6),
                cmap="coolwarm",
                vmax=1,
                vmin=-1,
            )
            plt.savefig(
                module_dir / "top_TF_LR_correlations.png", dpi=300, bbox_inches="tight"
            )
            plt.close()
    else:
        logger.warning(
            f"TF data not found at {tf_adata_path}. Skipping Generalized Bivariate analysis."
        )

    logger.info("LIANA+ Spatial CCC module completed successfully.")
