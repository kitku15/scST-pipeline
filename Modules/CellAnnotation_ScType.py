"""Cell Annotation Functions for ScType called in Annotation module."""

import warnings
from logging import getLogger

import pandas as pd
import numpy as np
from sctype_py import sctype_score

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def gene_sets_prepare(path_to_db_file, tissue_type=None):
    """
    Function that gets marker genes from ScType database or
    user's own marker list
    """
    if tissue_type:
        logger.info(f"Loading {tissue_type} markers")

    df = pd.read_excel(path_to_db_file)

    if tissue_type:  # filter by tissue type if provided
        df = df[df["tissueType"] == tissue_type]

    gs_positive, gs_negative = {}, {}
    for _, row in df.iterrows():
        cell_name = row["cellName"]
        pos_genes = str(row["geneSymbolmore1"]).replace(" ", "")
        gs_positive[cell_name] = pos_genes.split(",") if pos_genes != "nan" else []
        neg_genes = str(row["geneSymbolmore2"]).replace(" ", "")
        gs_negative[cell_name] = neg_genes.split(",") if neg_genes != "nan" else []
    logger.info(f"Loaded markers for {len(gs_positive)} cell types!")
    return {"gs_positive": gs_positive, "gs_negative": gs_negative}


def run_ScType(adata, cluster_col, tissue_type=None, custom_db_path=None):
    ScType_dir = "ScType_db"
    db_path = f"{ScType_dir}/ScTypeDB_full.xlsx"

    if custom_db_path:
        custom_db_path = f"{ScType_dir}/{custom_db_path}"

    logger.info(f"\n--- Processing {cluster_col} ---")

    # Pseudobulk logic
    unique_clusters = adata.obs[cluster_col].unique()
    cluster_means = {
        cluster: np.asarray(
            adata.X[adata.obs[cluster_col] == cluster].mean(axis=0)
        ).flatten()
        for cluster in unique_clusters
    }
    scRNAseqData_pb = pd.DataFrame(cluster_means, index=adata.var_names)

    if custom_db_path:
        logger.info("Using custom marker gene database.")
        gs_list = gene_sets_prepare(path_to_db_file=custom_db_path, tissue_type=None)
    else:
        logger.info(f"Using ScType's marker gene database (Tissue: {tissue_type}).")
        gs_list = gene_sets_prepare(path_to_db_file=db_path, tissue_type=tissue_type)

    # Run ScType
    es_max = sctype_score(
        scRNAseqData=scRNAseqData_pb,
        scaled=True,
        gs=gs_list["gs_positive"],
        gs2=gs_list["gs_negative"],
        check_genesymbols=False,
    )

    # Assign types
    cluster_to_type = {
        c: (es_max[c].idxmax() if es_max[c].max() >= 0 else "Unknown")
        for c in scRNAseqData_pb.columns
    }
    new_col_name = f"sctype_{cluster_col}"
    adata.obs[new_col_name] = (
        adata.obs[cluster_col].map(cluster_to_type).astype("category")
    )

    return adata, new_col_name
