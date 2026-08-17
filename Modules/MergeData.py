"""Module to merge multiple Spatial Transcriptomics Slides (datasets)."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_merge(input_files, slide_names, module_dir):
    """
    Merges multiple AnnData objects into a single AnnData object.
    Prefixes cell IDs with slide names to ensure global uniqueness.
    """
    module_dir = Path(module_dir)
    module_dir.mkdir(parents=True, exist_ok=True)

    adatas = []

    # Prep each dataset
    for i, file_path in enumerate(input_files):
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Cannot find input file: {file_path}")

        logger.info(f"Loading dataset {i + 1}/{len(input_files)}: {file_path}")
        adata = sc.read_h5ad(file_path)
        slide_id = slide_names[i]

        # Add metadata column
        adata.obs["slide_id"] = slide_id

        # Make Cell IDs globally unique
        adata.obs_names = [f"{slide_id}_{cell_id}" for cell_id in adata.obs_names]

        # Make Sample IDs globally unique
        if "sample_id" in adata.obs.columns:
            adata.obs["sample_id"] = slide_id + "_" + adata.obs["sample_id"].astype(str)

        # Ensure gene names are strings
        adata.var_names = adata.var_names.astype(str)

        adatas.append(adata)

    logger.info("Concatenating datasets...")
    # join='inner' ensures only genes present in all slides are kept
    adata_merged = ad.concat(adatas, join="inner", merge="same")

    # Verify and log the merge
    logger.info(f"Total cells: {adata_merged.n_obs}")
    logger.info(f"Total genes: {adata_merged.n_vars}")
    logger.info("Cells per slide:")
    for slide, count in adata_merged.obs["slide_id"].value_counts().items():
        logger.info(f"  {slide}: {count}")

    # Save the merged dataset
    out_path = module_dir / "adata.h5ad"
    adata_merged.write_h5ad(out_path)
    logger.info(f"Merged dataset saved to {out_path}")

    return out_path
