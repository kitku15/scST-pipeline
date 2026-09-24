import logging
from pathlib import Path
import anndata as ad

logger = logging.getLogger(__name__)


def run_merge(input_files: list[str], slide_names: list[str], module_dir: Path):
    """Merges multiple Spatial AnnDatas in memory."""
    module_dir = Path(module_dir)
    module_dir.mkdir(parents=True, exist_ok=True)

    adatas_to_merge = {}

    for file_path, slide_id in zip(input_files, slide_names):
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Cannot find input file: {path}")

        logger.info(f"Loading dataset for merge: {slide_id}")

        # We load into memory (removed backed="r" to avoid anndata slicing bugs)
        adata = ad.read_h5ad(path)

        # Safety check: If QC filtered out ALL cells, skip this slide
        if adata.n_obs == 0:
            logger.warning(f"Slide {slide_id} has 0 cells after QC! Skipping.")
            continue

        logger.info(
            f"Slide {slide_id} loaded: {adata.n_obs} cells, {adata.n_vars} genes."
        )

        # Ensure sample IDs are globally unique across slides
        if "sample_id" in adata.obs.columns:
            adata.obs["sample_id"] = f"{slide_id}_" + adata.obs["sample_id"].astype(str)

        adata.var_names = adata.var_names.astype(str)

        adatas_to_merge[slide_id] = adata

    if not adatas_to_merge:
        raise ValueError(
            "No valid datasets remaining to merge! All datasets had 0 cells."
        )

    logger.info("Concatenating datasets...")

    # anndata automatically uses the dictionary keys (slide_id) as the categories.
    adata_merged = ad.concat(
        adatas_to_merge, join="inner", merge="same", label="slide_id", index_unique="_"
    )

    out_path = module_dir / "adata.h5ad"
    logger.info(f"Writing merged matrix to disk: {out_path}")
    adata_merged.write_h5ad(out_path)

    logger.info(
        f"Merged dataset successfully saved: {adata_merged.n_obs} total cells, {adata_merged.n_vars} genes"
    )
    return out_path
