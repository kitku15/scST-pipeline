import gc
import os
import logging
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from spatialdata_io import cosmx, xenium

logger = logging.getLogger(__name__)

# 1. Capture the true original function BEFORE any patching happens
_true_pandas_read_csv = pd.read_csv


def _safe_read_csv(filepath_or_buffer, *args, **kwargs):
    """Internal patched function to sanitize CosMx inputs."""

    # 2. Call the true original function to get a REAL DataFrame
    df = _true_pandas_read_csv(filepath_or_buffer, *args, **kwargs)

    conflict_cols = [c for c in ["cell_id"] if c in df.columns]
    if conflict_cols:
        logger.debug(f"Dropping redundant columns: {conflict_cols}")
        df = df.drop(columns=conflict_cols)

    if isinstance(filepath_or_buffer, (str, Path)) and "exprMat_file.csv" in str(
        filepath_or_buffer
    ):
        string_cols = df.select_dtypes(include=["object", "string"]).columns
        extra_drops = [c for c in string_cols if c not in ["cell_ID", "fov"]]
        if extra_drops:
            logger.debug(f"Dropping extra string cols from exprMat: {extra_drops}")
            df = df.drop(columns=extra_drops)

    return df


def convert_to_zarr(
    data_type: str, dataset_path: Path, dataset_id: str, zarr_path: Path
):
    if data_type not in ["CosMx", "Xenium"]:
        raise ValueError(f"Unsupported data type: {data_type}.")

    logger.info(f"Formatting {data_type} data from {dataset_path}")
    abs_dataset_path = os.path.abspath(dataset_path)

    try:
        if data_type == "CosMx":
            # Context manager: Patch only lives inside this 'with' block
            with patch("pandas.read_csv", side_effect=_safe_read_csv):
                sdata = cosmx(path=abs_dataset_path, dataset_id=dataset_id)
        elif data_type == "Xenium":
            sdata = xenium(abs_dataset_path)

        logger.info(f"Writing to Zarr: {zarr_path}")
        sdata.write(zarr_path, overwrite=True)

    except FileNotFoundError as err:
        logger.error(f"File not found during formatting: {err}")
        raise
    finally:
        if "sdata" in locals():
            del sdata
        gc.collect()
