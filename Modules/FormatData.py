"""Module for formatting CosMx/Xenium data into Zarr format."""

import gc
import os
import warnings
from logging import getLogger
from pathlib import Path
import pandas as pd

from spatialdata_io import cosmx, xenium

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def convert_to_zarr(
    data_type: str, dataset_path: Path, dataset_id: str, zarr_path: Path
):
    """Convert Xenium/CosMx data to Zarr format."""

    if data_type not in ["CosMx", "Xenium"]:
        raise ValueError(
            f"Unsupported data type: {data_type}. Expected 'CosMx' or 'Xenium'."
        )

    if data_type == "CosMx":
        try:
            # Load CosMx data
            logger.info("Reading CosMx data...")
            logger.info(f"Dataset path: {dataset_path}, Dataset ID: {dataset_id}")
            abs_dataset_path = os.path.abspath(dataset_path)

            # Intercept and clean the CSV for spatialdata_io (newer CosMx has a string column in the expression matrix)
            orig_read_csv = pd.read_csv
            
            def patched_read_csv(filepath_or_buffer, *args, **kwargs):

                df = orig_read_csv(filepath_or_buffer, *args, **kwargs)
                
                # address spatialdata validation error
                # Drop 'cell_id' from ANY loaded file to avoid case-variant conflicts with 'cell_ID'
                conflict_cols = [c for c in ['cell_id'] if c in df.columns]
                if conflict_cols:
                    logger.info(f"⚠️ FIX: Dropping redundant columns to avoid case-variant errors: {conflict_cols}")
                    df = df.drop(columns=conflict_cols)
                    
                # address scipy sparse error: If it's the expression matrix, drop ANY remaining string columns
                if isinstance(filepath_or_buffer, (str, Path)) and "exprMat_file.csv" in str(filepath_or_buffer):
                    string_cols = df.select_dtypes(include=['object', 'string']).columns
                    
                    # Protect 'cell_ID' and 'fov' just in case they were somehow read as strings
                    extra_drops = [c for c in string_cols if c not in ['cell_ID', 'fov']]
                    
                    if extra_drops:
                        logger.info(f"⚠️ FIX: Dropping extra string columns from exprMat: {extra_drops}")
                        df = df.drop(columns=extra_drops)
                        
                return df
            
            # Temporarily replace pandas read_csv with patched version
            pd.read_csv = patched_read_csv
            
            try:
                sdata = cosmx(path=abs_dataset_path, dataset_id=dataset_id)
            finally:
                # restore the original pandas function so we don't break the rest of the pipeline
                pd.read_csv = orig_read_csv

        except FileNotFoundError as err:
            logger.info(f"File not found: {err}")
            raise err
    elif data_type == "Xenium":
        try:
            # Load Xenium data
            logger.info("Reading Xenium data...")
            sdata = xenium(dataset_path)
        except FileNotFoundError as err:
            logger.info(f"File not found: {err}")
            raise err

    try:
        # Write to Zarr format
        logger.info("Writing to Zarr...")
        sdata.write(zarr_path, overwrite=True)
        del sdata  # Free up memory after writing to Zarr
        gc.collect()
    except ValueError as err:
        logger.info(f"Failed writing to Zarr: {err}")
        raise err
