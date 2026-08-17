"""Module for formatting CosMx/Xenium data into Zarr format."""

import gc
import os
import warnings
from logging import getLogger
from pathlib import Path

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
            # cwd = os.getcwd()
            # logger.info("Current Working Directory:", cwd)

            # exit()
            abs_dataset_path = os.path.abspath(dataset_path)
            sdata = cosmx(path=abs_dataset_path, dataset_id=dataset_id)
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
