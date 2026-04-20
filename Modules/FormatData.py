"""Module for formatting CosMx/Xenium data into Zarr format."""

from spatialdata_io import cosmx, xenium
from pathlib import Path 
from config import settings

def convert_to_zarr(data_type: str, dataset_path: Path, dataset_id: str, zarr_path: Path):
    """Convert Xenium/CosMx data to Zarr format."""

    if data_type not in ["CosMx", "Xenium"]:
        raise ValueError(f"Unsupported data type: {data_type}. Expected 'CosMx' or 'Xenium'.")
    
    if data_type == "CosMx":
        try:
            # Load CosMx data
            print("Reading CosMx data...")
            sdata = cosmx(path=dataset_path, dataset_id=dataset_id)
        except FileNotFoundError as err:
            print(f"File not found: {err}")
            raise err
    elif data_type == "Xenium":
        try:
            # Load Xenium data
            print("Reading Xenium data...")
            sdata = xenium(dataset_path)
        except FileNotFoundError as err:
            print(f"File not found: {err}")
            raise err

    try:
        # Write to Zarr format
        print("Writing to Zarr...")
        sdata.write(zarr_path, overwrite=True)
    except ValueError as err:
        print(f"Failed writing to Zarr: {err}")
        raise err


if __name__ == "__main__":

    dataset_path = settings['io']['dataset_dir']
    dataset_id = settings['io']['dataset_id']
    zarr_path = settings['io']['zarr_dir']
    data_type = settings['project']['data_type']

    convert_to_zarr(data_type, dataset_path, dataset_id, zarr_path)