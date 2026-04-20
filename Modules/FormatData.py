"""Module for formatting CosMx data into Zarr format."""

from spatialdata_io import cosmx
from pathlib import Path 


def convert_CosMx_to_zarr(cosmx_path: Path, dataset_id: str, zarr_path: Path):
    """Convert CosMx data to Zarr format."""
    try:
        # Load CosMx data
        print("Reading CosMx data...")
        sdata = cosmx(path=cosmx_path, dataset_id=dataset_id)
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

    # path=r"C:\Users\bunga\python\Project2\CosMx\Lung13+SMI+Flat+data\Lung13\Lung13-Flat_files_and_images"
    path=r"C:\Users\bunga\python\Project2\CosMx\Kitam"
    dataset_id="Quarter"
    zarr_path = Path("Kitam.zarr")
    
    convert_CosMx_to_zarr(path, dataset_id, zarr_path)