# Module 0: Data Formatting

Parses raw datasets into SpatialData format.

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `base_raw_dir` | `String` | The root directory containing your raw vendor output folders. |
| `base_zarr_dir` | `String` | The directory where `.zarr` files will be saved. |
| `dataset_id` | `String` | *(Optional, CosMx only)* Prefix of the CosMx flat files. |

## Example Config 
These settings are defined in the `[io]` block of your `config.toml` file.

```toml
[io]
base_raw_dir = "data/raw_data"
base_zarr_dir = "data_zarrs"
```

## Outputs
* **`{slide_name}.zarr/`**: A SpatialData Zarr directory containing tables, shapes (boundaries), and points (transcripts).

##  Core Libraries & Functions Used
* **[SpatialData-IO](https://spatialdata.scverse.org/projects/io/en/latest/)**: 
    * `spatialdata_io.cosmx()`: Parses CosMx flat files.
    * `spatialdata_io.xenium()`: Parses Xenium `.parquet` and `.csv.gz` files.
* **[SpatialData](https://spatialdata.scverse.org/en/latest/)**: 
    * `SpatialData.write()`: Saves the object as a `.zarr` store.