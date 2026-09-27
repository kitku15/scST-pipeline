# Supported Inputs & Data Formats

!!! abstract "Overview"
    scSpatial-Kit natively supports flat-file outputs from the major sub-cellular spatial transcriptomics platforms. Module 0 (`FormatData`) automatically handles parsing these vendor-specific formats into a unified format for downstream analysis.

### NanoString CosMx SMI

The pipeline natively parses CosMx output directories. Point it to your root dataset directory containing:

* `*exprMat_file.csv` (Expression matrix)
* `*metadata_file.csv` (Cell metadata and centroids)
* `*polygons.csv` (Cell boundary coordinates)

### 10x Genomics Xenium

The pipeline natively supports standard Xenium output bundles. Point it to the directory containing:

* `cell_feature_matrix.h5`
* `cells.csv.gz`
* `cell_boundaries.parquet` (or `cell_boundaries.csv.gz`)

---

## What if my dataset is from another platform? 

wait for an update I guess

<!-- If you have data from Vizgen MERSCOPE, 10x Visium, Resolve Biosciences, or a public dataset that is *already processed*, you can still use scSpatial-Kit by skipping Module 0 and using the **Bring Your Own Data (BYOD)** feature.

### The BYOD Approach (Standard AnnData)
The core of scSpatial-Kit runs on standard Python `AnnData` objects (`.h5ad`). To use an unsupported platform:

1. **Format your data manually:** Create an `AnnData` object using `Scanpy`.
2. **Ensure spatial coordinates exist:** Your object *must* have a `.obsm` key containing physical coordinates (e.g., `adata.obsm['spatial']`).
3. **Skip Module 0:** In your configuration file, remove `0_FormatData` from the `modules` list.
4. **Point to your custom file:** Use the `custom_input_adata` parameter in your `[io]` config block and set your `entry_point`.

```toml
[io]
entry_point = "1" # Start directly at Quality Control
custom_input_adata = "path/to/my_custom_formatted_data.h5ad"
```

!!! tip "Migrating from Seurat (R)?"
    If your data is currently a Seurat object in R, you can use the [`SeuratDisk`](https://mojaveazure.github.io/seurat-disk/) library in R to convert your `.rds` file to an `.h5ad` file, and then inject it into scSpatial-Kit! -->
