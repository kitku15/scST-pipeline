# Module 1 & 1b: Quality Control & Merging

!!! abstract "Overview"
    **Module 1** calculates quality control metrics, filters out low-quality 'cells' based on user-defined thresholds and log normalizes transcript counts. It can optionally append external metadata (e.g. condition, annotations) or updated cell geometries (currently works with Proseg). 
    
    **Module 1b** optionally merges multiple QC'd slides into a single `adata.h5ad` dataset, making cell and sample IDs globally unique.

## Parameters

### Filtering Thresholds {#hide-me}
| Parameter | Description | Why adjust this? |
|-----------|-------------|------------------|
| `min_counts` | Minimum total transcripts required to keep a cell. | Filters out debris, or poorly segmented areas. |
| `min_genes` | Minimum unique genes required to keep a cell. | Removes artifacts expressing too few genes. |
| `min_area` / `max_area` | Allowed cell area (µm²). | Drops segmentation errors (e.g., massive fused cells or tiny fragments). |
| `min_dapi` | *(CosMx only)* Minimum Mean DAPI signal. | Ensures the segmented shape actually contains a nucleus. |
| `min_cells` | Min. cells a gene must be expressed in. | Drops rare genes to reduce matrix dimensions. |

 <mark>Go over DAPI filtering best practices for Xenium and CosMx </mark>


### Integration Parameters {#hide-me}
| Parameter | Description |
|-----------|-------------|
| `fov_metadata_path` | *(Optional)* Path to a CSV file to map external metadata (e.g. condition, annotations) to cells via a common column. |
| `metadata_join_col` | *(Optional unless appenidng metadata )*The column name that exists in both your CSV and your data, which the pipeline will use to link them together.|
| `proseg_zarr_path`  | *(Optional)* Path to an alternate SpatialData Zarr containing improved segmentation boundaries to replace the default ones. |

## Example Config

```toml
[modules.QualityControl]
min_counts = 20
min_cells = 3
min_genes = 10
min_area = 10.0
max_area = 200.0
min_dapi = 100.0
fov_metadata_path = "data/metadata_{slide_name}.csv" 
```
<mark>I should change fov_metadata_path to something else </mark>

## Outputs

Outputs are organized per slide under `1_QualityControl/` and merged centrally under `1b_MergeData/`:

```text
{analysis_name}_analysis/
├── 1_QualityControl/
│   ├── {slide_1_name}/
│   │   ├── adata.h5ad                                # Cleaned, filtered, and log-normalized AnnData for Slide 1
│   │   ├── Area_scatter.png                          # Spatial tissue scatter plots colored by Area
│   │   ├── Mean.DAPI_scatter.png                     # Spatial tissue scatter plots colored by mean DAPI
│   │   ├── cell_summary_histograms.png               # QC distributions (counts, genes, area, DAPI) with cutoff lines
│   │   ├── qc_metrics.csv                            # Per-cell calculated QC statistics table
│   │   └── qc_thresholds.json                        # Applied filtering threshold parameters (for web visualizer)
│   ├── {slide_2_name}/
│   │   ├── ...
│   └── ...

└── 1b_MergeData/
    └── adata.h5ad                                    # AnnData object containing all merged slides
```

---

### AnnData Updates {#hide-me}

#### Module 1  {#hide-me}

* `1_QualityControl/{slide_name}/adata.h5ad`
* **`adata.X`**: Library-size normalized and $\log_2$-transformed expression matrix (Base 2 log).
* **`adata.layers['counts']`**: Preserved raw, unnormalized integer transcript counts (strictly required for downstream scVI and PyDESeq2).
* **`adata.obs`**:
    * `'n_counts'` / `'total_counts'`: Total transcript molecules detected per cell.
    * `'n_genes'` / `'n_genes_by_counts'`: Total distinct genes detected per cell.
    * `'cell_area'`: Segmented cell surface area ($\mu\text{m}^2$).
    * `'mean_dapi'` / `'dapi_intensity'`: Mean nuclear stain intensity *(CosMx only)*.
    * *(External Metadata)*: Any sample-level or FOV-level clinical attributes joined from `fov_metadata_path`.
* **`adata.var`**:
    * `'n_cells'` / `'n_cells_by_counts'`: Number of cells expressing each gene.
    * `'control'`: Boolean flag indicating if the gene is a negative control probe. The pipeline automatically detects these using platform-specific regex (`^NegPrb|^SystemControl` for CosMx, `control_probe|control_codeword` for Xenium).
* **`adata.obs`**: *(Additional)*
    * `'control_counts'`: Total transcripts matching the negative control probes.
* **`adata.uns`**:
    * `'qc_thresholds'`: Dictionary of the filtering bounds applied.

#### Module 1b  {#hide-me}

* `1b_MergeData/adata.h5ad`
* **`adata.obs_names`**: Barcodes are prefixed to ensure global uniqueness across slides (e.g., `Slide_1_cell_001`, `Slide_2_cell_001`).
* **`adata.obs`**:
    * `'{batch_key}'`: Slide/batch identifier (e.g., `"Slide_1"`, `"Slide_2"`).
    * `'{sample_key}'`: Tissue region or sample identifier (e.g., `"Slide_1_Healthy_1"`).
* **`adata.var`**: Harmonized gene index across all merged slides. **Note:** Merging uses a strict `join="inner"` strategy. Only genes that survived QC across *every single slide* are retained in the final dataset.

##  Core Libraries & Functions Used
* **[Scanpy](https://scanpy.readthedocs.io/en/stable/)**:
    * [`sc.pp.calculate_qc_metrics()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html): Calculates total counts and unique genes per cell.
    * [`sc.pp.filter_genes()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.filter_genes.html): Removes genes expressed in too few cells.
    * [`sc.pp.normalize_total()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html): Normalizes cells to a target sum (10,000).
    * [`sc.pp.log1p()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html): Log-transforms the normalized counts (Base 2).
* **[Squidpy](https://squidpy.readthedocs.io/en/stable/)**:
    * [`sq.pl.spatial_scatter()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_scatter.html): Plots QC metrics mapped directly onto the tissue coordinate space.