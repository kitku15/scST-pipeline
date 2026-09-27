# Module 10: Web Visualizer Export

!!! abstract "Overview"
    This module aggregates and compresses outputs from all upstream modules (Modules 1 through 9) into a `.tar` file for interactive exploration in the **[Spatial-VisKit](https://github.com/kitku15/Spatial-VisKit)** web application. 
    

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `primary_annotation` | `String` | `"Broad_Celltype"` | The default cell-type column in `adata.obs` displayed upon loading the visualizer. |
| `microenv_col` | `String` | `"spatial_microenvironment"` | Column in `adata.obs` defining spatial microenvironments, used for Vitessce hierarchical grouping and filtering. |
| `annotation_columns` | `List[String]` | `["leiden", "CellTypist", "sctype", "cluster"]` | Substring match list. Any column in `adata.obs` containing one of these patterns is exported to the web application dropdown menus. |

## Example Config

```toml
[modules.WebVisPrep]
primary_annotation = "CellTypist_majorityvoting_leiden_n10_r1.0"
microenv_col = "spatial_microenvironment"
annotation_columns = ["leiden", "CellTypist", "ScType", "Final_Annotation", "Broad_Celltype"]
```

## Outputs

Outputs are prepared in the `{analysis_name}_analysis/10_WebVis/` directory.

```text
10_WebVis/
├── aux_data/                                           # JSON/CSV summaries for frontend visualisations
│   ├── causal_ccc/                                     # Condition-specific LIANA+ causal network data
│   ├── conditions_de_analysis/                         # Sub-setted differential expression (Condition vs Condition)
│   ├── de_analysis/                                    # Standard cluster marker DE results
│   ├── qc/                                             # QC metric histograms and thresholds
│   ├── segmentations/                                  # Cell boundary polygons scaled and aligned per sample/slide
│   ├── spatial_stats/                                  # Morans I, centrality, and MuSpAn metric exports
│   ├── cpdb_edges.json                                 # CellPhoneDB significant interactions
│   └── tf_heatmap_data.json                            # Decoupler Transcription Factor activity heatmap data
├── adata_{analysis_name}_web.zarr/                     # Chunked, web-optimized primary AnnData Zarr store
├── adata_{analysis_name}_tf_web.zarr/                  # Chunked AnnData Zarr containing TF activities (if available)
├── dataset_config.json                                 # Initialisation schema for backend
└── {module_name}.tar                                   # Compressed archive of the entire folder ready for upload
```

## Core Libraries & Functions Used
* **[SpatialData](https://spatialdata.scverse.org/)**: Reads boundary polygons (`.parquet` or `.csv`) from the slide-level Zarr stores, scales coordinate axes, and exports downsampled polygon contours.
* **[Zarr-Python](https://zarr.readthedocs.io/)** & **[Anndata](https://anndata.readthedocs.io/)**:
    * `adata.write_zarr()`: Serializes processed cell metadata, gene expression layers, and dimensionality reduction matrices into structured chunked storage.
* **Standard Python (`json`, `tarfile`)**: Encodes communication graphs and statistical distributions into JSON schemas and builds the transfer-ready `.tar` file.