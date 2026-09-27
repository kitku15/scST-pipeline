# Module 4: Spatial Visualization

!!! abstract "Overview"
    This module provides targeted, high-resolution visualizations of your spatial dataset. It maps specific genes of interest directly onto the tissue coordinates, highlights Field of View (FOV) boundaries, and generates strict Regions of Interest (ROI) grids for Xenium datasets.

## Parameters 

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chosen_cluster` | `String` | `"Final_Annotation"` | The primary cell type or cluster column in `adata.obs` to color the spatial tissue plot by. |
| `embedding_key` | `String` | `"X_umap"` | Specific embedding in `adata.obsm` to plot (e.g., `X_scVIVA`, `X_umap`, or `X_umap_n30_X_scANVI`). |
| `umap_color_columns` | `List[String]` | `["total_counts"]` | List of metadata columns (continuous metrics like `total_counts` or categorical labels) to display on the UMAP embedding. |
| `gene_list` | `List[String]` | `[]` | Specific genes to map onto physical spatial coordinates. Each gene generates an individual spatial expression subplot per sample. |
| `n_grid_x` | `Integer` | `10` | *(Xenium Only)* Number of horizontal grid divisions used to slice large tissue slides into smaller bounding box regions. |
| `n_grid_y` | `Integer` | `10` | *(Xenium Only)* Number of vertical grid divisions used to slice large tissue slides into smaller bounding box regions. |

## Example Config

```toml
[modules.ViewImages]
chosen_cluster = "CellTypist_majorityvoting_leiden_n10_r1.0"
embedding_key = "X_umap"
umap_color_columns = ["total_counts", "CellTypist_majorityvoting_leiden_n10_r1.0"]
gene_list = ["HLA-DRB", "KRT8", "TMSB4X", "TFF3", "IGHA1", "COX1"]
n_grid_x = 10
n_grid_y = 10
```

## Outputs

Outputs are saved in the `{analysis_name}_analysis/4_ViewImages/` directory:
<mark>Need to apply this folder structure: </mark>


```text
4_ViewImages/
├── adata.h5ad                                      # Updated AnnData object (with ROI grid bounding box assignments if Xenium)
├── gene_expression/
│   ├── gene_expression_{sample_id}.png            # Multi-panel spatial scatter plots showing gene expression across tissue
│   └── ...
├── embedding_plots/
│   ├── embedding_{embedding_key}_{column}.png      # UMAP or latent embeddings colored by metadata/annotations
│   └── ...
├── FOV_mapping_{sample_id}.png                     # (CosMx Only) Spatial plots displaying indexed FOV boundaries
└── ROI_grids/                                      # (Xenium Only) Artificial grid partitions for subsetting
    ├── Xenium_ROI_grid_{sample_id}.png             # Visual spatial plot of the generated X by Y grid
    └── Xenium_ROI_grid_{sample_id}.csv             # Coordinates table (min/max X and Y) for each grid box ID
```

### AnnData Updates {#hide-me}
Running this module processes visualizations and optionally updates spatial partition metadata in `adata.h5ad`:

* **`adata.obs`**:
    * `'grid_box_id'`: *(Xenium Only)* Categorical ID assigned to each cell indicating which artificial grid rectangle it falls within (used for downstream ROI-targeted analyses such as MuSpAn).
* **`adata.uns`**:
    * `'spatial_grid_dimensions'`: *(Xenium Only)* Stores the coordinate limits and resolution of the defined grid partitions.

## Core Libraries & Functions Used 
* **Note on Matplotlib:** The pipeline forces the `matplotlib.use("Agg")` backend. This ensures image generation does not crash on headless HPC cluster nodes that lack graphical displays.
* **[Squidpy](https://squidpy.readthedocs.io/en/stable/)**:
    * [`sq.pl.spatial_scatter()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_scatter.html): Generates spatial scatter plots on tissue coordinates for continuous transcript levels and categorical labels.
* **[Scanpy](https://scanpy.readthedocs.io/en/stable/)**:
    * [`sc.pl.embedding()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.embedding.html): Plots dimension reduction coordinates colored by specified cell-level metadata.
* **[Matplotlib](https://matplotlib.org/)**: 
    * `matplotlib.patches.Rectangle`: Draws FOV boundaries and artificial grid boxes directly onto physical coordinate coordinate plots.