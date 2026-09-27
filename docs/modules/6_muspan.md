# Module 6: MuSpAn (Morphometrics & Spatial Networks)

!!! abstract "Overview"
    **MuSpAn** (Multi-Scale Spatial Analysis) performs continuous-space spatial modeling and single-cell morphometric profiling. Unlike standard point-based coordinate methods, MuSpAn uses segmentation boundaries to treat cells as geometric polygons. 
    
    This module enables targeted extraction of specific **Regions of Interest (ROIs)** or Field of Views (FOVs) to evaluate shape morphometrics (area, circularity, principal axis), build boundary-to-boundary proximity networks, and compute continuous spatial statistics including the Cross-Pair Correlation Function (Cross-PCF).

## Parameters\

### General & ROI Selection Settings {#hide-me}

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `selection_names` | `List[String]` | `["Whole_Sample"]` | User-defined identifiers for the analyzed tissue regions (e.g., `["myeloid_niche"]`). |
| `chosen_cluster` | `String` | `"CellTypist_majorityvoting_leiden_n10_r1.0"` | Cell type annotation column in `adata.obs` used to group cells for proximity and correlation analysis. |
| `condition_key` | `String` | `""` | *(Optional)* Metadata column in `adata.obs` defining biological cohorts (e.g., `"TreatmentResponse"`). Aggregates contact composition and morphometrics across conditions. |
| `reference_condition` | `String` | `""` | *(Optional)* Baseline cohort identifier (e.g., `"Healthy"`) used for differential comparison boxplots. |
| `transcripts` | `List[String]` | `[]` | Specific transcripts to project as point coordinates overlaid directly onto cell boundary polygons. |
| `[modules.MuSpan.selected_fovs]` | `Dict` | `{}` | Mapping table defining which FOVs (CosMx) or artificial Grid Box IDs (Xenium) to subset per slide/sample. |

* **Custom ROI Bypass (Module 6a):** If you wish to manually define a custom cell selection rather than using FOVs/Grid Boxes, you can place a CSV named `{selection_name}.csv` inside the target output folder (`6_MuSpan/{sample_id}/{selection_name}/`). The pipeline will automatically detect it, bypass Module 6a, and construct the MuSpAn domain using your specific cells.

### Network & Proximity Settings {#hide-me}

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cell_types` | `List[String]` | `[]` | Specifies two cell types from `chosen_cluster` for targeted pairwise Cross-PCF computation. |
| `network_celltypes` | `List[String]` | `[]` | Cell types to highlight individually in stacked contact composition charts; unlisted types are grouped into "Other". |
| `max_distance` | `Integer` | `200` | Maximum Euclidean physical distance threshold (µm) for defining cell-cell proximity networks. |
| `distance_list` | `List[Int]` | `[50, 90, 120]` | Metric distances (µm) evaluated for multi-scale proximity network construction. |
| `k_list` | `List[Int]` | `[2, 5, 10, 15]` | List of $k$ nearest neighbors evaluated during KNN graph generation. |
| `min_edge_distance` | `Integer` | `0` | Minimum edge length threshold (µm) for centroid-based networks. |
| `max_edge_distance` | `Integer` | `200` | Maximum edge length cutoff (µm) for centroid-based networks. |
| `min_edge_distance_shape` | `Integer` | `0` | Minimum edge distance threshold between polygon boundaries. |
| `max_edge_distance_shape` | `Integer` | `1` | Maximum edge distance cutoff between polygon boundaries (default `1` captures direct physical contact). |
| `cellboundary_label` | `String` | `"Cell boundaries"` | Key in the raw spatial store containing the cell segmentation polygon collections. |

* **Note on CosMx Physical Scaling:** During extraction, the pipeline explicitly multiplies all CosMx pixel coordinates by a hardcoded physical scaling factor (`COSMX_PIXEL_SIZE_UM = 0.12`) to ensure MuSpAn calculates Area and Perimeter accurately in micrometers (µm).

## Example Config

```toml
[modules.MuSpan]
selection_names = ["myeloid_niche"]
chosen_cluster = "CellTypist_majorityvoting_leiden_n10_r1.0"
condition_key = "TreatmentResponse"
reference_condition = "Healthy"
transcripts = ["HLA-DRB", "KRT8", "TMSB4X", "TFF3", "IGHA1", "COX1"]
cell_types = ["Macrophage", "CD8 T cell"]
network_celltypes = ["Macrophage", "CD8 T cell", "Epithelial cell"]

# Network & Graph thresholds
max_distance = 200
distance_list = [50, 90, 120]
k_list = [2, 5, 10, 15]
min_edge_distance = 0
max_edge_distance = 200
min_edge_distance_shape = 0
max_edge_distance_shape = 1
cellboundary_label = "Cell boundaries"

[modules.MuSpan.selected_fovs]
Slide_1 = ["1", "2", "3", "4"]
```

## Outputs

Outputs are organized per ROI under `{analysis_name}_analysis/6_MuSpan/{sample_id}/{selection_name}/`:
<mark>Need to apply this folder structure: </mark>

```text
6_MuSpan/
├── Aggregated_Results/                                   # (Generated only if condition_key is provided)
│   ├── Aggregated_Area.png                               # Boxplots comparing area between conditions
│   ├── Aggregated_Circularity.png                        # Boxplots comparing circularity between conditions
│   ├── Consensus_Contacts_{cond}.png                     # Heatmap of mean cell-cell contacts per condition
│   └── Diff_Contacts_{cond}_vs_{reference_condition}.png # Delta contact heatmap against reference condition
└── {sample_id}/
    └── {selection_name}/
        ├── domain/
        │   ├── muspan_domain_visualization.png           # Polygon boundaries of extracted cells and transcripts
        │   ├── muspan_cell_centroids.png                 # (Xenium only) Centroid representation of cells
        │   └── muspan_cell_centroids_n_boundaries.png    # (Xenium only) Boundaries overlaid with centroids
        ├── networks/
        │   ├── muspan_knn.png                            # KNN connectivity graphs across values of k
        │   ├── muspan_delaunay.png                       # Continuous Delaunay triangulation mesh
        │   ├── muspan_proximity_point.png                # Centroid distance-threshold networks
        │   └── muspan_proximity_shape.png                # True polygon boundary-to-boundary contact networks
        ├── spatial_stats/
        │   ├── cross_pair_correlation_function_all.png                       # Pairwise Cross-PCF curves for all cell types
        │   ├── cross_pair_correlation_function_{cell_type_1}_{cell_type_2}.png # Targeted Cross-PCF plot for selected cell_types
        │   ├── visualize_{cell_type_1}_{cell_type_2}.png                     # Spatial plot highlighting targeted cell types
        │   └── cross_pcf_all.json                                            # Raw Cross-PCF values and distance intervals
        ├── shape_analysis/
        │   ├── ms_shapes.png                             # Spatial maps of Area, Perimeter, Convexity, and Circularity
        │   ├── ms_praxis.png                             # Spatial map showing principal axis orientation angles
        │   └── morphometrics.csv                         # Quantified single-cell geometric measurements table
        └── proximity_analysis/
            ├── contact_statistics.png                    # Stacked bar charts showing cell-cell boundary contact compositions
            ├── zoomed_khop_neighborhood.png              # High-resolution boundary contact plot around target populations
            ├── global_proximity_network.png              # Full ROI map showing connectivity edges between cells
            └── contact_composition.csv                   # Raw matrix of neighboring contact frequencies per cell type
```

## Core Libraries & Functions Used
* **[MuSpAn](https://github.com/ImperialCollegeLondon/MuSpAn)**:
    * `ms.geometry`: Parses boundary polygon coordinates and quantifies geometric properties (Area, Perimeter, Convexity, Principal Axis).
    * `ms.networks.generate_network()`: Builds network architectures, including $k$-nearest neighbors, Delaunay triangulation, centroid proximity, and polygon boundary intersection graphs.
    * `ms.spatial_statistics.cross_pair_correlation_function()`: Quantifies spatial point-pattern associations (attraction vs. repulsion) between distinct cell types as a continuous function of distance.