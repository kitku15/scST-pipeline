# Module 5: Spatial Statistics (Squidpy)

!!! abstract "Overview"
    This module analyzes the physical architecture of the tissue. It builds Delaunay neighborhood graphs to determine which cells are physically touching, and then calculates statistical metrics like **Neighborhood Enrichment**, **Co-occurrence Probabilities**, **Centrality Scores**, and **Moran's I** spatial autocorrelation. 
    
    If multiple conditions (e.g., Healthy vs Disease) are provided, the module aggregates these metrics across samples to compute consensus and differential spatial interaction signatures.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chosen_cluster` | `String` | `"CellTypist_majorityvoting_leiden_n10_r1.0"` | The cell type or cluster column in `adata.obs` used to compute Neighborhood Enrichment, Centrality, and Co-occurrence. |
| `condition_key` | `String` | `""` | *(Optional)* Metadata column in `adata.obs` defining experimental cohorts (e.g., `"TreatmentResponse"`). Aggregates stats across samples sharing a condition. |
| `reference_condition` | `String` | `""` | *(Optional)* The baseline condition (e.g., `"Healthy"`). Used to generate difference ($\Delta$) heatmaps against other conditions. |
| `sample_key` | `String` | *(Required)* | Metadata column in `adata.obs` denoting individual samples/slides to iterate over. |
| `skip_compute` | `Boolean` | `False` | If `True`, activates "Rescue Mode" to bypass heavy computation and reload previously generated JSONs/CSVs for aggregation. |

*Note: To optimize performance, the pipeline automatically subsamples cells by 50% (`fraction=0.5`) for computationally heavy metrics like Co-occurrence and Moran's I.*

## Example Config

```toml
[modules.SpatialStat]
chosen_cluster = "CellTypist_majorityvoting_leiden_n10_r1.0"
condition_key = "TreatmentResponse"
reference_condition = "Healthy"
n_jobs = -1
sample_key = "SampleID"
skip_compute = false

```

## Outputs

Outputs are organized per sample and, if cohorts are defined, aggregated across experimental conditions:

```text
5_SpatialStat/
├── adata.h5ad                                           # Updated AnnData object containing spatial graphs and metrics
├── {sample_id}/                                         # Per-sample spatial statistics
│   ├── nhood_enrichment_{sample_id}.png / .json         # Heatmap and raw matrix of enrichment z-scores
│   ├── co_occurrence_{sample_id}.png / .json            # Conditional co-occurrence probability curves and raw data
│   ├── centrality_scores_{sample_id}.png / .json        # Closeness, degree, clustering plots and raw data
│   └── moranI_results_{sample_id}.csv                   # Table of Moran's I autocorrelation statistics per gene
└── Aggregated_Results/                                  # (Generated if condition_key is provided)
    ├── Consensus_Nhood_{condition}.png                  # Mean neighborhood enrichment across replicates in a condition
    ├── Diff_Nhood_{condition}_vs_{reference}.png        # Delta (Δ) heatmap showing gained or lost interactions
    ├── Boxplot_Degree_Centrality.png                    # Degree centrality distribution across cohorts
    ├── Aggregated_Centrality.csv                        # Merged centrality metrics across all samples
    ├── Boxplot_Morans_I.png                             # Distribution of spatial autocorrelation for top variable genes
    └── Aggregated_MoransI.csv                           # Merged Moran's I statistics across all samples
```

## Core Libraries & Functions Used
* **[Squidpy](https://squidpy.readthedocs.io/en/stable/)**:
    * [`sq.gr.spatial_neighbors()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html): Builds Delaunay triangulation graphs based on cell centroid spatial coordinates (`spatial_key='spatial'`).
    * [`sq.gr.nhood_enrichment()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.nhood_enrichment.html): Tests whether pairs of cell types are spatial neighbors more (or less) frequently than expected under random label permutations.
    * [`sq.gr.co_occurrence()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.co_occurrence.html): Computes conditional co-occurrence probabilities as a function of increasing Euclidean radius.
    * [`sq.gr.centrality_scores()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.centrality_scores.html): Quantifies cell type network topology metrics (degree, closeness centrality, clustering coefficient).
    * [`sq.gr.spatial_autocorr()`](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_autocorr.html): Calculates Moran's I global spatial autocorrelation to identify non-randomly distributed transcripts across tissue space.