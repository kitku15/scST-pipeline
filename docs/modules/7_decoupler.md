# Module 7: TF Enrichment & Pseudobulking (Decoupler)

!!! abstract "Overview"
    This module infers **Transcription Factor (TF) activity** from spatial transcriptomics data. Because spatial targeted gene panels often capture low or variable counts for transcription factors directly, [Decoupler](https://decoupler-py.readthedocs.io/en/latest/) leverages Gene Regulatory Networks (GRNs) to estimate TF activity from the enriched expression of downstream target genes.
    
    Additionally, this module aggregates single-cell transcriptomes into **pseudobulk** profiles stratified by cell type and biological sample, which serve as inputs for downstream differential expression (Module 9) and causal communication modeling (Module 8c).

## Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `organism` | `String` | `"human"` | Organism model (`"human"`, `"mouse"`, or `"rat"`). Handles gene nomenclature and capitalization mapping across GRNs automatically. |
| `grn` | `String` | `"collectri"` | Gene Regulatory Network database to query (`"collectri"` or `"dorothea"`). |
| `celltype_key` | `String` | `"CellTypist_majorityvoting_leiden_n10_r1.0"` | Column in `adata.obs` defining cell identities or microenvironments used for cluster-level TF contrast and pseudobulking. |
| `active_tfs_file_name` | `String` | `"active_tf.txt"` | Filename for the exported text file storing statistically active TFs (used downstream by CellPhoneDB v5 CellSign). |
| `dorothea_levels` | `List[String]` | `["A", "B"]` | Confidence tiers included when using DoRothEA (`"A"` through `"D"`, with `"A"` indicating the highest empirical support). Ignored if `grn = "collectri"`. |

## Example Config

```toml
[modules.Decoupler]
organism = "human"
grn = "collectri"
celltype_key = "CellTypist_majorityvoting_leiden_n10_r1.0"
active_tfs_file_name = "active_tf.txt"
dorothea_levels = ["A", "B", "C"]
```

## Outputs

Outputs are saved under the `{analysis_name}_analysis/7_Decoupler/` directory:

```text
7_Decoupler/
├── tf_activity_scores.h5ad                   # Transformed AnnData: cells × TF activity estimates (ULM scores)
├── pseudobulk.h5ad                           # Pseudobulk-aggregated AnnData object (samples × cell types)
├── active_tf.txt                             # Text file list of active TFs for CellPhoneDB CellSign
├── top3tfs_heatmap.png                       # Clustered heatmap of top 3 TF activities per cell group
├── tf_heatmap_data.json                      # JSON export of clustered heatmap data for the web interface
└── microenvironment_tfs_{sample}.png         # Spatial scatter plots mapping the top TFs onto tissue coordinates
```

### AnnData Objects Generated {#hide-me}
Instead of appending data to the main AnnData object, this module extracts the results into two specialized, standalone `.h5ad` files for downstream use:

* **`tf_activity_scores.h5ad`**: 
    * A transformed AnnData object where the main `.X` matrix contains the continuous **TF activity estimates** (ULM scores) rather than gene expression. 
    * Dimensions are `Cells × Transcription Factors`. 
    * Retains the original single-cell `.obs` metadata (cell types, sample IDs, etc.) for spatial plotting and clustering.
* **`pseudobulk.h5ad`**: 
    * An aggregated AnnData object where raw gene counts are summed across cells of the same type within the same biological sample. 
    * Dimensions are `(Samples × Cell Types) × Genes`. 
    * The `.obs` dataframe includes custom filtering metrics calculated during generation: `'manual_total_counts'` (minimum 1000) and `'manual_n_cells'` (minimum 10).

## Core Libraries & Functions Used
* **[Decoupler](https://decoupler-py.readthedocs.io/en/latest/)**:
    * `dc.op.collectri()` / `dc.op.dorothea()`: Retrieves curated regulon gene sets linking transcription factors to mode-of-regulation target genes.
    * `dc.mt.ulm()`: Fits Univariate Linear Models to evaluate consensus regulator activity scores per single cell.
    * `dc.pp.pseudobulk()`: Sums raw integer transcript counts across cell types and biological samples to construct replicated pseudobulk matrices.
* **[Scanpy](https://scanpy.readthedocs.io/en/stable/)**:
    * [`sc.pl.matrixplot()`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.matrixplot.html): Generates hierarchically clustered matrix plots of regulatory scores across cell types.
* **[Squidpy](https://squidpy.readthedocs.io/en/stable/)**:
    * `sq.pl.spatial_scatter()`: Maps the top transcription factor activities directly onto spatial tissue coordinates per sample.