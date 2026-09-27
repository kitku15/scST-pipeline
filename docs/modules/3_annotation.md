# Module 3: Cell Annotation & Differential Expression

!!! abstract "Overview"
    This module assigns biological identities to the clusters generated in Module 2. It supports two automated cell-type annotation frameworks: **CellTypist** (machine-learning based using logistic regression) and **ScType** (marker-gene score driven), as well as a **PreAnnotated** mode for pre-existing or scANVI annotations.
    
    After assigning identities, it performs comprehensive **Differential Expression (DE) Analysis** across clusters using Wilcoxon rank-sum tests to extract statistically significant marker genes, generating summary tables, dotplots, per-sample spatial diagnostic plots, and per-cluster tables for downstream web visualization.

---

## Parameters

### General Settings {#hide-me}
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chosen_cluster` | `String` | `"leiden_n10_r1.0"` | Which clustering resolution to use (e.g., `leiden_n10_r1.0`). If mode is set to `"All"`, this is used as fallback or ignored. Set to `"All"` in PreAnnotated mode to run DE on all clusterings. |
| `plot` | `Boolean` | `true` | Whether to generate diagnostic UMAP, spatial plots, cell-type composition bar charts, and dotplots. |

### CellTypist (Machine Learning) {#hide-me}
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `CellTypist_anno` | `Boolean` | `true` | Set to `true` to enable CellTypist annotation. |
| `CellTypist_mode` | `String` | `"All"` | `"one"` runs on `chosen_cluster` only. `"All"` loops through and annotates every clustering resolution generated in Module 2. |
| `CellTypist_model` | `String` | `"Cells_Intestinal_Tract"` | Name of the pre-trained model to download/load (e.g., `Immune_All_Low`, `Cells_Intestinal_Tract`). |
| `CellTypist_custom_model` | `String` | `""` | *(Optional)* Path to a locally saved custom `.pkl` model file. |
| `CellTypist_train` | `Boolean` | `false` | If `true`, trains a new CellTypist model using an external reference dataset. |
| `CellTypist_train_data` | `String` | `""` | *(Optional)* Path to reference `.h5ad` file for training. Required if `CellTypist_train = true`. |
| `CellTypist_train_labels` | `String` | `""` | *(Optional)* Column name in reference `.obs` containing ground truth labels. Required if `CellTypist_train = true`. |

### ScType (Marker Gene Driven) {#hide-me}
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ScType_anno` | `Boolean` | `true` | Set to `true` to enable ScType annotation. |
| `ScType_mode` | `String` | `"All"` | `"one"` runs on `chosen_cluster` only. `"All"` loops through and annotates every clustering resolution generated in Module 2. |
| `ScType_tissue` | `String` | `"Intestine"` | The specific tissue type from the ScType database (e.g., `Immune system`, `Intestine`, `Brain`). |
| `ScType_custom_db` | `String` | `""` | *(Optional)* Path to a custom `.xlsx` marker gene database file. |

### Differential Expression Strictness (`DE_params`) {#hide-me}
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pval_adj` | `Float` | `0.05` | Maximum Benjamini-Hochberg adjusted p-value for a gene to be considered significant. |
| `logfoldchange` | `Float` | `0.5` | Minimum Log2 Fold Change required for a gene to be considered a marker. |
| `min_expr_frac` | `Float` | `0.25` | Minimum fraction of cells inside the cluster expressing the gene (`pct_nz_group`). |
| `specificity_margin` | `Float` | `0.1` | Minimum difference in non-zero fraction between the target cluster and background (`pct_nz_group - pct_nz_reference`). |
| `junk_prefixes` | `List[String]` | `["MT-", "RPS", "RPL", "MALAT1"]` | Gene prefixes ignored during marker selection (e.g., mitochondrial, ribosomal). |

* **Safety Filters:** Clusters containing fewer than 2 cells are automatically excluded from Differential Expression testing to prevent matrix instability. Integer-based cluster names are automatically prefixed (e.g., `Cluster_0`) to ensure compatibility with downstream JSON parsers.

---

## Example Config

```toml
[modules.Annotate]
chosen_cluster = "leiden_n10_r1.0"
plot = true

# CellTypist Settings
CellTypist_anno = true
CellTypist_mode = "All"
CellTypist_model = "Cells_Intestinal_Tract"
CellTypist_custom_model = ""
CellTypist_train = false
CellTypist_train_data = ""
CellTypist_train_labels = ""

# ScType Settings
ScType_anno = true
ScType_mode = "All"
ScType_tissue = "Intestine"
ScType_custom_db = ""

[modules.Annotate.DE_params]
pval_adj = 0.05
logfoldchange = 0.5
min_expr_frac = 0.25
specificity_margin = 0.1
junk_prefixes = ["MT-", "RPS", "RPL", "MALAT1"]
```

---

## Outputs

### Folder Structure {#hide-me}

Outputs are saved in the `{analysis_name}_analysis/3_Annotate/` directory. Each method (`CellTypist`, `ScType`, or `PreAnnotated`) creates subdirectories per clustering resolution (`{cluster_col}`):

```text
3_Annotate/
├── adata.h5ad                                            # Updated AnnData object containing all annotations
└── {Method}/                                             # "CellTypist", "ScType", or "PreAnnotated"
    └── {cluster_col}/                                    # e.g., "leiden_n10_r1.0"
        ├── celltypecomp_{annotation_col}.png             # Stacked bar chart of predictions per cluster (CellTypist)
        ├── celltypecomp_{annotation_col}.txt             # Text breakdown of major/minor predicted cell types (CellTypist)
        ├── umapspatialscatter_{annotation_col}_{sample}.png  # Side-by-side Spatial scatter and UMAP (one per sample)
        └── DE_analysis/
            ├── markers_{annotation_col}.xlsx             # Filtered list of statistically significant marker genes
            ├── top_DEgenes_{annotation_col}.csv          # Summary table of the top 10 marker genes per cluster
            ├── dotplot_{annotation_col}.png              # Scanpy dot plot for top 5 marker genes per cluster
            └── DEgenes/                                  # Individual cluster files for Module 10 (Web Backend)
                ├── cluster_{cluster_name_1}_data.csv
                ├── cluster_{cluster_name_2}_data.csv
                └── ...
```

*(Note: `{annotation_col}` corresponds to the cell type label evaluated, e.g. `CellTypist_majorityvoting_leiden_n10_r1.0`, `sctype_leiden_n10_r1.0`, or `leiden_n10_r1.0`)*

---

### AnnData Updates {#hide-me}

Running Module 3 adds the following metadata and layers to `adata.h5ad`:

* **`adata.obs`**:
    * `'CellTypist_predictedlabels'`: Raw per-cell label predictions from CellTypist (if run).
    * `'CellTypist_confidence'`: Prediction confidence scores from CellTypist (if run).
    * `'CellTypist_majorityvoting_{cluster}'`: Harmonized cell type labels derived from majority voting across clusters.
    * `'sctype_{cluster}'`: Cell type assignments per cluster computed by ScType (if run).
    * `'named_{cluster}'`: Renamed cluster column (e.g., `Cluster_0`, `Cluster_1`) if original clusters were integer IDs.
* **`adata.obsm`**:
    * `'CellTypist_probabilities'`: Prediction probability matrix across all candidate reference classes (if CellTypist run).
* **`adata.uns`**:
    * `'CellTypist_probability_columns'`: List of cell type names corresponding to columns of `adata.obsm['CellTypist_probabilities']`.
    * `'rank_genes_groups'`: Wilcoxon rank-sum test metrics generated by Scanpy's `sc.tl.rank_genes_groups`.
    * `'{annotation_col}_colors'`: Hex color palettes generated dynamically for consistent visualization across plots.

## Core Libraries & Functions Used
