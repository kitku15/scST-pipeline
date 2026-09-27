# Module 8: Cell-Cell Communication & Signaling Networks

!!! abstract "Overview"
    This module maps ligand-receptor (LR) cross-talk across the tissue using three complementary approaches:
    
    * **Module 8 (CellPhoneDB v5)**: Evaluates pairwise LR interactions constrained within physically defined spatial microenvironments, with optional downstream transcription factor validation via the CellSign module.
    * **Module 8b (LIANA+ Spatial)**: Computes continuous spatial bivariate cross-correlations (spatially-weighted expression decay) between ligands and receptors, and uses Non-Negative Matrix Factorization (NMF) to extract coordinate-anchored communication programs.
    * **Module 8c (LIANA+ Causal via Corneto)**: Constructs condition-specific, directed intracellular signaling cascades (Ligand $\rightarrow$ Receptor $\rightarrow$ Kinase $\rightarrow$ TF) by formulating an optimization problem on the OmniPath protein-protein interactome.

## Parameters

### Module 8: CellPhoneDB Settings (`[modules.Cellphonedb]`) {#hide-me}

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chosen_cluster` | `String` | `"Final_Annotation"` | The cell type column in `adata.obs` to analyze interactions between. |
| `microenv_resolution` | `Float` | `0.5` | Leiden resolution used to cluster the tissue into physical spatial microenvironments prior to running CPDB. |
| `target_sample` | `String` | `""` | Target a specific sample ID to isolate for detailed plotting. |
| `celltypes` | `List[String]` | `[]` | Subset of cell types from `chosen_cluster` to generate dedicated dot and chord plots for. |
| `target_genes` | `List[String]` | `[]` | Limits analysis to candidate ligands and receptors of interest. |
| `gene_family` | `String` or `List` | `""` | Filter interactions by functional functional family (e.g., `"chemokines"`). |
| `cpdb_version` | `String` | `"v5.0.0"` | CellPhoneDB database version release tag. |
| `human` | `Boolean` | `true` | If `true`, assumes human gene nomenclature. If `false`, mouse symbols are converted to human orthologs. |
| `use_cellsign` | `Boolean` | `false` | Enables the CellSign module. Requires active TF inputs from Module 7 to validate receptor-to-TF downstream compatibility. |
| `target_microenv` | `List[String]` | `[]` | Restricts output plotting to designated physical microenvironments. |

### Module 8b: LIANA+ Spatial Settings (`[modules.LIANA]`) {#hide-me}

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bandwidth` | `Integer` | `400` | Gaussian distance decay bandwidth ($\mu\text{m}$) used to construct spatial interaction weight matrices. |
| `cutoff` | `Float` | `0.1` | Weight threshold below which spatial neighbor influence drops to zero. |
| `n_nmf_components` | `Integer` | `5` | Number of Non-negative Matrix Factorization (NMF) latent factors to deconvolve into distinct spatial signaling programs. |
| `nz_prop` | `Float` | `0.05` | Minimum proportion of non-zero expressing cells required for an LR pair to be evaluated. |
| `resource_name` | `String` | `"consensus"` | Ligand-receptor interaction database resource pulled from OmniPath. |

### Module 8c: LIANA+ Causal Settings (`[modules.LIANA_Causal]`) {#hide-me}

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `treatment_col` | `String` | `"TreatmentResponse"` | Column in `adata.obs` defining biological conditions/cohorts. |
| `celltype_col` | `String` | `"Final_Annotation"` | Cell identity column used for defining interacting pairs. |
| `comparisons` | `List[List[String]]` | `[]` | Condition contrasts formatted as `["Test_Condition", "Reference_Condition"]`. |
| `cell_type_pairs` | `List[List[String]]` | `[]` | Directed communication pairs to model, formatted as `["Source_Cell", "Target_Cell"]`. |
| `solver` | `String` | `"scipy"` | Optimization solver library used by Corneto (e.g., `"HIGHS"`, `"SCIPY"`, `"CLARABEL"`). |
| `num_top_receptors` | `Integer` | `10` | Limits causal inference to top $N$ differentially active receptors to restrict network size. |
| `num_top_tfs` | `Integer` | `5` | Limits causal inference to top $N$ active transcription factors derived from Module 7. |
| `curation_effort` | `Integer` | `5` | Minimum literature citation cutoff for OmniPath protein-protein interaction edges. |
| `expr_prop` | `Float` | `0.1` | Minimum fraction of expressing cells per population for a node to be included. |
| `node_cutoff` | `Float` | `0.1` | Node penalty score threshold applied during graph optimization. |
| `max_runs` | `Integer` | `50` | Maximum optimization cycles Corneto executes to find an optimal network topology. |
| `stable_runs` | `Integer` | `10` | Early-stopping condition: required number of identical network solutions. |

* **Note on Internal DEA:** Before inferring the causal network, Module 8c automatically executes an internal PyDESeq2 Differential Expression workflow on the source/target cell pairs, applying LFC shrinkage and filtering out genes with less than 10 total counts across at least 5 expressing cells.

## Example Config

```toml
# --- Module 8: CellPhoneDB ---
[modules.Cellphonedb]
chosen_cluster = "CellTypist_majorityvoting_leiden_n10_r1.0"
microenv_resolution = 0.5
target_sample = "Slide_1_Healthy_1"
celltypes = ["Macrophage", "CD8 T cell"]
target_genes = []
gene_family = ""
cpdb_version = "v5.0.0"
human = true
use_cellsign = false
target_microenv = []

# --- Module 8b: LIANA+ Spatial ---
[modules.LIANA]
bandwidth = 200
cutoff = 0.1
n_nmf_components = 5
nz_prop = 0.05
resource_name = "consensus"

# --- Module 8c: LIANA+ Causal ---
[modules.LIANA_Causal]
celltype_col = "CellTypist_majorityvoting_leiden_n10_r1.0"
treatment_col = "TreatmentResponse"
comparisons = [["CPIc", "Healthy"]]
cell_type_pairs = [["Macrophage", "CD8 T cell"]]
solver = "HIGHS"
num_top_receptors = 10
num_top_tfs = 10
curation_effort = 5
expr_prop = 0.1
node_cutoff = 0.1
max_runs = 20
stable_runs = 10
```

## Outputs

Outputs are saved under the respective module directories:
<mark>Need to apply this folder structure: </mark>

```text
{analysis_name}_analysis/
├── 8_Cellphonedb/
│   ├── adata.h5ad                                       # Updated AnnData with physical microenvironment assignments
│   ├── cpdb_out/                                        # Raw CellPhoneDB outputs
│   │   ├── deconvoluted.txt                             # Mean expression per complex subunit
│   │   ├── means.txt                                    # Uncorrected average interaction strengths
│   │   ├── pvalues.txt                                  # Empirical permutation test significance
│   │   └── significant_means.txt                        # Interaction scores for significant pairs (p < 0.05)
│   ├── chord_plots/
│   │   └── chordplot_{celltype}.png                     # Circos communication diagrams per selected cell type
│   └── dot_plots/
│       └── dotplot_{celltype}.png                       # Dot plot of significant LR pairs across microenvironments
│
├── 8b_LIANA/
│   ├── lrdata.h5ad                                      # Spatially-weighted ligand-receptor cross-correlation object
│   ├── nmf_adata.h5ad                                   # Deconvolved NMF signaling factors matrix object
│   ├── bdata_tf_lr.h5ad                                 # Generalized cross-talk object (TF vs LR pairs)
│   └── plots/
│       └── NMF_spatial_factors_{sample_id}.png          # Spatial maps showing active communication programs
│
└── 8c_LIANA_Causal/
    ├── {Contrast_Name}/
    │   ├── causal_net_{contrast}_{source}_to_{target}.csv # Directed molecular cascade nodes
    │   ├── dea_{contrast}.csv                           # Condition-specific differential expression results
    │   ├── liana_lr_{contrast}.csv                      # Condition-specific ligand-receptor interactions
    │   ├── tf_estimates_{contrast}.csv                  # Inferred downstream TF activities
    │   └── plots/
    │       ├── interaction_stat_histogram.png           # Wald statistic distributions
    │       └── liana_tileplot.png                       # Top LR interaction significance tileplot
```

## Core Libraries & Functions Used
* **[CellPhoneDB](https://github.com/ventolab/CellphoneDB)**:
    * `cpdb_statistical_analysis_method.call()`: Permutation-based scoring of multi-subunit ligand-receptor pairings constrained within spatial microenvironments.
* **[ktplots-py](https://github.com/zktuong/ktplots-py)**:
    * `kpy.plot_cpdb_chord()` & `kpy.plot_cpdb()`: Renders chord diagrams and dot plots for interaction results.
* **[LIANA+](https://liana-py.readthedocs.io/en/latest/)**:
    * **Sample-Aware Graphs:** The pipeline utilizes a custom wrapper to build spatial connectivity matrices *sample-by-sample*. This prevents artificial "ghost" communication edges from forming across different physical tissue slices stored in the same AnnData object.
    * `li.mt.bivariate()`: Evaluates spatially lagged bivariate correlation between partner proteins across continuous tissue coordinates.
    * `li.mt.nmf()`: Matrix decomposition identifying spatial co-expression programs.
* **[Corneto](https://github.com/saezlab/corneto)**:
    * `corneto.methods`: Solves the integer linear programming problem to infer optimal, parsimonious causal signaling paths connecting ligand-bound receptors to active transcription factors.