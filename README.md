# General Single-cell Spatial Transcriptomics Pipeline

An end-to-end Python pipeline for the processing, analysis, and visualization of Spatial Transcriptomics data. Built specifically to handle **NanoString CosMx** and **10x Genomics Xenium** datasets, this pipeline uses modern spatial data frameworks (spatialdata, scanpy, squidpy), advanced spatial statistics (muspan), downstream TF/CCC Analysis, and Causal Network Inference. 

More visual indepth explanation of the pipeline presented [here](https://docs.google.com/presentation/d/1k17lRxf43-NQRZWIfLpKvocqIi1PN-YfwAcD-urLCCc/edit?usp=sharing). 

## 🔑 Key Features

- Non linear Dimension reduction with [scVI](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scvi.html#), [scANVI](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scanvi.html), and [scVIVA](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scviva.html)
- Supports both Machine Learning-based  ([CellTypist](https://www.celltypist.org/)) and Marker-based ([ScType](https://github.com/kris-nader/sc-type-py)) cell type annotation.
- Generates Delaunay, KNN, and Proximity graphs, along with cross-Pair Correlation Functions (PCF), and cell morphological metrics via the [SquidPy](https://squidpy.readthedocs.io/en/stable) and [MuSpAn](https://www.muspan.co.uk/) library.
- Transcription factor activity inference via [DecoupleR](https://decoupler.readthedocs.io/en/latest/index.html) and targeted pseudobulk Differential Expression using PyDESeq2.
- Cell-Cell Communication (CCC) & Causal Networks - Broad microenvironment CCC via [CellphoneDB](https://cellphonedb.readthedocs.io/en/latest), continuous spatial CCC via LIANA+, and downstream causal signaling cascade inference via Corneto.
- Prepares and packages all results (Zarr, JSON) into a `.tar` archive for local exploration using the [Spatial-VisKit](https://github.com/kitku15/Spatial-VisKit) React web application.
- Run pipeline reproducibly on HPC clusters (or locally) using Apptainer/Singularity containers, controlled via TOML configuration files.

## ⭐ Pipeline Structure
The workflow is divided into two main execution phases to allow slide-specific quality control before integration:

### 1️⃣ Data Formating and Quality Control 
Modules:
1. `0_format`
2. `1_QualityControl`

*Details:* Uses one config file per slide/batch so you can adjust QC thresholds independently. Supports parallel processing on HPC array jobs.

### 2️⃣ Downstream Analysis
Modules:

- `1b_MergeData`: Merges all slides into a single AnnData object.
- `2_DimensionReduction`: scVI, scANVI, scVIVA, UMAP, and Leiden clustering.
- `3_Annotate`: CellTypist, ScType, and marker gene identification.
- `4_ViewImages`: Spatial embedding plots, ROI grid generation, and gene expression mapping.
- `5_SpatialStat`: Squidpy (Centrality Scores, Co-occurrence, Neighborhood Enrichment, Moran's I).
- `6_MuSpan`: Cross-PCF, Graph Construction, Shape/Morphology, and Cell Proximity Analysis.
- `7_Decoupler`: TF activity inference (CollecTRI/DoRothEA) and pseudobulk generation.
- `8_Cellphonedb`: Standard microenvironment-based CCC.
- `8b_LIANA`: Bivariate spatial CCC and NMF signaling programs.
- `8c_LIANA_Causal`: Condition-specific causal networks connecting LR pairs to TFs (via Corneto).
- `9_DEAnalysis`: Targeted pairwise pseudobulk differential expression.
- `10_WebVis`: Compiles and optimizes outputs (Zarr, JSON) for Spatial-VisKit.

Additional:
- Merges all slides into one anndata object
- One config file
- slide/batch correction handled by scVI and scVIVA

## Configuration Files
The pipeline is driven by `.toml` files. You need two main types of configurations:

### Config File 1: Data Formatting and Quality Control
File name: `config_myanalysis{slide_no}.toml`. Each slide/ batch needs to have one of these as you want to be able to adjust QC thresholds according to each slide. Example:

Filename: `config_tyler1.toml`
```toml
log_level = "INFO"
seed = 21122023

[project]
analysis_name = "tyler_project" # analysis/output folder
slide_name = "Slide_1" # or Slide_2, 3, 4 
data_type = "CosMx" # or Xenium

[io]
base_dir = "."
dataset_id = "TylerWooldridgeSlide1" # only used for CosMx dataset 
dataset_dir = "data/tyler_1" # store your dataset folder inside the data dir
zarr_dir = "data/tyler_1.zarr" # keep this the same name as your dataset name

[pipeline] # dont change anything here 
modules = ["0_format", "1_QualityControl"]

[modules.QualityControl]
min_counts = 50 # threshold- minimum gene count for a cell 
min_cells = 5 # threshold- minimum cell count for a gene
min_dapi = 100 # threshold- minimum DAPI for a cell
fov_metadata_path = "data/tyler_1/TylerWooldridgeSlide1_fov_cohort.csv" # metadata that contains columns to add to adata.obs (will need to make this customizable but currently expects these columns: DiseaseType, TreatmentResponse, sample_id) If None, no additional metadata collumns will be added. 

```

### Config File 2: Downstream Analysis
File name: `config_downstream.toml`. Only one configuration file for your merged dataset (merged all slides).
```toml
log_level = "INFO"
seed = 21122023

[project]
analysis_name = "tyler_project" # THIS MUST BE THE SAME AS WHAT YOU SET IN 01 CONFIG
data_type = "CosMx" # or Xenium
batch_key = "slide_id" 
sample_key = "sample_id"

[io]
base_dir = "."

[io.raw_data]
"Slide_1" = { dataset_dir = "data/tyler_1", zarr_dir = "data/tyler_1.zarr", proseg_zarr_dir = "data/tyler1_proseg.zarr" }
"Slide_2" = { dataset_dir = "data/tyler_2", zarr_dir = "data/tyler_2.zarr", proseg_zarr_dir = "data/tyler2_proseg.zarr" }

[pipeline]
modules = [
    "1b_MergeData", "2_DimensionReduction", "3_Annotate", "4_ViewImages", 
    "5_SpatialStat", "6_MuSpan", "7_Decoupler", "8_Cellphonedb", 
    "8b_LIANA", "8c_LIANA_Causal", "9_DEAnalysis", "10_WebVis"
    ]

[modules.MergeData] # details your individual slides to be merged
input_files = [
    "analysis_tyler_project/1_QualityControl/Slide_1/adata.h5ad",
    "analysis_tyler_project/1_QualityControl/Slide_2/adata.h5ad",
    "analysis_tyler_project/1_QualityControl/Slide_3/adata.h5ad",
    "analysis_tyler_project/1_QualityControl/Slide_4/adata.h5ad"
]
slide_names = ["Slide_1", "Slide_2", "Slide_3", "Slide_4"] # everything you set for slide_name in all your slide configs

[modules.DimensionReduction]
use_scviva = false
n_neighbors = [5, 10, 30, 50] # list of no. of neighbours parameter for neighbourhood graph   
resolution = [0.5, 1.0, 1.5, 2.0] # list of resolution parameter for leiden clustering 
cluster_name = "leiden" # dont change this
run_pca = true # do you want to see a PCA plot? (true/false)
n_comps = 50 # if yes how many PCs do you want plotted?
cluster_name = "leidenpca"
run_pca = true
n_comps = 50     
umap_latent = "X_pca"
scvi_epochs = 400

scviva_layer = "counts"
scviva_batch_key = "slide_id" 
scviva_spatial_knn = 10
scviva_epochs = 400

[modules.Annotate]
chosen_cluster = "leiden_n10_r1.0" # your favorite clustering column from the DR module

# ScType Config:
ScType_anno = true # do you want to run ScType? (true/false)
ScType_mode = "All" # do you want ScType Annotation to be done for all ("All") clustering columns or just the one you chosen_cluster? ("one")
ScType_custom_db = "" # you can insert your own marker db here 
ScType_tissue = "Intestine" # ScType Tissue - look inside ScType_db/ScTypeDB_full.xlsx

# CellTypist Config:
CellTypist_anno = true # do you want to run CellTypist? (true/false)
CellTypist_mode = "All" # do you want CellTypist Annotation to be done for all ("All") clustering columns or just the one you chosen_cluster? ("one")
CellTypist_model = "Cells_Intestinal_Tract" # CellTypist model of choice 

[modules.ViewImages]
# list of genes you want to view spatially on each sample
gene_list = ["HLA-DRB", "KRT8", "TMSB4X", "TFF3", "IGHA1", "COX1"]

[modules.MuSpan]
selection_name = "myMuSpan" # make a selection name for MuSpan
cell_types = [0, 1] # select 2 cell types you want analysis on

# list of genes to visualize on MuSpan plots 
transcripts = ["HLA-DRB", "KRT8", "TMSB4X", "TFF3", "IGHA1", "COX1"]

# MuSpAn Spatial Graphs config
min_edge_distance = 0
max_edge_distance = 200
distance_list = [50, 90, 120]
min_edge_distance_shape = 0
max_edge_distance_shape = 1
k_list = [2, 5, 10, 15]

# MuSpAn Cell Proximity Analysis config
cellboundary_label = 'Cell boundaries' 
max_distance = 200 # Threshold for contact
network_celltypes = ["0", "1", "2", "3"]

[modules.Decoupler]
organism = "human" # or mouse
grn = "dorothea" # choice of GRN db to use 'dorothea' or 'collectri'
dorothea_levels = ["A", "B", "C"] 

[modules.Cellphonedb]
no_microenv = 20 # how many physical microenvironments to split the slide on
chosen_cluster = 'CellTypist_majorityvoting_leiden_n10_r1.0' # cell type annotation column to use
cpdb_version = 'v5.0.0'
human = true # if false, mouse is assummed 
celltypes = ["IgA plasma cell", "BEST2+ Goblet cell"] # cell types you want dot and chord plots for


[modules.LIANA_Causal]
celltype_col = "Final_Annotation"
treatment_col = "TreatmentResponse"
comparisons = [["CPIc", "Healthy"], ["IFX_NR", "Healthy"]]
cell_type_pairs = [["Macrophage", "CD8 T cell"], ["B Cells", "T Cells"]]

[modules.DEAnalysis]
celltype_col = "Final_Annotation"
treatment_col = "TreatmentResponse"
comparisons = [["IFX_NR", "IFX_R"], ["CPIc", "Healthy"]]

[modules.WebVisPrep]
primary_annotation = "Final_Annotation"
microenv_col = "spatial_microenvironment"
annotation_columns = ["leiden", "Final_Annotation", "Broad_Celltype"]

```
## ⭐ Inputs
These dataset folders need to be in the data directory
### 🔘 Xenium 10X

Example below showing one slide. Need to refine.

```toml
# check this below if everything is really needed
MyXeniumDataset/
├── morphology_focus/
├── analysis_summary.html
├── analysis.tar.gz
├── analysis.zarr.zip
├── aux_outputs.tar.gz
├── cell_boundaries.csv.gz
├── cell_boundaries.parquet
├── cell_feature_matrix.h5
├── cell_feature_matrix.tar.gz
├── cell_feature_matrix.zarr.zip
├── cells.csv.gz
├── cells.parquet
├── cells.zarr.zip
├── experiment.xenium
├── gene_panel.json
├── metrics_summary.csv
├── morphology.ome.tif
├── nucleus_boundaries.csv.gz
├── nucleus_boundaries.parquet
├── transcripts.csv.gz
├── transcripts.parquet
└── transcripts.zarr.zip
```

### 🔘 CosMx NanoString
Example below showing 2 slides. `dataset_id` (for 1st config file) for MyCosMxDataset_1 is MyCosMxDatasetSlide1. `dataset_id` for MyCosMxDataset_2 is MyCosMxDatasetSlide2.

```toml
MyCosMxDataset_1/
├── CellComposite/ # these could be empty 
├── CellLabels/ # these could be empty 
├── MyCosMxDatasetSlide1_exprMat_file.csv
├── MyCosMxDatasetSlide1_fov_cohort.csv
├── MyCosMxDatasetSlide1_fov_positions_file.csv
├── MyCosMxDatasetSlide1_metadata_file.csv
├── MyCosMxDatasetSlide1_tx_file.csv
└── MyCosMxDatasetSlide1-polygons.csv
MyCosMxDataset_2/
├── CellComposite/ # these could be empty 
├── CellLabels/ # these could be empty 
├── MyCosMxDatasetSlide2_exprMat_file.csv
├── MyCosMxDatasetSlide2_fov_cohort.csv
├── MyCosMxDatasetSlide2_fov_positions_file.csv
├── MyCosMxDatasetSlide2_metadata_file.csv
├── MyCosMxDatasetSlide2_tx_file.csv
└── MyCosMxDatasetSlide2-polygons.csv
```
⚠️ **WARNING!**: make sure your `*metadata_file.csv` does not have both `cell_id` and `cell_ID` or any other case variation of `cell_id`. If this is the case, a simple python script can just rename one of them. 

## ⭐ Getting Started

### 1. Run Data Formatting and QC

Adjust settings in your job script (`01_run_qc.sh`), and set parameters in your config file. The script leverages PBS arrays to process multiple slides in parallel.
```bash
# Inside 01_run_qc.sh
CONFIG_NAME="config_tyler${PBS_ARRAY_INDEX}.toml" # path to your 1st config files
SIF_IMAGE="$(readlink -f kitku.sif)" # singularity image

apptainer run --writable-tmpfs -W "$APPTAINER_WORKDIR" \
  --bind "$PBS_O_WORKDIR:/app" \
  "$SIF_IMAGE" \
  "$CONFIG_NAME" --modules 0 1
```
Then to queue the job, run on the command line:
```bash
qsub 01_run_qc.sh
```
### 2. Run Downstream Analysis

Once the QC jobs finish, configure `config_downstream.toml` and submit the downstream script (`02_run_downstream.sh`). You can run all modules at once, or incrementally. 

```bash
# Inside 02_run_downstream.sh
CONFIG_NAME="config_downstream.toml" # path to your 2nd config file
SIF_IMAGE="$(readlink -f kitku.sif)" # singularity image

# You can pass multiple modules separated by space:
apptainer run --writable-tmpfs -W "$APPTAINER_WORKDIR" \
  --bind "$PBS_O_WORKDIR:/app" \
  "$SIF_IMAGE" \
  "$CONFIG_NAME" --modules 1b 2 3 4 5 6 7 8 8b 8c 9 10
```
Then to queue the job, run on the command line:
```bash
qsub 02_run_downstream.sh
```
## ⭐ Outputs
The pipeline generates a highly organized output directory. Module 10 packages the necessary files into a `.tar` archive optimized for the [Spatial-VisKit](https://github.com/kitku15/Spatial-VisKit) Web App.
```text
{analysis_name}_analysis/
├── 1_QualityControl/
│   ├── Slide_1/             # QC plots, histograms, and filtered adata.h5ad
│   └── Slide_2/  
├── 1b_MergeData/            # Merged dataset containing all slides
├── 2_DimensionReduction/    # PCA, UMAP, and scVIVA embeddings
├── 3_Annotate/              # CellTypist/ScType annotations & marker genes
├── 4_ViewImages/            # Spatial tissue plots and ROI grid definitions
├── 5_SpatialStat/           # Squidpy neighborhood enrichments & Moran's I
├── 6_MuSpan/                # MuSpAn graphs, PCF, and Morphometrics per ROI
├── 7_Decoupler/             # TF activities and Pseudobulk generation
├── 8_Cellphonedb/           # CellPhoneDB dot/chord plots
├── 8b_LIANA/                # Bivariate spatial CCC & NMF components
├── 8c_LIANA_Causal/         # Corneto causal networks & pathway CSVs
├── 9_DEAnalysis/            # PyDESeq2 condition-specific DE outputs
├── 10_WebVisPrep/           # Web-optimized JSONs, Vitessce Zarrs, and Sankeys
│
└── logs/                    # Comprehensive logs and runtime trackers (.csv)
```


## ⭐ About
This 22-week project was completed as part of the MRes in Bioinformatics and Theoretical Systems Biology at Imperial College London. It builds upon the [The ReCoDe-spatial-transcriptomics repository](https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics). The work was supervised by Tamas Korcsmaros and Balazs Bohar.

- 😸 Author: Bunga Tiasyaira Hutasuhut (Syaii) 
- 📩 Academic Email: bth22@ic.ac.uk
- 📮 Personal Email: bungatiasyaira@outlook.com

