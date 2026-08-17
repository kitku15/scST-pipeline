# General Single-cell Spatial Transcriptomics Pipeline

An end-to-end Python pipeline for the processing, analysis, and visualization of Spatial Transcriptomics data. Built specifically to handle **NanoString CosMx** and **10x Genomics Xenium** datasets, this pipeline uses modern spatial data frameworks (spatialdata, scanpy, squidpy), advanced spatial statistics (muspan), and downstream TF and CCC Analysis. More visual indepth explanation of the pipeline presented [here](https://docs.google.com/presentation/d/1k17lRxf43-NQRZWIfLpKvocqIi1PN-YfwAcD-urLCCc/edit?usp=sharing). 

## 🔑 Key Features

- Non linear Dimension reduction with [scVI](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scvi.html#) and [scVIVA](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scviva.html)
- Supports both Machine Learning-based  ([CellTypist](https://www.celltypist.org/)) and Marker-based ([ScType](https://github.com/kris-nader/sc-type-py)) cell type annotation.
- Generates Delaunay, KNN, and Proximity graphs, along with cross-Pair Correlation Functions (PCF) via the [SquidPy](https://squidpy.readthedocs.io/en/stable) and [MuSpAn](https://www.muspan.co.uk/) library.
- Transcription factor analysis with [DecoupleR](https://decoupler.readthedocs.io/en/latest/index.html)
- Cell-Cell communication analysis with [CellphoneDB](https://cellphonedb.readthedocs.io/en/latest)
- Reproducible runs controlled entirely via a several TOML configuration files.
- Run on HPC with Singularity/Apptainer! (works nicely :D)
- In development: a web visualization tool ([Spatial-VisKit](https://github.com/kitku15/Spatial-VisKit)) to interact and explore outputs of this pipeline

## ⭐ Pipeline Structure
Below is a high level overview of what the pipeline does. It is divided into 2 parts:
### 1️⃣ Data Formating and Quality Control 
Modules:
1. **Format Data**
2. **Quality Control**

Additional:
- One config file each batch / slide 
- Process multiple slides in parallel 

### 2️⃣ Downstream Analysis
Modules:

3. **Dimension Reduction**: scVI, scVIVA, UMAP, Leiden clustering
4. **Cell Type Annotation**: CellTypist, ScType, DE Analysis
5. **Visualization**: Gene visualization
6. **Spatial Stats (Squidpy)**: Centrality Scores, Co-occurrence Probability, Neighborhood Enrichment, Moran's I 
7. **Spatial Stats (MuSpan)**: Cross-Pair Correlation Function, Graph Construction, Shape & Morphology Analysis, Cell-Cell Proximity & Neighborhoods
8. **Transcription Factor Analysis (Decoupler)**: CollecTRI or DoRothEA Gene Regulatory Network database 
9. **Cell-Cell Communication Analysis (CellPhoneDB)**: Split by microenvironment

Additional:
- Merges all slides into one anndata object
- One config file
- slide/batch correction handled by scVI and scVIVA

## Configuration Files
To configure analysis parameters specific to a dataset, sample, or slide, we use a `config.toml` file. For the pipeline, you need 2 main types of config files. 

### Config File 1: Data Formating and Quality Control
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

[io]
base_dir = "."

[pipeline]
modules = [ # dont change anything here 
    "1b_MergeData",
    "2_DimensionReduction",
    "3_Annotate",
    "4_ViewImages",
    "5_SpatialStat",
    "6_MuSpan",
    "7_Decoupler",
    "8_Cellphonedb"
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
n_neighbors = [5, 10, 30, 50] # list of no. of neighbours parameter for neighbourhood graph   
resolution = [0.5, 1.0, 1.5, 2.0] # list of resolution parameter for leiden clustering 
cluster_name = "leiden" # dont change this
run_pca = true # do you want to see a PCA plot? (true/false)
n_comps = 50 # if yes how many PCs do you want plotted?


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

How to get started with the pipeline on the HPC with Singularity. First, we run the first part: Data Formating and Quality Control. Adjust settings in `01_run_qc.sh`, set parameters in your config file. 
```bash
CONFIG_NAME="config_tyler${PBS_ARRAY_INDEX}.toml" # path to your 1st config files
SIF_IMAGE="$(readlink -f kitku.sif)" # singularity image
```
Then to queue the job, run on the command line:
```bash
qsub 01_run_qc.sh
```
When that job finishes, we will run the second part: Downstream Analysis. Adjust settings in `02_run_downstream.sh`, set parameters in your config file.
```bash
CONFIG_NAME="config_downstream.toml" # path to your 2nd config file
SIF_IMAGE="$(readlink -f kitku.sif)" # singularity image

# find this section in the file
# this is where you can select which modules to run 
# you can run all, or one at a time, but not a later module without running all earlier ones first!
"$CONFIG_NAME" --modules 1b 2 3 4 5 6 7 8
```
Then to queue the job, run on the command line:
```bash
qsub 02_run_downstream.sh
```
## ⭐ Outputs
Outputs will be structured like shown below.

```bash
{analysis_name}_analysis/
├── 1_QualityControl/
    ├── Slide_1/ # QC plots and adata for each slide 
    └── Slide 2/  
├── 1b_MergedData/
    └── merged adata.h5ad # merged adata
├── 2_DimensionReduction/
    ├── n10_r0.1/ # plots for each DR parameter
    └── n10_r0.5/
├── 3_Annotate/
    ├──CellTypist/ # plots for each CT annotation for each DR parameter
    └──ScType/
├── 4_ViewImages/ 
├── 5_SpatialStat/
├── 6_MuSpan/
    └── selection_name/ 
├── 7_Decoupler/
├── 8_Cellphonedb/
    └── CellType A # for cell type specific dot & chord plots
    └── CellType B
    └── cpdb_out # cellphonedb outputs
└── logs/
```


## ⭐ Additional
This 22-week project was completed as part of the MRes in Bioinformatics and Theoretical Systems Biology at Imperial College London. It builds upon the [The ReCoDe-spatial-transcriptomics repository](https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics). The work was supervised by Tamas Korcsmaros and Balazs Bohar.

- 😸 Author: Bunga Tiasyaira Hutasuhut (Syaii) 
- 📩 Academic Email: bth22@ic.ac.uk
- 📮 Personal Email: bungatiasyaira@outlook.com

