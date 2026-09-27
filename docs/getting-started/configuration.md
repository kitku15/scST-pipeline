# Configuration (TOML)

scSpatial-Kit is controlled entirely through a `config.toml` files (more about them [here](https://toml.io/en/)). Below is a breakdown of the three core configuration blocks. Module-specific configurable settings are detailed on their respective pages under the [Pipeline Modules](../../modules/0_format/) tab.

## 1. `[project]` {#hide-me}
Defines global metadata for your analysis.
```toml
[project]
analysis_name = "Liver_Cancer_Atlas"
data_type = "CosMx"        # Must be "CosMx" or "Xenium"
batch_key = "slide_id"     # Column name for slide/batch correction
sample_key = "sample_id"   # Column name for distinct biological samples
```

## 2. `[io]` {#hide-me}
Defines where your raw data lives and where intermediate/Zarr files should be stored.

### Batch Mode (Multiple Slides) {#hide-me}
If you have multiple slides, point `base_raw_dir` to the parent folder containing all your slide directories. The pipeline will automatically iterate through them.
```toml
[io]
base_raw_dir = "data/raw_data/"
base_zarr_dir = "data/zarr_stores/"
```

### Use Your Own Data {#hide-me}
If you are bypassing Module 0 and bringing an already processed `.h5ad` file, use `custom_input_adata`:
```toml
[io]
entry_point = "3" # e.g., Start at Cell Annotation
custom_input_adata = "data/my_preprocessed_dataset.h5ad"
```

## 3. `[pipeline]` {#hide-me}
Defines exactly which modules the pipeline should execute. You can comment out modules you wish to skip.  <mark>check this again</mark>

```toml
[pipeline]
modules = [
    "0_FormatData",
    "1_QualityControl",
    "1b_MergeData",
    "2_DimensionReduction",
    "3_Annotate",
    "4_ViewImages",
    "5_SpatialStat",
    "6_MuSpan",
    "7_Decoupler",
    "8_Cellphonedb",
    "8b_LIANA",
    "8c_LIANA_Causal",
    "9_DEAnalysis",
    "10_WebVisPrep"
]
```

!!! tip "Partial Execution"
    You don't have to run everything at once! A common workflow is to run Modules 0, 1, and 1b first to check QC metrics, and then modify the TOML to run Modules 2 through 10 once you are happy with the cell filtering thresholds.
