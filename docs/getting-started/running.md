# Running the Pipeline

scSpatial-Kit supports two modes of execution: **Snakemake** for automated workflow management (`run_pipeline.py`), and a Direct Python Orchestrator (`__main__.py`) for manual HPC cluster job arrays.

## Basic Execution
Once your `config.toml` is set up, run the pipeline by passing the config file to the orchestrator:

```bash
python run_pipeline.py my_config.toml --cores 16
```
* **`--cores`**: Defines the maximum number of CPU threads the pipeline is allowed to use. 

## Dry Run Mode
Always use the `--dry-run` flag before running a huge job. This tells the pipeline to calculate the execution graph and check for missing files without actually executing the data processing steps.

```bash
python run_pipeline.py my_config.toml --dry-run
```

## Advanced Execution
If you are running the pipeline on an HPC cluster (like PBS, Slurm, or LSF) and prefer to bypass Snakemake to submit parallel array jobs manually, you can interact directly with the pipeline's core entry point: `__main__.py`.

```bash
python Modules/__main__.py my_config.toml --modules 1 2 3
```

### Useful CLI Overrides
The `__main__.py` script accepts several overrides that are particularly useful for cluster computing:

| Flag | Example | Description |
|------|---------|-------------|
| `--modules` | `--modules 1 1b 2` | Overrides the `[pipeline.modules]` list in the TOML, allowing you to force specific steps to run. |
| `--sample_index` | `--sample_index 2` | *(Useful for Slurm/PBS Arrays)* Tells the pipeline to only process the $N^{th}$ slide folder in your `base_raw_dir`. |
| `--sample_name` | `--sample_name "Slide_A"` | Manually targets a specific sample for processing (useful for re-running a failed sample). |
| `--comp_index` | `--comp_index 0` | Isolates a specific Differential Expression comparison for parallel processing. |

 <mark>check above again? </mark>

## Runtime Tracking & Logs
Every time the pipeline runs, it tracks its execution time and success/failure states.

* **Console Logs:** Detailed progress logs are printed to your terminal and saved to `analysis_{name}/logs/scSpatial-Kit-{date}.log`.
* **Runtime Profiling:** The duration of every module is recorded in `analysis_{name}/logs/pipeline_runtimes.csv`.  <mark>fix this / check again</mark>
