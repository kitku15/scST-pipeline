# tests/test_annotate.py
import pandas as pd
import scanpy as sc
import numpy as np
from Annotate import run_annotate


def test_run_annotate(dummy_adata, tmp_path):
    # 1. SETUP: Pretend Module 2 just finished.
    # Add a fake 'leiden' cluster column and save the dummy data.
    dummy_adata.obs["leiden"] = np.random.choice(
        ["0", "1", "2"], size=dummy_adata.n_obs
    )

    prev_module_dir = tmp_path / "Module2"
    prev_module_dir.mkdir()
    dummy_adata.write_h5ad(prev_module_dir / "adata.h5ad")

    out_dir = tmp_path / "Module3"

    # 2. ACTION: Run the annotation module
    run_annotate(
        module_dir=out_dir,
        module_name="3_Annotate",
        cluster_name="leiden",
        new_clusters="cell_type",
        prev_module_dir=prev_module_dir,
    )

    # 3. ASSERTIONS
    # Check if files were created
    assert (out_dir / "markers.xlsx").exists(), "Markers Excel file missing"
    assert (out_dir / "top_differentially_expressed_genes.csv").exists(), (
        "Top genes CSV missing"
    )

    # Check if cluster specific CSVs were created
    cluster_0_csv = out_dir / "cluster_diff_genes" / "cluster_0_data.csv"
    assert cluster_0_csv.exists(), "Cluster specific CSV missing"

    # Check the contents of the generated CSV
    top_genes_df = pd.read_csv(out_dir / "top_differentially_expressed_genes.csv")
    assert "Cluster Number" in top_genes_df.columns, (
        "CSV is missing the Cluster Number column"
    )
    assert "Top Genes" in top_genes_df.columns, "CSV is missing the Top Genes column"

    # Check if the AnnData object was updated and saved
    processed_adata = sc.read_h5ad(out_dir / "adata.h5ad")
    assert "cell_type" in processed_adata.obs.columns, (
        "New cluster name column not added to adata"
    )

    # Verify the mapping worked (e.g., '0' became 'Cluster_0')
    assert "Cluster_0" in processed_adata.obs["cell_type"].values
