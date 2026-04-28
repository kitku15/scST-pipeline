# tests/test_spatial.py
import scanpy as sc
import numpy as np
from ViewImages import run_view_images
from SpatialStat import run_spatial_statistics


def test_view_images(dummy_adata, tmp_path):
    # Setup
    dummy_adata.obs["leiden"] = np.random.choice(["0", "1"], size=dummy_adata.n_obs)
    prev_dir = tmp_path / "Module3"
    prev_dir.mkdir()
    dummy_adata.write_h5ad(prev_dir / "adata.h5ad")

    out_dir = tmp_path / "Module4"

    # Action
    gene_list = [dummy_adata.var_names[0]]  # Just use the first gene in dummy data
    run_view_images(
        data_type="CosMx",
        prev_module_dir=prev_dir,
        module_dir=out_dir,
        gene_list=gene_list,
    )

    # Assertions
    assert (out_dir / "leiden_clusters.png").exists()
    assert (out_dir / "gene_expression.png").exists()
    assert (out_dir / "adata.h5ad").exists()


def test_spatial_statistics(dummy_adata, tmp_path):
    # Setup
    dummy_adata.obs["leiden"] = np.random.choice(["0", "1"], size=dummy_adata.n_obs)
    # Squidpy requires categorical data for clusters
    dummy_adata.obs["leiden"] = dummy_adata.obs["leiden"].astype("category")

    prev_dir = tmp_path / "Module4"
    prev_dir.mkdir()
    dummy_adata.write_h5ad(prev_dir / "adata.h5ad")

    out_dir = tmp_path / "Module5"

    # Action
    run_spatial_statistics(module_dir=out_dir, prev_module_dir=prev_dir)

    # Assertions
    # Did Squidpy compute the spatial connectivities?
    processed_adata = sc.read_h5ad(out_dir / "adata.h5ad")
    assert "spatial_connectivities" in processed_adata.obsp

    # Were the plots and CSVs generated?
    assert (out_dir / "centrality_scores.png").exists()
    assert (out_dir / "co_occurrence.png").exists()
    assert (out_dir / "moranI_results.csv").exists()
