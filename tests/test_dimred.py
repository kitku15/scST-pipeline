import scanpy as sc
from DimensionReduction import run_dimension_reduction


def test_dimension_reduction(dummy_adata, tmp_path):
    # Setup: Pretend QC already ran
    sc.pp.normalize_total(dummy_adata)
    sc.pp.log1p(dummy_adata)
    dummy_adata.write(tmp_path / "adata.h5ad")  # Save as previous step output

    out_dir = tmp_path / "Module2"

    # Action
    run_dimension_reduction(
        data_type="CosMx",
        prev_module_dir=tmp_path,
        module_dir=out_dir,
        module_name="2_DimRed",
        n_comps=15,  # Keep it small for the test!
        n_neighbors=5,
        resolution=0.5,
        cluster_name="leiden",
    )

    # Assertions
    processed_adata = sc.read_h5ad(out_dir / "adata.h5ad")
    assert "X_pca" in processed_adata.obsm  # PCA ran
    assert "X_umap" in processed_adata.obsm  # UMAP ran
    assert "leiden" in processed_adata.obs  # Clustering ran
    assert (out_dir / "umap_2_DimRed.png").exists()
