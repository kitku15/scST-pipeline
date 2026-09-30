import anndata as ad
from Modules.MergeData import run_merge


def test_run_merge_combines_slides_and_prefixes_barcodes(mock_adata, tmp_path):
    """Test that run_merge correctly concatenates datasets and makes IDs unique."""

    # 1. Setup temporary slide directories and write mock data to them
    slide1_path = tmp_path / "Slide_1_adata.h5ad"
    slide2_path = tmp_path / "Slide_2_adata.h5ad"

    # Write the identical mock data to both paths
    mock_adata.write_h5ad(slide1_path)
    mock_adata.write_h5ad(slide2_path)

    input_files = [str(slide1_path), str(slide2_path)]
    slide_names = ["Slide_1", "Slide_2"]

    # 2. Run your merge function
    output_dir = tmp_path / "merged_output"
    out_path = run_merge(input_files, slide_names, output_dir)

    # 3. Read the result and assert
    merged_adata = ad.read_h5ad(out_path)

    # Since we merged two 10-cell datasets, we should have 20 cells
    assert merged_adata.n_obs == 20, "Merged dataset has the wrong number of cells!"

    # Check that the batch key ('slide_id') was created
    assert "slide_id" in merged_adata.obs.columns, "slide_id column missing!"

    # Check that sample_ids and obs_names got prefixed so they are globally unique
    # e.g., 'cell_0' from Slide_1 becomes 'Slide_1_cell_0'
    assert "cell_0_Slide_1" in merged_adata.obs_names, (
        "Barcodes were not suffixed correctly!"
    )

    sample_ids = merged_adata.obs["sample_id"].unique()
    assert "Slide_1_Sample_A" in sample_ids, "Sample IDs were not prefixed!"
