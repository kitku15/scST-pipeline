from Modules.QualityControl import filter_and_normalize
from Modules.config import QCConfig


def test_filter_and_normalize_removes_bad_cells(mock_adata):
    """Test that cells below thresholds are properly dropped."""

    # 1. Sabotage our mock data!
    # Cell 0 has too few counts
    mock_adata.obs.loc["cell_0", "total_counts"] = 5
    # Cell 1 is too small (debris)
    mock_adata.obs.loc["cell_1", "Area"] = 2.0

    # 2. Setup the configuration
    # Note: min_cells=0 so we don't accidentally drop all genes in our tiny 10-cell mock dataset
    qc_config = QCConfig(min_counts=50, min_area=10.0, max_area=200.0, min_cells=0)
    cfg = {"has_dapi": True, "area_col": "Area", "nucleus_col": "Mean.DAPI"}

    # 3. Run your function
    filtered_adata = filter_and_normalize(mock_adata, qc_config, cfg)

    # 4. Assertions
    # We started with 10 cells, sabotaged 2, so we should have exactly 8 left.
    assert filtered_adata.n_obs == 8, "Failed to drop the correct number of cells."
    assert "cell_0" not in filtered_adata.obs_names, "Low count cell wasn't filtered!"
    assert "cell_1" not in filtered_adata.obs_names, "Low area cell wasn't filtered!"


def test_filter_and_normalize_preserves_raw_counts(mock_adata):
    """Test the data contract that .layers['counts'] is saved before normalization."""

    qc_config = QCConfig(min_cells=0)
    cfg = {"has_dapi": True, "area_col": "Area", "nucleus_col": "Mean.DAPI"}

    filtered_adata = filter_and_normalize(mock_adata, qc_config, cfg)

    # Check if the layers contract was fulfilled
    assert "counts" in filtered_adata.layers, "Raw counts layer was not preserved!"
    assert filtered_adata.raw is not None, ".raw was not populated!"
