# tests/test_muspan.py
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from MuSpan import run_muspan

# Skip this entire file if the user hasn't pip installed muspan yet
pytest.importorskip("muspan")


@patch("MuSpan.sd.read_zarr")  # Mocks the Zarr reader
@patch("MuSpan.io.imread")  # Mocks the Tiff reader
def test_run_muspan_cosmx(mock_imread, mock_read_zarr, dummy_adata, tmp_path):
    # 1. SETUP: Mocking the heavy spatial data
    # Fake an image matrix for Tiff reading
    import numpy as np

    fake_mask = np.zeros((100, 100), dtype=np.uint16)
    fake_mask[10:30, 10:30] = 1  # A square "cell" with ID 1
    mock_imread.return_value = fake_mask

    # Fake a Zarr object response
    mock_sdata = MagicMock()
    mock_pts = MagicMock()

    mock_pts.__getitem__.return_value = (
        mock_pts  # So that sdata.points["1_points"] works
    )
    # Fake dataframe for transcripts
    mock_pts.compute.return_value = pd.DataFrame(
        {"x_global_px": [10, 20], "y_global_px": [10, 20], "target": ["GeneA", "GeneB"]}
    )
    mock_sdata.points = {"1_points": mock_pts}
    mock_sdata.tables = {"table": dummy_adata}
    mock_read_zarr.return_value = mock_sdata

    # 1. Create the CellLabels directory
    label_dir = tmp_path / "CellLabels"
    label_dir.mkdir()

    # 2. Create a dummy file so glob doesn't return an empty list
    # The code looks for F001 because fov_id="1" is padded to 3 digits
    dummy_mask = label_dir / "sample_F001_labels.tif"
    dummy_mask.write_text("dummy content")

    # Fake previous module
    dummy_adata.obs["cell_type"] = "Cluster_0"
    prev_dir = tmp_path / "Module5"
    prev_dir.mkdir()
    dummy_adata.write_h5ad(prev_dir / "adata.h5ad")

    out_dir = tmp_path / "Module6"

    # 2. ACTION
    # We expect this to run through to the end using our faked data
    domain = run_muspan(
        dataset_type="CosMx",
        module_dir=str(out_dir),
        prev_module_dir=str(prev_dir),
        domain_name="Test_Domain",
        cluster_labels="cell_type",
        transcripts_of_interest=["GeneA"],
        zarr_path="fake_path.zarr",
        flat_files_dir=str(tmp_path),
        fov_id="1",
    )

    # 3. ASSERTIONS
    assert domain is not None
    assert domain.name == "Test_Domain"
    # Check if the visualization image saved
    assert (out_dir / "muspan_fov_1_visualization.png").exists()
