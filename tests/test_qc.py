# tests/test_qc.py
import scanpy as sc
from QualityControl import run_qc
from unittest.mock import patch, MagicMock


@patch("QualityControl.sd.read_zarr")
def test_quality_control(mock_read_zarr, dummy_adata, tmp_path):
    # 1. SETUP
    # Mock the Zarr reader so it returns our tiny dummy_adata instead of a 10GB file
    mock_sdata = MagicMock()
    mock_sdata.tables = {"table": dummy_adata}
    mock_read_zarr.return_value = mock_sdata

    out_dir = tmp_path / "Module1"

    # 2. ACTION: Run QC with strict cutoffs so we know cells get filtered
    run_qc(
        data_type="CosMx",
        module_dir=out_dir,
        zarr_path="fake.zarr",
        min_counts=5,
        min_cells=1,
        min_dapi=5,
    )

    # 3. ASSERTIONS
    processed_adata = sc.read_h5ad(out_dir / "adata.h5ad")

    # Was the raw data saved?
    assert processed_adata.raw is not None

    # Were the plots created?
    assert (out_dir / "cell_summary_histograms.png").exists()

    # Check if highly variable genes were calculated
    assert "highly_variable" in processed_adata.var.columns
