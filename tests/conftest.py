# tests/conftest.py

import matplotlib

matplotlib.use("Agg")  # Forces matplotlib to run without opening GUI windows
import os
import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# ==========================================
# 1. PRE-TEST SETUP (Runs before imports)
# ==========================================
# We must create a dummy config file and set the env var immediately,
# otherwise config.py will crash with FileNotFoundError when pytest
# tries to import the scripts.

DUMMY_CONFIG_PATH = Path("pytest_dummy_config.toml")
if not DUMMY_CONFIG_PATH.exists():
    DUMMY_CONFIG_PATH.write_text("""
[project]
analysis_name = "test_run"
data_type = "CosMx"

[pipeline]
modules = ["1_QualityControl", "2_DimensionReduction"]

[io]
dataset_dir = "fake/path"
dataset_id = "fake_id"
zarr_dir = "fake/path.zarr"
    """)

# Force config.py to use this file
os.environ["RECODE_CONFIG"] = str(DUMMY_CONFIG_PATH)


# ==========================================
# 2. FIXTURES (Dummy data for tests)
# ==========================================
@pytest.fixture
def dummy_adata():
    """Generates a tiny AnnData object with fake spatial data."""
    n_cells = 100
    n_genes = 200

    # Fake counts matrix
    X = np.random.poisson(1.5, size=(n_cells, n_genes)).astype(np.float32)

    # Fake observation (cell) data
    obs = pd.DataFrame(
        {
            "cell_ID": [f"cell_{i}" for i in range(n_cells)],
            "Area": np.random.uniform(50, 150, n_cells),
            "Mean.DAPI": np.random.uniform(10, 100, n_cells),
            "total_counts": np.sum(X, axis=1),
            "n_genes_by_counts": np.sum(X > 0, axis=1),
            "fov": np.random.choice(["1", "2"], n_cells),
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    # Fake spatial coordinates
    spatial_coords = np.random.uniform(0, 1000, size=(n_cells, 2))

    # Fake var (gene) data (include a negative probe for QC tests)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    var.index.values[0] = "Negative_probe_1"

    adata = sc.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = spatial_coords
    adata.obsm["global"] = spatial_coords

    return adata


# Cleanup after all tests finish
def pytest_sessionfinish(session, exitstatus):
    """Remove the dummy config file when tests are done."""
    if DUMMY_CONFIG_PATH.exists():
        DUMMY_CONFIG_PATH.unlink()
