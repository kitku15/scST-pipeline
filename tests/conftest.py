import os
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

dummy_config_path = Path(__file__).parent / "test_data" / "test_config.toml"
os.environ["scSpatial-Kit"] = str(dummy_config_path)


@pytest.fixture
def mock_adata():
    """Generates a tiny, 10-cell x 20-gene AnnData object for testing."""
    n_cells, n_genes = 10, 20

    # Sparse expression matrix (raw counts)
    X = sparse.csr_matrix(np.random.poisson(2, size=(n_cells, n_genes)))

    # Metadata (Make all cells "good" by default)
    obs = pd.DataFrame(
        {
            "sample_id": ["Sample_A"] * n_cells,
            "total_counts": [500] * n_cells,
            "n_genes_by_counts": [15] * n_cells,
            "Area": [100.0] * n_cells,
            "Mean.DAPI": [50.0] * n_cells,
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    # Spatial coordinates
    spatial_coords = np.random.rand(n_cells, 2) * 1000

    adata = ad.AnnData(X=X, obs=obs, var=var, obsm={"global": spatial_coords})

    return adata
