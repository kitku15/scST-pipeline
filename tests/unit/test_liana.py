import numpy as np
import pandas as pd
import anndata as ad
from Modules.LIANA_CCC import _build_sample_aware_spatial_graph


def test_build_sample_aware_spatial_graph_prevents_ghost_edges():
    """Ensure cells in different samples but identical coordinates do not connect."""

    # 1. We need >100 cells per sample to satisfy LIANA's internal KNN optimization.
    n_dummy = 105

    # Create our 4 test cells PLUS the dummy cells
    sample_ids = (
        ["Sample_A", "Sample_A", "Sample_B", "Sample_B"]
        + ["Sample_A"] * n_dummy
        + ["Sample_B"] * n_dummy
    )

    obs = pd.DataFrame({"sample_id": sample_ids})
    obs.index = [f"cell_{i}" for i in range(len(sample_ids))]

    # 2. Coordinates!
    # Our critical test cells
    test_coords = [
        [0.0, 0.0],  # Index 0: Sample A, Cell 1
        [5.0, 5.0],  # Index 1: Sample A, Cell 2 (Close enough to connect)
        [0.0, 0.0],  # Index 2: Sample B, Cell 1 (Identical to A1, MUST NOT connect)
        [5.0, 5.0],  # Index 3: Sample B, Cell 2
    ]

    # Dummy cells: Put them far away (1000+ coordinates) so they don't interfere with our test
    dummy_coords = [[1000.0 + i, 1000.0 + i] for i in range(n_dummy * 2)]
    spatial_coords = np.array(test_coords + dummy_coords)

    # Create the AnnData object
    adata = ad.AnnData(obs=obs, obsm={"spatial": spatial_coords})

    # 3. Run your custom graph building function
    _build_sample_aware_spatial_graph(
        adata=adata,
        sample_key="sample_id",
        spatial_key="spatial",
        bandwidth=10,
        cutoff=0.01,
    )

    # 4. Extract the resulting adjacency matrix as a standard dense array
    conn = adata.obsp["spatial_connectivities"].toarray()

    # --- ASSERTIONS ---

    # Index 0 and 1 are in Sample A and are 5 units apart (bandwidth is 10). They SHOULD connect.
    assert conn[0, 1] > 0, (
        "Intra-sample connection failed! Cells in the same sample didn't connect."
    )

    # Index 0 is Sample A, Index 2 is Sample B. They share exact (0,0) coordinates.
    # They MUST NOT be connected because the graph is sample-aware.
    assert conn[0, 2] == 0, "GHOST EDGE DETECTED! Cells in different samples connected."
