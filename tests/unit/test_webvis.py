import json
import pandas as pd
from Modules.WebVisPrep import export_de_analysis_for_web


def test_export_de_analysis_for_web_creates_valid_json(tmp_path):
    """Test that the WebVis DE exporter creates JSONs with exactly the keys the frontend expects."""

    # 1. Setup a fake Module 3 directory
    mod3_dir = tmp_path / "3_Annotate"
    mod3_dir.mkdir()

    # Create a fake top DE genes file
    annotation_col = "leiden_n10"
    (mod3_dir / f"top_DEgenes_{annotation_col}.csv").touch()

    # Create a fake DEgenes subdirectory
    cluster_dir = mod3_dir / "DEgenes"
    cluster_dir.mkdir()

    # Create a fake cluster data CSV (mimicking Scanpy's output)
    fake_scanpy_df = pd.DataFrame(
        {
            "names": ["GeneA", "GeneB"],
            "logfoldchanges": [1.5, -0.5],
            "pvals_adj": [0.01, 0.04],
        }
    )
    fake_csv_path = cluster_dir / "cluster_0_data.csv"
    fake_scanpy_df.to_csv(fake_csv_path, index=False)

    # 2. Run the WebVis export function
    out_dir = tmp_path / "10_WebVis"
    export_de_analysis_for_web(mod3_dir, out_dir)

    # 3. Assertions
    # Did it create the output directory?
    target_de_dir = out_dir / "de_analysis"
    assert target_de_dir.exists(), "Output directory was not created!"

    # Did it create the JSON file for the cluster?
    json_path = target_de_dir / f"{annotation_col}_cluster_0.json"
    assert json_path.exists(), "Cluster JSON was not generated!"

    # Open the JSON and verify the "Data Contract" for the Web Frontend
    with open(json_path, "r") as f:
        data = json.load(f)

    assert "names" in data, "Missing 'names' key required by Spatial-VisKit!"
    assert "logfc" in data, "Missing 'logfc' key required by Spatial-VisKit!"
    assert "pvals" in data, "Missing 'pvals' key required by Spatial-VisKit!"

    # Check that the data was actually populated
    assert data["names"] == ["GeneA", "GeneB"], "Gene names were not parsed correctly!"
    assert data["logfc"] == [1.5, -0.5], "LogFC values were not parsed correctly!"
