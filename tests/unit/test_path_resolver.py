import pytest
from unittest.mock import patch
from Modules.path_resolver import PathResolver


@pytest.fixture
def mock_settings():
    """A mock settings dictionary for the PathResolver."""
    return {"io": {"entry_point": "0", "custom_input_adata": None}}


# We use @patch to "fake" the get_module function so we don't rely on global states
@patch("Modules.path_resolver.get_module")
def test_linear_resolution(mock_get_module, mock_settings, tmp_path):
    """Test that Module 3 asks for data from Module 2."""

    # Setup our "fake" Module 2 directory and create a fake adata.h5ad in it
    fake_mod2_dir = tmp_path / "2_DimensionReduction"
    fake_mod2_dir.mkdir()
    fake_adata = fake_mod2_dir / "adata.h5ad"
    fake_adata.touch()  # Creates an empty file

    # Tell our fake get_module to return this directory when asked for Module 2
    mock_get_module.return_value = ("2_DimensionReduction", fake_mod2_dir)

    # Initialize the resolver
    resolver = PathResolver(analysis_dir=tmp_path, settings=mock_settings)

    # Ask the resolver: "Where should Module 3 get its data?"
    resolved_path = resolver.get_adata_for_module("3")

    # Assert it correctly routed to Module 2's output
    assert resolved_path == fake_adata


def test_byod_custom_input(tmp_path):
    """Test the Bring Your Own Data (BYOD) logic."""

    # Create a fake custom dataset on disk
    custom_data_path = tmp_path / "my_custom_data.h5ad"
    custom_data_path.touch()

    settings_byod = {
        "io": {
            "entry_point": "4",  # Entering at Module 4
            "custom_input_adata": str(custom_data_path),
        }
    }

    resolver = PathResolver(analysis_dir=tmp_path, settings=settings_byod)

    # Ask the resolver: "Where should Module 4 get its data?"
    resolved_path = resolver.get_adata_for_module("4")

    # Assert it grabbed the custom data instead of looking for Module 3
    assert resolved_path == custom_data_path
