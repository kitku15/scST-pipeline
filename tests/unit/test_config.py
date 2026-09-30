import pytest
from pydantic import ValidationError
from Modules.config import ProjectConfig, IOConfig


def test_project_config_valid():
    """Test that a valid project configuration is accepted."""
    config_data = {
        "analysis_name": "Test_Run",
        "data_type": "CosMx",
        "batch_key": "slide_id",
        "sample_key": "sample_id",
    }
    # This should succeed without errors
    config = ProjectConfig(**config_data)
    assert config.analysis_name == "Test_Run"
    assert config.data_type == "CosMx"


def test_project_config_invalid_datatype():
    """Test that an unsupported data_type raises an error immediately."""
    config_data = {
        "analysis_name": "Test_Run",
        "data_type": "Visium",  # Invalid! Only CosMx and Xenium are allowed
    }

    # We expect Pydantic to raise a ValidationError
    with pytest.raises(ValidationError) as exc_info:
        ProjectConfig(**config_data)

    # Verify the error message complains about the pattern we set in config.py
    assert "String should match pattern '^(CosMx|Xenium)$'" in str(exc_info.value)


def test_io_config_defaults():
    """Test that the IO configuration sets the correct defaults if left empty."""
    io_config = IOConfig()

    # Check that default entry_point is "0"
    assert io_config.entry_point == "0"
    assert io_config.base_zarr_dir == "data_zarrs"
