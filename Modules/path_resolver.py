"""Intelligent Path Resolution for BYOD and DAG execution."""

from pathlib import Path
from logging import getLogger
from config import get_module

logger = getLogger(__name__)


class PathResolver:
    def __init__(self, analysis_dir: Path, settings: dict):
        self.analysis_dir = analysis_dir

        # Safely get io settings (using dict or object access)
        io_settings = settings.get("io", {})
        if hasattr(io_settings, "entry_point"):
            self.entry_point = str(io_settings.entry_point)
            self.custom_adata = (
                Path(io_settings.custom_input_adata)
                if io_settings.custom_input_adata
                else None
            )
        else:
            self.entry_point = str(io_settings.get("entry_point", "0"))
            custom_path = io_settings.get("custom_input_adata")
            self.custom_adata = Path(custom_path) if custom_path else None

    def get_adata_for_module(self, target_module_prefix: str) -> Path:
        """Resolves where the input AnnData should come from based on the DAG."""
        target = str(target_module_prefix)

        # 1. If entering AT this specific module, use custom file
        if self.entry_point == target:
            if not self.custom_adata or not self.custom_adata.exists():
                raise FileNotFoundError(
                    f"Entry point is set to '{target}', but custom_input_adata "
                    f"'{self.custom_adata}' is missing or invalid."
                )
            logger.info(
                f"BYOD Active: Injecting {self.custom_adata} into Module {target}"
            )
            return self.custom_adata

        # 2. Parallel Fan-Out: Modules 4, 5, 6, 7, 8, 8b, 8c, 9, 10 ALL read from Module 3 (Annotate)
        if target in ["4", "5", "6", "7", "8", "8b", "8c", "9", "10"]:
            _, mod3_dir = get_module(3)
            return self._find_adata(mod3_dir)

        # 3. Linear Upstream Steps
        if target == "3":
            _, mod_dir = get_module(2)
        elif target == "2":
            try:
                _, mod_dir = get_module("1b")
                # If 1b folder doesn't exist (no merge happened), fallback to 1
                if not mod_dir.exists():
                    _, mod_dir = get_module(1)
            except ValueError:
                _, mod_dir = get_module(1)
        else:
            raise ValueError(
                f"Module '{target}' input logic undefined in PathResolver."
            )

        return self._find_adata(mod_dir)

    def _find_adata(self, directory: Path) -> Path:
        file_path = directory / "adata.h5ad"
        if file_path.exists():
            return file_path

        # Check subdirectories (for single-slide QC folders)
        sub_files = list(directory.glob("*/adata.h5ad"))
        if len(sub_files) == 1:
            return sub_files[0]
        elif len(sub_files) > 1:
            raise FileNotFoundError(
                f"Multiple .h5ad files in {directory}. Did you forget to run Module 1b (Merge)?"
            )

        raise FileNotFoundError(
            f"Expected input adata in {directory} but found none. Did the previous step finish?"
        )
