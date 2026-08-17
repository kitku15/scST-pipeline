import os
import sys
from pathlib import Path

# 1. Compatibility handling
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

# 2. Load the config dynamically
# Use an environment variable set by __main__.py, or default to CosMx
CONFIG_PATH_STR = os.getenv("RECODE_CONFIG", "config_CosMx.toml")
CONFIG_PATH = Path(CONFIG_PATH_STR)

if not CONFIG_PATH.exists():
    raise FileNotFoundError(f"Configuration file not found: {CONFIG_PATH}")

with open(CONFIG_PATH, "rb") as f:
    settings = tomllib.load(f)

# 3. Pull analysis_name from TOML and build the directory path
analysis_name = settings["project"]["analysis_name"]
analysis_dir = Path(f"analysis_{analysis_name}")

# 4. Generate the MODULES Dictionary
MODULES = {
    name: {"name": name, "dir": analysis_dir / name}
    for name in settings["pipeline"]["modules"]
}

# Create the base analysis directory if it doesn't exist
analysis_dir.mkdir(parents=True, exist_ok=True)


def get_module(index):
    """Returns the name and dir for a module starting with 'n_'"""
    key = next((k for k in MODULES if k.startswith(f"{index}_")), None)

    if key:
        return MODULES[key]["name"], MODULES[key]["dir"]
    raise ValueError(f"Module starting with {index}_ not found.")
