"""Generates selection CSV file based on user criteria"""

import pandas as pd
import numpy as np
import warnings
import logging
from logging import getLogger
import scanpy as sc
from config import settings, get_module
from pathlib import Path

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)
logger = getLogger(__name__)


def cosmx_csv(
    module_dir,
    prev_module_dir,
    selection_name,
    cluster_col,
    selected_fovs,
    selected_celltypes,
):
    """
    Generates a CSV file with cell_ID, fov, and cluster_name for
    cells that match the selected fovs and cell types (clusters).
    For CosMx data.
    """
    adata_path = Path(prev_module_dir) / "adata.h5ad"
    logger.info(f"Loading Annotated AnnData from {adata_path}...")
    adata = sc.read_h5ad(adata_path)

    # Create selection DataFrame
    df_selection = pd.DataFrame(
        {
            "cell_ID": adata.obs["cell_ID"],
            "fov": adata.obs["fov"],
            "cluster_name": adata.obs[cluster_col].astype(str),
        }
    )

    # Apply filtering based on cell type (cluster/ cell annotation) and fov
    if selected_fovs is not None:
        df_selection = df_selection[df_selection["fov"].isin(selected_fovs)]

    if selected_celltypes is not None:
        df_selection = df_selection[
            df_selection["cluster_name"].isin(selected_celltypes)
        ]

    # save selection using a custom selection name
    output_csv_path = f"{module_dir}/{selection_name}_CosMxselection.csv"
    output_dir = Path(module_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df_selection.to_csv(output_csv_path, index=False)

    logger.info(f"Successfully created {output_csv_path}")
    logger.info(df_selection.head())
    logger.info(f"Number of cells exported: {len(df_selection)}")

    return output_csv_path


def xenium_csv(
    module_dir,
    prev_module_dir,
    selection_name,
    genes_of_interest,
    cluster_col,
    box_ids=None,
    grid_csv_path=None,
):
    """
    Generates a CSV file with cell_ID, cluster, transcript count for
    cells that match the selected clusters. For Xenium data.

    Parameters:
    ...
    box_ids : int, str, or list, optional
        The ID(s) of the box(es) from the generated spatial grid to filter by.
        Can be a single integer/string or a list of integers/strings.
    grid_csv_path : str or Path, optional
        Path to the CSV file containing Box_ID, x_min, x_max, y_min, y_max coordinates.
    """

    adata_path = Path(prev_module_dir) / "adata.h5ad"
    logger.info(f"Loading Annotated AnnData from {adata_path}...")
    adata = sc.read_h5ad(adata_path)

    # Optional Spatial filtering based on Box IDs and Grid CSV
    if box_ids is not None and grid_csv_path is not None:
        logger.info(f"Reading grid coordinates from {grid_csv_path}...")

        # Standardize box_ids to a list of integers
        if isinstance(box_ids, (int, str)):
            box_ids = [int(box_ids)]
        else:
            box_ids = [int(b) for b in box_ids]

        # Load the coordinate CSV
        df_grid = pd.read_csv(grid_csv_path)

        # Filter the dataframe for the requested Box_IDs
        valid_boxes = df_grid[df_grid["Box_ID"].isin(box_ids)]

        if valid_boxes.empty:
            logger.error(
                f"None of the provided Box_IDs {box_ids} were found in the CSV."
            )
            raise ValueError(f"Provided Box_IDs {box_ids} are invalid.")

        # Warn if some IDs weren't found
        missing_boxes = set(box_ids) - set(valid_boxes["Box_ID"])
        if missing_boxes:
            logger.warning(
                f"The following Box_IDs were not found and will be ignored: {missing_boxes}"
            )

        # Extract coordinates from AnnData
        if "spatial" in adata.obsm.keys():
            x_coords = adata.obsm["spatial"][:, 0]
            y_coords = adata.obsm["spatial"][:, 1]
        else:
            # Fallback: Sometimes Xenium coordinates are stored in obs
            x_coords = adata.obs["x_centroid"]
            y_coords = adata.obs["y_centroid"]

        # Initialize an all-False mask (meaning no cells are selected yet)
        combined_spatial_mask = np.zeros(adata.n_obs, dtype=bool)

        # Loop through each valid box and add its cells to the mask
        for _, row in valid_boxes.iterrows():
            x_min, x_max = row["x_min"], row["x_max"]
            y_min, y_max = row["y_min"], row["y_max"]

            logger.info(
                f"Adding Box {int(row['Box_ID'])} to selection (X: {x_min} to {x_max}, Y: {y_min} to {y_max})."
            )

            # Mask for current box
            current_box_mask = (
                (x_coords >= x_min)
                & (x_coords <= x_max)
                & (y_coords >= y_min)
                & (y_coords <= y_max)
            )

            # Combine using logical OR (keeps cells that are in ANY of the boxes)
            combined_spatial_mask = combined_spatial_mask | current_box_mask

        # Subset the AnnData object using the combined mask
        adata = adata[combined_spatial_mask].copy()
        logger.info(
            f"Filtered by Box IDs {list(valid_boxes['Box_ID'])}. {adata.n_obs} cells remaining."
        )

    # Only try to slice genes that actually exist in your adata to prevent errors
    valid_genes = [g for g in genes_of_interest if g in adata.var_names]
    adata_subset = adata[:, valid_genes]

    counts_matrix = adata_subset.layers["counts"]

    if hasattr(counts_matrix, "todense"):
        custom_transcript_counts = np.array(counts_matrix.sum(axis=1)).flatten()
    else:
        # Fallback if it's not a sparse matrix
        custom_transcript_counts = np.array(counts_matrix.sum(axis=1))

    df_selection = pd.DataFrame(
        {
            "Cell ID": adata.obs["cell_id"]
            if "cell_id" in adata.obs
            else adata.obs.index,
            "Cluster": "Cluster " + adata.obs[cluster_col].astype(str),
            "Transcripts": custom_transcript_counts,
            "Area (µm^2)": adata.obs["cell_area"],
        }
    )

    # turn into integer format
    df_selection["Transcripts"] = df_selection["Transcripts"].astype(int)
    # Round the area to 2 decimal places
    df_selection["Area (µm^2)"] = df_selection["Area (µm^2)"].round(2)

    # Write to CSV
    total_selection_area = df_selection["Area (µm^2)"].sum()
    output_dir = Path(module_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv_path = f"{output_dir}/{selection_name}_XeniumSelection.csv"

    with open(output_csv_path, "w") as f:
        f.write(f"#Selection name: {selection_name}\n")
        f.write(f"#Area (µm^2): {total_selection_area:.2f}\n")

    df_selection.to_csv(output_csv_path, mode="a", index=False)

    logger.info(f"Successfully created {output_csv_path}")
    logger.info(df_selection.head())
    logger.info(f"Number of cells exported: {len(df_selection)}")

    return output_csv_path


if __name__ == "__main__":
    data_type = settings["project"]["data_type"]

    _, module_5_dir = get_module(5)
    module_6_name, module_6_dir = get_module(6)
    ms_settings = settings["modules"]["MuSpan"]
    cluster_name = settings["modules"]["Annotate"]["chosen_cluster"]
    cell_types = ms_settings["cell_types"]
    transcript_list = ms_settings["transcripts"]

    selection_name = ms_settings["selection_name"]
    selected_fovs = ms_settings["selected_fovs"]
    selected_celltypes = ms_settings["selected_celltypes"]

    box_ids = ms_settings["box_ids"]
    grid_csv_path = "path"

    if data_type == "CosMx":
        logger.info("Creating Selection CSVs...")
        cosmx_csv(
            module_dir=module_6_dir,
            prev_module_dir=module_5_dir,
            selection_name=selection_name,
            cluster_col=cluster_name,
            selected_fovs=selected_fovs,
            selected_celltypes=selected_celltypes,
        )

    elif data_type == "Xenium":
        logger.info("Creating Selection CSVs...")
        # not final yet, need more filtering steps
        cell_selection_csv = xenium_csv(
            selection_name=selection_name,
            genes_of_interest=transcript_list,
            cluster_col=cluster_name,
            box_ids=box_ids,
            grid_csv_path=grid_csv_path,
        )
