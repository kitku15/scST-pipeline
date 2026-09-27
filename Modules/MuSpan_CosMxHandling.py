"""MuSpan CosMx handling module."""

import gc
import logging
import warnings
from pathlib import Path
from typing import Tuple, List, Dict, Set, Optional, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.validation import make_valid

try:
    import muspan as ms
except ModuleNotFoundError as err:
    raise ModuleNotFoundError(
        "Could not load MuSpAn. Install it before running this script."
    ) from err

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# Hardware constants
COSMX_PIXEL_SIZE_UM = 0.12


class MuSpanGeometryError(Exception):
    """Custom exception for errors parsing spatial geometries."""

    pass


def _load_cell_selection(
    csv_path: str,
) -> Tuple[Dict[str, str], Dict[str, Set[str]], List[str]]:
    """Loads and standardizes cell selection criteria."""
    logger.info("Loading cell selection from CSV...")
    df_sel = pd.read_csv(
        csv_path, dtype={"cell_ID": str, "fov": str, "cluster_name": str}
    )

    cell_id_to_type = {
        (str(row.fov), str(row.cell_ID)): str(row.cluster_name)
        for row in df_sel.itertuples()
    }

    selected_cells_by_fov = {}
    fov_list = df_sel["fov"].unique().tolist()

    for f in fov_list:
        cells = set(df_sel[df_sel["fov"] == f]["cell_ID"].tolist())
        selected_cells_by_fov[f] = cells
        logger.debug(f"FOV {f}: Identified {len(cells)} valid local cells.")

    if not fov_list:
        raise ValueError("No FOV provided via fov_id or CSV.")

    return cell_id_to_type, selected_cells_by_fov, fov_list


def _process_proseg_zarr(
    proseg_zarr_path: str,
) -> Tuple[Any, Dict[Tuple[str, str], int]]:
    """Extracts boundaries and metadata from Proseg Zarr."""
    try:
        import spatialdata

        logger.info(f"Loading Proseg data from {proseg_zarr_path}")
        proseg_sdata = spatialdata.read_zarr(proseg_zarr_path)

        proseg_gdf = proseg_sdata.shapes.get("cell_boundaries")
        if hasattr(proseg_gdf, "compute"):
            proseg_gdf = proseg_gdf.compute()

        proseg_id_map = {}
        if "table" in proseg_sdata.tables:
            proseg_obs = proseg_sdata.tables["table"].obs
            if "fov" in proseg_obs.columns and "cell_ID" in proseg_obs.columns:
                for idx, row in proseg_obs.iterrows():
                    proseg_id_map[
                        (str(row["fov"]), str(row["cell_ID"]).replace(".0", ""))
                    ] = idx
            elif "original_cell_id" in proseg_obs.columns:
                for idx, row in proseg_obs.iterrows():
                    parts = str(row["original_cell_id"]).split("_")
                    if len(parts) >= 2:
                        proseg_id_map[(parts[-2], parts[-1].replace(".0", ""))] = idx

        return proseg_gdf, proseg_id_map
    except Exception as e:
        logger.error(f"Failed to load Proseg Zarr: {e}")
        return None, {}


def _get_polygon_coords(poly: Polygon) -> Optional[np.ndarray]:
    """Safely extracts exterior coordinates from a Shapely polygon."""
    if not poly.is_valid:
        poly = make_valid(poly)
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda a: a.area)
        elif poly.geom_type == "GeometryCollection":
            polys = [geom for geom in poly.geoms if geom.geom_type == "Polygon"]
            if not polys:
                return None
            poly = max(polys, key=lambda a: a.area)

    if poly.area <= 0:
        return None

    coords = np.array(poly.exterior.coords)
    # Ensure Counter-Clockwise orientation (Shoelace formula)
    x, y = coords[:, 0], coords[:, 1]
    if 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) < 0:
        coords = coords[::-1]

    return coords


def CosMx_to_domain(
    sdata: Any,
    domain: Any,
    transcripts_of_interest: List[str],
    cluster_labels: str,
    cell_selection_csv: str,
    flat_files_dir: str,
    out_dir: Path,
    proseg_zarr_path: Optional[str] = None,
) -> Tuple[Any, Any, Any, Any]:
    """Main function to populate MuSpAn domain with CosMx features."""

    cell_id_to_type, selected_cells_by_fov, fov_list = _load_cell_selection(
        cell_selection_csv
    )

    all_trans_coords, all_trans_targets = [], []
    all_cell_coords, all_cell_ids, all_cell_types = [], [], []
    all_geometries, all_boundary_types = [], []

    proseg_gdf, proseg_id_map = None, {}
    if proseg_zarr_path:
        proseg_gdf, proseg_id_map = _process_proseg_zarr(proseg_zarr_path)

    # Fallback to CSV
    df_poly, poly_x_col, poly_y_col = None, None, None
    if proseg_gdf is None or not proseg_id_map:
        poly_files = list(Path(flat_files_dir).glob("*polygons*.csv"))
        if poly_files:
            logger.info(f"Loading flat polygons from {poly_files[0].name}...")
            df_poly = pd.read_csv(poly_files[0])
            df_poly.rename(columns={"cellID": "cell_ID"}, inplace=True, errors="ignore")

            if "x_global_px" in df_poly.columns:
                poly_x_col, poly_y_col = "x_global_px", "y_global_px"
            elif "x_local_px" in df_poly.columns:
                poly_x_col, poly_y_col = "x_local_px", "y_local_px"

    for current_fov in fov_list:
        logger.info(f"--- Processing FOV {current_fov} ---")
        local_selected_cells = selected_cells_by_fov.get(current_fov)

        points_key = f"{current_fov}_points"
        if points_key not in sdata.points:
            continue

        # 1. Memory-Safe Transcript Parsing
        df_all_pts = sdata.points[points_key].compute()
        if local_selected_cells is not None and "cell_ID" in df_all_pts.columns:
            pts_cell_ids = (
                df_all_pts["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
            )
            df_pts_plot = df_all_pts[pts_cell_ids.isin(local_selected_cells)]
        else:
            df_pts_plot = df_all_pts

        if transcripts_of_interest:
            df_pts_plot = df_pts_plot[
                df_pts_plot["target"].isin(transcripts_of_interest)
            ]

        if not df_pts_plot.empty:
            all_trans_coords.append(
                df_pts_plot[["x_global_px", "y_global_px"]].values * COSMX_PIXEL_SIZE_UM
            )
            all_trans_targets.append(df_pts_plot["target"].astype(str).values)

        # Clear heavy transcript dataframe immediately
        del df_all_pts, df_pts_plot
        gc.collect()

        # 2. Extract Valid Cells and True Centers
        sdata_table = sdata.tables["table"]
        adata_fov = sdata_table[sdata_table.obs["fov"].astype(str) == str(current_fov)]
        if adata_fov.n_obs == 0:
            continue

        valid_local_ids = set(
            adata_fov.obs["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
        )

        # THE FIX: Use enumerate to get integer position index for the .obsm array
        global_coords = adata_fov.obsm["global"]
        true_centers = {}
        for pos_idx, c_id in enumerate(adata_fov.obs["cell_ID"]):
            clean_id = str(c_id).replace(".0", "")
            true_centers[clean_id] = (
                global_coords[pos_idx, 0] * COSMX_PIXEL_SIZE_UM,
                global_coords[pos_idx, 1] * COSMX_PIXEL_SIZE_UM,
            )

        # 3. Process Boundaries
        if proseg_gdf is not None and proseg_id_map:
            for cid_str in valid_local_ids:
                if local_selected_cells and cid_str not in local_selected_cells:
                    continue

                proseg_idx = proseg_id_map.get((str(current_fov), cid_str))
                if proseg_idx is None:
                    continue

                try:
                    geom = proseg_gdf.loc[int(proseg_idx)].geometry
                except KeyError:
                    geom = proseg_gdf.loc[str(proseg_idx)].geometry

                poly_coords = _get_polygon_coords(geom)
                if poly_coords is None:
                    continue

                poly_coords *= COSMX_PIXEL_SIZE_UM
                gcx, gcy = np.mean(poly_coords[:, 0]), np.mean(poly_coords[:, 1])
                c_type = cell_id_to_type.get((str(current_fov), cid_str), "Unknown")

                all_geometries.append(poly_coords)
                all_cell_coords.append([gcx, gcy])
                all_cell_ids.append(f"{current_fov}_{cid_str}")
                all_cell_types.append(c_type)
                all_boundary_types.append(c_type)

        elif df_poly is not None and poly_x_col:
            fov_poly = df_poly[df_poly["fov"].astype(str) == str(current_fov)]
            for cid, group in fov_poly.groupby("cell_ID"):
                cid_str = str(cid).replace(".0", "")
                if cid_str not in valid_local_ids or (
                    local_selected_cells and cid_str not in local_selected_cells
                ):
                    continue
                if cid_str not in true_centers:
                    continue

                pt_x, pt_y = true_centers[cid_str]
                raw_x, raw_y = group[poly_x_col].values, group[poly_y_col].values

                # Shift and scale
                shifted_x = pt_x + ((raw_x - np.mean(raw_x)) * COSMX_PIXEL_SIZE_UM)
                shifted_y = pt_y + ((raw_y - np.mean(raw_y)) * COSMX_PIXEL_SIZE_UM)

                poly_coords = _get_polygon_coords(
                    Polygon(np.column_stack((shifted_x, shifted_y)))
                )
                if poly_coords is None:
                    continue

                c_type = cell_id_to_type.get((str(current_fov), cid_str), "Unknown")

                all_geometries.append(poly_coords)
                all_cell_coords.append(
                    [np.mean(poly_coords[:, 0]), np.mean(poly_coords[:, 1])]
                )
                all_cell_ids.append(f"{current_fov}_{cid_str}")
                all_cell_types.append(c_type)
                all_boundary_types.append(c_type)

        else:
            # Fallback to centroids
            for cid_str in valid_local_ids:
                if local_selected_cells is None or cid_str in local_selected_cells:
                    if cid_str in true_centers:
                        all_cell_coords.append(list(true_centers[cid_str]))
                        all_cell_ids.append(f"{current_fov}_{cid_str}")
                        all_cell_types.append(
                            cell_id_to_type.get(
                                (str(current_fov), str(cid_str)), "Unknown"
                            )
                        )

    logger.info("Adding collected objects to MuSpAn domain...")
    qTrans, qCells, qBoundaries = None, None, None

    if all_trans_coords:
        domain.add_points(np.vstack(all_trans_coords), collection_name="Transcripts")
        qTrans = ms.query.query(domain, ("collection",), "is", "Transcripts")
        domain.add_labels(
            "Transcript ID", np.concatenate(all_trans_targets), add_labels_to=qTrans
        )

    if all_cell_coords:
        domain.add_points(np.vstack(all_cell_coords), collection_name="Cell centroids")
        qCells = ms.query.query(domain, ("collection",), "is", "Cell centroids")
        domain.add_labels("Cell ID", all_cell_ids, add_labels_to=qCells)
        domain.add_labels(cluster_labels, all_cell_types, add_labels_to=qCells)

    if all_geometries:
        domain.add_shapes(shapes=all_geometries, collection_name="Cell boundaries")
        qBoundaries = ms.query.query(domain, ("collection",), "is", "Cell boundaries")
        domain.add_labels(cluster_labels, all_boundary_types, add_labels_to=qBoundaries)

    return domain, qTrans, qCells, qBoundaries


def cosmx_initial_plotting(
    domain: Any,
    qTrans: Any,
    qCells: Any,
    qBoundaries: Any,
    cluster_labels: str,
    out_dir: Path,
    color_dict: Optional[dict] = None,
) -> None:
    """Renders the initial visualization figures."""
    if color_dict is not None:
        try:
            domain.update_colors(
                color_dict, colors_to_update="labels", label_name=cluster_labels
            )
            logger.info("Successfully synced Squidpy colors to MuSpAn CosMx Domain.")
        except Exception as e:
            logger.warning(f"Failed to apply custom colors to MuSpAn domain: {e}")

    fig, ax = plt.subplots(figsize=(20, 10), nrows=1, ncols=2)

    if qTrans is not None:
        ms.visualise.visualise(
            domain,
            color_by=("label", "Transcript ID"),
            ax=ax[0],
            objects_to_plot=qTrans,
            marker_size=0.5,
        )
        ax[0].set_title("Filtered Transcripts")
    else:
        ax[0].set_title("No Transcripts Found")

    if qBoundaries is not None:
        ms.visualise.visualise(
            domain,
            color_by=("label", cluster_labels),
            ax=ax[1],
            objects_to_plot=qBoundaries,
        )
        ax[1].set_title(f"Cell Boundaries by {cluster_labels}")
    elif qCells is not None:
        ms.visualise.visualise(
            domain,
            color_by=("label", cluster_labels),
            ax=ax[1],
            objects_to_plot=qCells,
            marker_size=20,
        )
        ax[1].set_title(f"Cell Centroids by {cluster_labels} (No Boundaries)")

    plt.tight_layout()
    domain_dir = out_dir / "domain"
    domain_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(domain_dir / "muspan_cosmx_visualization.png", dpi=300)
    plt.close()
