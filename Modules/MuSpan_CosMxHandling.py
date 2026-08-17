"""MuSpan CosMx handling module."""

import logging
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.validation import make_valid

try:
    import muspan as ms
except ModuleNotFoundError as err:
    raise ModuleNotFoundError(
        "Could not load necessary MuSpAn package. "
        "Please install it before running this script."
    ) from err

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)
logger = getLogger(__name__)


def CosMx_to_domain(
    sdata,
    domain,
    transcripts_of_interest,
    cluster_labels,
    cell_selection_csv,
    flat_files_dir,
    out_dir,
    proseg_zarr_path=None,
):
    PIXEL_SIZE = 0.12

    logger.info("Loading cell selection from CSV...")
    df_sel = pd.read_csv(cell_selection_csv)

    # normalize types
    df_sel["cell_ID"] = df_sel["cell_ID"].astype(str)
    df_sel["fov"] = df_sel["fov"].astype(str)
    df_sel["cluster_name"] = df_sel["cluster_name"].astype(str)

    # mapping: (fov, local_cell_id) -> cluster
    cell_id_to_type = {
        (fov, cell): cluster
        for fov, cell, cluster in zip(
            df_sel["fov"], df_sel["cell_ID"], df_sel["cluster_name"]
        )
    }

    selected_cells_by_fov = {}
    fov_list = df_sel["fov"].astype(str).unique().tolist()

    for f in fov_list:
        cells_in_fov = df_sel[df_sel["fov"].astype(str) == str(f)]["cell_ID"].astype(
            str
        )
        cells_in_fov = set(cells_in_fov)
        selected_cells_by_fov[f] = cells_in_fov
        logger.info(f"FOV {f}: Identified {len(cells_in_fov)} valid local cells.")

    if not fov_list:
        raise ValueError("No FOV provided via fov_id or CSV.")

    all_trans_coords, all_trans_targets = [], []
    all_cell_coords, all_cell_ids, all_cell_types = [], [], []
    all_geometries, all_boundary_types = [], []

    # -------------------------------------------------------------
    # 1. OPTIONAL: LOAD PROSEG GEOMETRIES & CREATE MAPPING
    # -------------------------------------------------------------
    proseg_gdf = None
    proseg_id_map = {}

    if proseg_zarr_path:
        try:
            import spatialdata

            logger.info(f"Loading Proseg data from {proseg_zarr_path}")
            proseg_sdata = spatialdata.read_zarr(proseg_zarr_path)

            if "cell_boundaries" in proseg_sdata.shapes:
                proseg_gdf = proseg_sdata.shapes["cell_boundaries"]
                if hasattr(proseg_gdf, "compute"):
                    proseg_gdf = proseg_gdf.compute()
                logger.info(f"Loaded {len(proseg_gdf)} Proseg cell boundaries.")

            if "table" in proseg_sdata.tables:
                proseg_obs = proseg_sdata.tables["table"].obs
                if "fov" in proseg_obs.columns and "cell_ID" in proseg_obs.columns:
                    for idx, row in proseg_obs.iterrows():
                        f = str(row["fov"])
                        c = str(row["cell_ID"]).replace(".0", "")
                        proseg_id_map[(f, c)] = idx
                    logger.info(
                        "Successfully built (FOV, cell_ID) -> Proseg Global Index mapping."
                    )
                elif "original_cell_id" in proseg_obs.columns:
                    # Parses Proseg's combined string format (e.g., 'c_1_28_2457' -> fov='28', cell_ID='2457')
                    for idx, row in proseg_obs.iterrows():
                        orig_id = str(row["original_cell_id"])
                        parts = orig_id.split("_")
                        if len(parts) >= 2:
                            f = parts[-2]
                            c = parts[-1].replace(".0", "")
                            proseg_id_map[(f, c)] = idx
                    logger.info(
                        "Successfully built (FOV, cell_ID) -> Proseg Global Index mapping using 'original_cell_id'."
                    )
                else:
                    logger.warning(
                        "Proseg table.obs is missing 'fov', 'cell_ID', or 'original_cell_id' columns."
                    )

        except Exception as e:
            logger.error(f"Failed to load Proseg Zarr: {e}")

    # -------------------------------------------------------------
    # 2. FALLBACK: LOAD CSV POLYGONS (If Proseg isn't used)
    # -------------------------------------------------------------
    df_poly = None
    if proseg_gdf is None or not proseg_id_map:
        flat_dir_path = Path(flat_files_dir)
        poly_files = list(flat_dir_path.glob("*polygons*.csv"))

        if poly_files:
            logger.info(
                f"Found boundary file: {poly_files[0].name}. Loading polygons into memory..."
            )
            df_poly = pd.read_csv(poly_files[0])

            if "cellID" in df_poly.columns:
                df_poly.rename(columns={"cellID": "cell_ID"}, inplace=True)

            if (
                "x_global_px" not in df_poly.columns
                or "y_global_px" not in df_poly.columns
            ):
                logger.warning(
                    f"Could not find global coordinates in {poly_files[0].name}. Falling back to centroids."
                )
                df_poly = None
        else:
            logger.warning(
                "No polygons.csv file found in dataset directory. Falling back to point centroids."
            )

    for current_fov in fov_list:
        logger.info(f"--- Processing FOV {current_fov} ---")
        local_selected_cells = selected_cells_by_fov.get(current_fov, None)

        points_key = f"{current_fov}_points"
        if points_key not in sdata.points:
            continue

        df_pts_dask = sdata.points[points_key]
        df_all_pts = df_pts_dask.compute()

        # 1. Process Transcripts
        if local_selected_cells is not None and "cell_ID" in df_all_pts.columns:
            df_pts_plot = df_all_pts[
                df_all_pts["cell_ID"].astype(str).isin(local_selected_cells)
            ]
        else:
            df_pts_plot = df_all_pts

        if transcripts_of_interest is not None:
            df_pts_plot = df_pts_plot[
                df_pts_plot["target"].isin(transcripts_of_interest)
            ]

        if not df_pts_plot.empty:
            all_trans_coords.append(
                df_pts_plot[["x_global_px", "y_global_px"]].values * PIXEL_SIZE
            )
            all_trans_targets.append(df_pts_plot["target"].astype(str).values)

        # 2. Get valid cells from AnnData
        sdata_table = sdata.tables["table"]
        mask = sdata_table.obs["fov"].astype(str) == str(current_fov)
        adata_fov = sdata_table[mask]

        if adata_fov.n_obs == 0:
            continue

        valid_local_ids = set(
            adata_fov.obs["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
        )

        # 3. Process Boundaries / Centroids
        # --- A. PROSEG GEOMETRIES LOGIC ---
        if proseg_gdf is not None and proseg_id_map:
            for cid_str in valid_local_ids:
                if (
                    local_selected_cells is not None
                    and cid_str not in local_selected_cells
                ):
                    continue

                # Fetch global integer index for this specific FOV/Cell combination
                proseg_idx = proseg_id_map.get((str(current_fov), cid_str))
                if proseg_idx is None:
                    continue

                # Extract shapely geometry
                try:
                    geom = proseg_gdf.loc[int(proseg_idx)].geometry
                except KeyError:
                    geom = proseg_gdf.loc[str(proseg_idx)].geometry

                if geom.geom_type == "MultiPolygon":
                    geom = max(geom.geoms, key=lambda a: a.area)
                if geom.geom_type != "Polygon":
                    continue

                poly_coords = np.array(geom.exterior.coords) * PIXEL_SIZE

                # Fix orientation
                x = poly_coords[:, 0]
                y = poly_coords[:, 1]
                signed_area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
                if signed_area < 0:
                    poly_coords = poly_coords[::-1]

                gcx, gcy = np.mean(poly_coords[:, 0]), np.mean(poly_coords[:, 1])

                unique_cell_id = f"{current_fov}_{cid_str}"
                c_type = cell_id_to_type.get((str(current_fov), cid_str), "Unknown")

                all_geometries.append(poly_coords)
                all_cell_coords.append([gcx, gcy])
                all_cell_ids.append(unique_cell_id)
                all_cell_types.append(c_type)
                all_boundary_types.append(c_type)

        # --- B. ORIGINAL CSV LOGIC ---
        elif df_poly is not None:
            fov_poly = df_poly[df_poly["fov"].astype(str) == str(current_fov)]

            if not fov_poly.empty:
                for cid, group in fov_poly.groupby("cell_ID"):
                    cid_str = str(cid)

                    if cid_str not in valid_local_ids:
                        continue
                    if (
                        local_selected_cells is not None
                        and cid_str not in local_selected_cells
                    ):
                        continue

                    # Get polygon coordinates
                    poly_coords = (
                        group[["x_global_px", "y_global_px"]].values * PIXEL_SIZE
                    )

                    try:
                        poly = Polygon(poly_coords)
                        if not poly.is_valid:
                            poly = make_valid(poly)

                            if poly.geom_type == "MultiPolygon":
                                poly = max(poly.geoms, key=lambda a: a.area)
                            elif poly.geom_type == "GeometryCollection":
                                polys = [
                                    geom
                                    for geom in poly.geoms
                                    if geom.geom_type == "Polygon"
                                ]
                                if not polys:
                                    continue
                                poly = max(polys, key=lambda a: a.area)

                        if poly.area <= 0:
                            continue

                        poly_coords = np.array(poly.exterior.coords)

                        x = poly_coords[:, 0]
                        y = poly_coords[:, 1]
                        signed_area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])

                        if signed_area < 0:
                            poly_coords = poly_coords[::-1]

                    except Exception as e:
                        logger.debug(
                            f"Skipping corrupted polygon for cell {cid_str}: {e}"
                        )
                        continue

                    gcx, gcy = np.mean(poly_coords[:, 0]), np.mean(poly_coords[:, 1])

                    unique_cell_id = f"{current_fov}_{cid_str}"
                    c_type = cell_id_to_type.get((str(current_fov), cid_str), "Unknown")

                    all_geometries.append(poly_coords)
                    all_cell_coords.append([gcx, gcy])
                    all_cell_ids.append(unique_cell_id)
                    all_cell_types.append(c_type)
                    all_boundary_types.append(c_type)
        else:
            # FALLBACK: If polygons missing, use pure transcript point centroids
            logger.warning(
                f"FOV {current_fov}: Using point centroids as boundaries are missing."
            )

            valid_pts = df_all_pts[df_all_pts["cell_ID"].astype(str) != "0"].copy()
            valid_pts["cell_ID_str"] = (
                valid_pts["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
            )
            global_centers = (
                valid_pts.groupby("cell_ID_str")[["x_global_px", "y_global_px"]].mean()
                * PIXEL_SIZE
            )

            for cid_str in valid_local_ids:
                if cid_str in global_centers.index:
                    if local_selected_cells is None or cid_str in local_selected_cells:
                        unique_cell_id = f"{current_fov}_{cid_str}"
                        c_type = cell_id_to_type.get(
                            (str(current_fov), str(cid_str)), "Unknown"
                        )

                        gcx = global_centers.loc[cid_str, "x_global_px"]
                        gcy = global_centers.loc[cid_str, "y_global_px"]

                        all_cell_coords.append([gcx, gcy])
                        all_cell_ids.append(unique_cell_id)
                        all_cell_types.append(c_type)

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
    domain, qTrans, qCells, qBoundaries, cluster_labels, out_dir, color_dict=None
):
    if color_dict is not None:
        try:
            domain.update_colors(
                color_dict, colors_to_update="labels", label_name=cluster_labels
            )
            logger.info("Successfully synced Squidpy colors to MuSpAn CosMx Domain.")
        except Exception as e:
            logger.warning(f"Failed to apply custom colors to MuSpAn domain: {e}")

    # Plotting (CosMx)
    fig, ax = plt.subplots(figsize=(20, 10), nrows=1, ncols=2)

    # Plot transcripts
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

    # Plot cell boundaries colored by cluster labels
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
    plt.savefig(out_dir / "muspan_cosmx_visualization.png", dpi=300)
    plt.close()
