"""MuSpan CosMx handling module."""
# needs to be reviewed in detail

import warnings
import logging
from logging import getLogger
from pathlib import Path
import pickle
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from skimage import io, measure

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
):
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

    # create dictionary with format {fov: set of selected local cell IDs}
    selected_cells_by_fov = {}

    # Define Columns
    cell_col = "cell_ID"
    fov_col = "fov"

    # get unique fovs involved in selection
    fov_list = df_sel[fov_col].astype(str).unique().tolist()

    for f in fov_list:
        cells_in_fov = df_sel[df_sel["fov"].astype(str) == str(f)][cell_col].astype(str)
        cells_in_fov = set(cells_in_fov)
        selected_cells_by_fov[f] = cells_in_fov
        logger.info(f"FOV {f}: Identified {len(cells_in_fov)} valid local cells.")

    if not fov_list:
        raise ValueError("No FOV provided via fov_id or CSV.")

    all_trans_coords, all_trans_targets = [], []
    all_cell_coords, all_cell_ids, all_cell_types = [], [], []
    all_geometries, all_boundary_types = [], []

    for current_fov in fov_list:
        logger.info(f"--- Processing FOV {current_fov} ---")
        local_selected_cells = selected_cells_by_fov.get(current_fov, None)

        points_key = f"{current_fov}_points"
        if points_key not in sdata.points:
            continue

        df_pts_dask = sdata.points[points_key]
        df_all_pts = df_pts_dask.compute()

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
            all_trans_coords.append(df_pts_plot[["x_global_px", "y_global_px"]].values)
            all_trans_targets.append(df_pts_plot["target"].astype(str).values)

        if "cell_ID" in df_all_pts.columns and "x_global_px" in df_all_pts.columns:
            valid_pts = df_all_pts[df_all_pts["cell_ID"].astype(str) != "0"].copy()
            valid_pts["cell_ID_str"] = (
                valid_pts["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
            )
            global_centers = valid_pts.groupby("cell_ID_str")[
                ["x_global_px", "y_global_px"]
            ].mean()
        else:
            global_centers = pd.DataFrame()

        sdata_table = sdata.tables["table"]
        mask = sdata_table.obs["fov"].astype(str) == str(current_fov)
        adata_fov = sdata_table[mask]

        if adata_fov.n_obs == 0:
            continue

        valid_local_ids = set(
            adata_fov.obs["cell_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
        )

        padded_fov = f"F{int(current_fov):03d}"
        label_dir = Path(flat_files_dir) / "CellLabels"
        mask_files = list(label_dir.glob(f"*{padded_fov}*.tif")) + list(
            label_dir.glob(f"*{padded_fov}*.png")
        )

        if mask_files:
            logger.info(
                "Extracting boundaries and auto-aligning to global coordinate space..."
            )

            cache_file = out_dir / f"fov_{current_fov}_geometry_local_cache.pkl"

            if cache_file.exists():
                with open(cache_file, "rb") as f:
                    cache_data = pickle.load(f)

                for geom, centroid, lid in zip(
                    cache_data["geometries"],
                    cache_data["centroids"],
                    cache_data["local_ids"],
                ):
                    if local_selected_cells is None or lid in local_selected_cells:
                        unique_cell_id = f"{current_fov}_{lid}"
                        c_type = cell_id_to_type.get(
                            (str(current_fov), str(lid)), "Unknown"
                        )

                        all_geometries.append(geom)
                        all_cell_coords.append(centroid)
                        all_cell_ids.append(unique_cell_id)
                        all_cell_types.append(c_type)
                        all_boundary_types.append(c_type)
            else:
                mask_img = io.imread(mask_files[0])
                props = measure.regionprops(mask_img)

                raw_data = {}
                matched_cx, matched_cy, matched_gx, matched_gy = [], [], [], []

                for prop in props:
                    cid_str = str(prop.label)  # Local ID (e.g. "1")
                    if cid_str not in valid_local_ids:
                        continue

                    cell_mask = np.pad(
                        prop.image, pad_width=1, mode="constant", constant_values=False
                    )
                    contours = measure.find_contours(cell_mask, 0.5)
                    if not contours:
                        continue
                    contour = max(contours, key=len)

                    min_r, min_c, _, _ = prop.bbox
                    local_y = contour[:, 0] - 1 + min_r
                    local_x = contour[:, 1] - 1 + min_c
                    cx, cy = np.mean(local_x), np.mean(local_y)

                    raw_data[cid_str] = {
                        "local_x": local_x,
                        "local_y": local_y,
                        "cx": cx,
                        "cy": cy,
                    }

                    if cid_str in global_centers.index:
                        matched_cx.append(cx)
                        matched_cy.append(cy)
                        matched_gx.append(global_centers.loc[cid_str, "x_global_px"])
                        matched_gy.append(global_centers.loc[cid_str, "y_global_px"])

                if len(matched_cx) >= 3:
                    corr_x = np.corrcoef(matched_cx, matched_gx)[0, 1]
                    corr_y = np.corrcoef(matched_cy, matched_gy)[0, 1]
                    flip_x = corr_x < 0 if not np.isnan(corr_x) else False
                    flip_y = corr_y < 0 if not np.isnan(corr_y) else False

                    ox = [
                        gx + cx if flip_x else gx - cx
                        for cx, gx in zip(matched_cx, matched_gx)
                    ]
                    oy = [
                        gy + cy if flip_y else gy - cy
                        for cy, gy in zip(matched_cy, matched_gy)
                    ]
                    fov_offset_x = np.median(ox)
                    fov_offset_y = np.median(oy)
                else:
                    flip_x, flip_y, fov_offset_x, fov_offset_y = False, False, 0, 0
                    logger.warning(
                        f"FOV {current_fov}: Not enough matched cells. Using local coords."
                    )

                cache_geoms, cache_cents, cache_local_ids = [], [], []
                (
                    current_fov_geometries,
                    current_fov_coords,
                    current_fov_ids,
                    current_fov_types,
                ) = [], [], [], []

                for cid_str, data in raw_data.items():
                    lx, ly, lcx, lcy = (
                        data["local_x"],
                        data["local_y"],
                        data["cx"],
                        data["cy"],
                    )

                    global_x = (-lx if flip_x else lx) + fov_offset_x
                    global_y = (-ly if flip_y else ly) + fov_offset_y
                    poly_coords = np.column_stack((global_x, global_y))

                    gcx = (-lcx if flip_x else lcx) + fov_offset_x
                    gcy = (-lcy if flip_y else lcy) + fov_offset_y

                    if not np.array_equal(poly_coords[0], poly_coords[-1]):
                        poly_coords = np.vstack((poly_coords, poly_coords[0]))

                    if len(poly_coords) >= 3:
                        cache_geoms.append(poly_coords)
                        cache_cents.append([gcx, gcy])
                        cache_local_ids.append(cid_str)

                        if (
                            local_selected_cells is None
                            or cid_str in local_selected_cells
                        ):
                            unique_cell_id = f"{current_fov}_{cid_str}"
                            c_type = cell_id_to_type.get(
                                (str(current_fov), str(cid_str)), "Unknown"
                            )

                            current_fov_geometries.append(poly_coords)
                            current_fov_coords.append([gcx, gcy])
                            current_fov_ids.append(unique_cell_id)
                            current_fov_types.append(c_type)

                with open(cache_file, "wb") as f:
                    pickle.dump(
                        {
                            "geometries": cache_geoms,
                            "centroids": cache_cents,
                            "local_ids": cache_local_ids,
                        },
                        f,
                    )

                if current_fov_geometries:
                    all_geometries.extend(current_fov_geometries)
                    all_cell_coords.extend(current_fov_coords)
                    all_cell_ids.extend(current_fov_ids)
                    all_cell_types.extend(current_fov_types)
                    all_boundary_types.extend(current_fov_types)

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
    domain, qTrans, qCells, qBoundaries, cluster_labels, out_dir
):
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

    # Plot cell centroids
    if qCells is not None:
        ms.visualise.visualise(
            domain,
            ax=ax[1],
            objects_to_plot=qCells,
            marker_size=2,
            color_by=("constant", "#000000"),
        )

    ax[1].set_title(f"Cell Boundaries by {cluster_labels}")

    plt.tight_layout()
    plt.savefig(out_dir / "muspan_cosmx_visualization.png", dpi=300)
    plt.close()
