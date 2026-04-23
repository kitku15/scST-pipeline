import warnings
import logging
from logging import getLogger
from pathlib import Path
import pickle

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

import spatialdata as sd
from skimage import io, measure
from config import settings, get_module

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


def map_cell_types_to_domain_xenium(adata, domain, adata_cell_id, cluster_labels):
    """Helper function to map AnnData cell types to a Muspan domain (Used for Xenium)"""
    domain_cell_ids_ordered = [
        str(cell_id) for cell_id in domain.labels["Cell ID"]["labels"]
    ]
    domain_cell_ids_unique = set(domain_cell_ids_ordered)

    logger.info(f"Number of unique cells in the domain: {len(domain_cell_ids_unique)}")

    # Filter adata to include only cells in the area of interest
    filt_adata = adata[adata.obs[adata_cell_id].isin(domain_cell_ids_unique)]

    logger.info(f"Filtered adata from {adata.n_obs} to {filt_adata.n_obs} cells")
    logger.info("Adding cell_type IDs to domain with cluster labels")

    cell_id_to_type = dict(
        zip(filt_adata.obs[adata_cell_id], filt_adata.obs[cluster_labels])
    )
    cell_types_ordered = [
        cell_id_to_type.get(cell_id, "Unknown") for cell_id in domain_cell_ids_ordered
    ]

    domain.add_labels(label_name=cluster_labels, labels=cell_types_ordered)
    return domain


def run_muspan(
    dataset_type: str,
    module_dir: str,
    prev_module_dir: str,
    domain_name: str,
    cluster_labels: str = "cell_type",
    transcripts_of_interest: list = None,
    # CosMx specific kwargs
    zarr_path: str = None,
    flat_files_dir: str = None,
    fov_id: str = None,
    # Xenium specific kwargs
    xenium_dir: str = None,
    area_path: str = None,
    adata_cell_id: str = "cell_id",
):
    """
    Run Muspan spatial analysis for either CosMx or Xenium data.
    """
    dataset_type = dataset_type.lower()
    if dataset_type not in ["cosmx", "xenium"]:
        raise ValueError("dataset_type must be either 'CosMx' or 'Xenium'")

    out_dir = Path(module_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load AnnData
    adata_path = Path(prev_module_dir) / "adata.h5ad"
    logger.info(f"Loading Annotated AnnData from {adata_path}...")
    adata = sc.read_h5ad(adata_path)

    # 2. COSMX
    if dataset_type == "cosmx":
        logger.info(f"Processing CosMx data for FOV {fov_id}...")
        domain = ms.domain(domain_name)
        sdata = sd.read_zarr(zarr_path)

        # Transcripts
        points_key = f"{fov_id}_points"
        logger.info(f"Loading transcripts from {points_key}...")

        # df_pts = sdata.points[points_key].compute()

        # if transcripts_of_interest is not None:
        #     df_pts['target'] = df_pts['target'].astype(str)
        #     df_pts = df_pts[df_pts['target'].isin(transcripts_of_interest)]

        #  filter the Dask dataframe before pulling it into memory.
        df_pts_dask = sdata.points[points_key]
        if transcripts_of_interest is not None:
            # Filter lazily, THEN compute
            df_pts = df_pts_dask[
                df_pts_dask["target"].isin(transcripts_of_interest)
            ].compute()
            df_pts["target"] = df_pts["target"].astype(str)
        else:
            df_pts = df_pts_dask.compute()

        trans_coords = df_pts[["x_global_px", "y_global_px"]].values
        domain.add_points(trans_coords, collection_name="Transcripts")

        qTrans = ms.query.query(domain, ("collection",), "is", "Transcripts")
        domain.add_labels(
            "Transcript ID", df_pts["target"].values, add_labels_to=qTrans
        )

        # Cell Centroids
        logger.info("Loading Cell Centroids...")
        sdata_table = sdata.tables["table"]
        mask = sdata_table.obs["fov"].astype(str) == str(fov_id)
        adata_fov = sdata_table[mask].copy()

        if "CenterX_global_px" in adata_fov.obs.columns:
            x_cents = adata_fov.obs["CenterX_global_px"].values
            y_cents = -adata_fov.obs["CenterY_global_px"].values  # Y-INVERSION
            cell_coords = np.column_stack((x_cents, y_cents))
        else:
            cell_coords = adata_fov.obsm["spatial"].copy()
            # cell_coords[:, 1] = -cell_coords[:, 1] # Fallback Y-inversion

        domain.add_points(cell_coords, collection_name="Cell centroids")
        qCells = ms.query.query(domain, ("collection",), "is", "Cell centroids")
        cell_ids = adata_fov.obs["cell_ID"].values.astype(str)
        domain.add_labels("Cell ID", cell_ids, add_labels_to=qCells)

        # Map Labels to Centroids
        logger.info("Mapping annotated cell types to centroids...")
        cell_id_to_type = dict(
            zip(adata.obs_names.astype(str), adata.obs[cluster_labels])
        )
        sdata_indices = adata_fov.obs_names.astype(str)
        cell_types_ordered = [
            cell_id_to_type.get(idx, "Filtered_in_QC") for idx in sdata_indices
        ]
        domain.add_labels(cluster_labels, cell_types_ordered, add_labels_to=qCells)

        # Trace Polygons from Image
        logger.info("Extracting boundaries from CellLabels images...")

        padded_fov = f"F{int(fov_id):03d}"
        label_dir = Path(flat_files_dir) / "CellLabels"
        mask_files = list(label_dir.glob(f"*{padded_fov}*.tif")) + list(
            label_dir.glob(f"*{padded_fov}*.png")
        )

        mask_img = io.imread(mask_files[0])
        fov_center_y = mask_img.shape[0] / 2.0

        cache_file = out_dir / f"fov_{fov_id}_geometry_cache.pkl"
        if cache_file.exists():
            with open(cache_file, "rb") as f:
                cache_data = pickle.load(f)
            geometries = cache_data["geometries"]
            poly_cell_ids = cache_data["poly_cell_ids"]
        else:
            geometries, poly_cell_ids = [], []
            unique_cells = np.unique(mask_img)

            for cid in unique_cells:
                if cid == 0:
                    continue

                contours = measure.find_contours((mask_img == cid), 0.5)
                if not contours:
                    continue

                contour = max(contours, key=len)
                local_y, local_x = contour[:, 0], contour[:, 1]
                local_y = (2 * fov_center_y) - local_y  # HORIZONTAL MIRROR FLIP

                poly_coords = np.column_stack((local_x, local_y))
                if not np.array_equal(poly_coords[0], poly_coords[-1]):
                    poly_coords = np.vstack((poly_coords, poly_coords[0]))

                if len(poly_coords) >= 3:
                    geometries.append(poly_coords)
                    poly_cell_ids.append(str(cid))

            with open(cache_file, "wb") as f:
                pickle.dump(
                    {"geometries": geometries, "poly_cell_ids": poly_cell_ids}, f
                )

        domain.add_shapes(shapes=geometries, collection_name="Cell boundaries")
        qBoundaries = ms.query.query(domain, ("collection",), "is", "Cell boundaries")

        raw_id_to_type = dict(zip(cell_ids, cell_types_ordered))
        boundary_types_ordered = [
            raw_id_to_type.get(cid, "Filtered_in_QC") for cid in poly_cell_ids
        ]
        domain.add_labels(
            cluster_labels, boundary_types_ordered, add_labels_to=qBoundaries
        )

        # Plotting (CosMx)
        logger.info(f"Visualizing the MuSpAn domain: {domain_name}")
        fig, ax = plt.subplots(figsize=(20, 10), nrows=1, ncols=2)

        ms.visualise.visualise(
            domain,
            color_by=("label", "Transcript ID"),
            ax=ax[0],
            objects_to_plot=qTrans,
            marker_size=0.5,
        )
        ax[0].set_title(f"Filtered Transcripts (FOV {fov_id})")

        ms.visualise.visualise(
            domain,
            color_by=("label", cluster_labels),
            ax=ax[1],
            objects_to_plot=qBoundaries,
        )
        ms.visualise.visualise(
            domain,
            ax=ax[1],
            objects_to_plot=qCells,
            marker_size=2,
            color_by=("constant", "#000000"),
        )
        ax[1].set_title(f"Cell Boundaries by {cluster_labels} (FOV {fov_id})")

        plt.tight_layout()
        plt.savefig(out_dir / f"muspan_fov_{fov_id}_visualization.png", dpi=300)

    # 3. XENIUM
    elif dataset_type == "xenium":
        logger.info("Processing Xenium data...")

        # Load Xenium Domain directly via MuSpan
        domain = ms.io.xenium_to_domain(
            path_to_xenium_data=str(xenium_dir),
            cells_from_selection_csv=str(area_path),
            domain_name=domain_name,
            load_transcripts=True,
            selected_transcripts=transcripts_of_interest,
            load_nuclei=True,
            load_cells_as_shapes=True,
            exclude_no_nuclei_cells=True,
        )

        # Map labels
        domain = map_cell_types_to_domain_xenium(
            adata, domain, adata_cell_id, cluster_labels
        )

        # Queries
        qCells = ms.query.query(domain, ("Collection",), "is", "Cell boundaries")
        qTrans = ms.query.query(domain, ("Collection",), "is", "Transcripts")
        qNuc = ms.query.query(domain, ("Collection",), "is", "Nucleus boundaries")

        # Plotting 1: Overviews (Xenium)
        fig, ax = plt.subplots(figsize=(20, 15), nrows=2, ncols=2)
        ms.visualise.visualise(domain, ax=ax[0, 0], marker_size=0.05)
        ax[0, 0].set_title("All objects")

        ms.visualise.visualise(
            domain,
            color_by=("label", cluster_labels),
            ax=ax[0, 1],
            objects_to_plot=qCells,
        )
        ax[0, 1].set_title("Cell type")

        ms.visualise.visualise(
            domain,
            color_by=("label", "Transcript ID"),
            ax=ax[1, 0],
            objects_to_plot=qTrans,
            marker_size=5,
        )
        ax[1, 0].set_title("Transcripts")

        ms.visualise.visualise(
            domain,
            color_by=("label", "Nucleus Area"),
            ax=ax[1, 1],
            objects_to_plot=qNuc,
        )
        ax[1, 1].set_title("Nuclei")

        plt.tight_layout()
        plt.savefig(out_dir / "muspan_domain_visualization.png")

        # Convert boundaries to centroids
        logger.info("Convert cell boundaries to cell centres (centroids)")
        domain.convert_objects(
            population=("Collection", "Cell boundaries"),
            object_type="point",
            conversion_method="centroids",
            collection_name="Cell centroids",
            inherit_collections=False,
        )

        # Plotting 2: Centroids (Xenium)
        plt.figure(figsize=(10, 6))
        ms.visualise.visualise(
            domain,
            objects_to_plot=("collection", "Cell centroids"),
            color_by=cluster_labels,
            ax=plt.gca(),
        )
        plt.tight_layout()
        plt.savefig(out_dir / "muspan_cell_centroids.png")

        # Plotting 3: Overlay (Xenium)
        plt.figure(figsize=(10, 6))
        ms.visualise.visualise(
            domain, objects_to_plot=qCells, color_by=cluster_labels, ax=plt.gca()
        )
        ms.visualise.visualise(
            domain,
            objects_to_plot=("collection", "Cell centroids"),
            color_by=("constant", "black"),
            ax=plt.gca(),
            marker_size=1,
            add_cbar=False,
        )
        plt.tight_layout()
        plt.savefig(out_dir / "muspan_cell_centroids_n_boundaries.png")

    # 4. SAVE DOMAIN
    ms.io.save_domain(
        domain,
        name_of_file=f"muspan_object_{domain_name.replace(' ', '_')}",
        path_to_save=str(out_dir),
    )
    logger.info(f"Domain saved successfully to {out_dir}")

    return domain


if __name__ == "__main__":
    data_type = settings["project"]["data_type"]

    module_5_name, module_5_dir = get_module(5)
    module_6_name, module_6_dir = get_module(6)

    zarr_path = settings["io"]["zarr_dir"]
    dataset_path = settings["io"]["dataset_dir"]

    fov = settings["modules"]["MuSpan"]["fov"]
    transcripts = settings["modules"]["MuSpan"]["transcripts"]
    muspan_object = f"muspan_object_fov_{fov}.muspan"
    cluster_labels = settings["modules"]["MuSpan"]["cluster_labels"]
    area_path = settings["modules"]["MuSpan"]["area_path"]

    if data_type == "CosMx":
        run_muspan(
            dataset_type=data_type,
            module_dir=module_6_dir,
            prev_module_dir=module_5_dir,
            domain_name=f"{data_type}_FOV_{fov}",
            cluster_labels=cluster_labels,
            transcripts_of_interest=transcripts,
            zarr_path=zarr_path,
            flat_files_dir=dataset_path,
            fov_id=fov,
        )
    elif data_type == "Xenium":
        run_muspan(
            dataset_type=data_type,
            module_dir=module_6_dir,
            prev_module_dir=module_5_dir,
            domain_name=f"{data_type}_domain",
            cluster_labels=cluster_labels,
            transcripts_of_interest=transcripts,
            xenium_dir=dataset_path,
            area_path=area_path,
            adata_cell_id="cell_id",
        )
