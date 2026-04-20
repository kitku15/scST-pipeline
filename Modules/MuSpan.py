import warnings
from logging import getLogger
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import spatialdata as sd
import numpy as np
import muspan as ms
from skimage import io, measure
import pickle

warnings.filterwarnings("ignore")
logger = getLogger(__name__)

import logging
logging.basicConfig(level=logging.INFO)

def run_cosmx_muspan(
    zarr_path: str, 
    flat_files_dir: str, 
    fov_id: str, 
    module_5_dir: str, 
    module_dir: str, 
    cluster_labels: str = "cell_type",
    transcripts_to_keep: list = None
):
    """Run Muspan on CosMx SpatialData for a specific FOV."""
    
    out_dir = Path(module_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Load Data
    logger.info("Loading CosMx SpatialData and Annotated AnnData...")
    sdata = sd.read_zarr(zarr_path)

    adata_path = Path(module_5_dir) / "adata.h5ad"
    external_adata = sc.read_h5ad(adata_path) 
    
    domain_name = f"CosMx_FOV_{fov_id}"
    domain = ms.domain(domain_name)


    # 2. Add Transcripts to Domain
    points_key = f"{fov_id}_points"
    logger.info(f"Loading transcripts from {points_key}...")
    
    df_pts = sdata.points[points_key].compute()
    
    if transcripts_to_keep is not None:
        logger.info(f"Filtering down to {len(transcripts_to_keep)} specific transcripts...")
        df_pts['target'] = df_pts['target'].astype(str)
        df_pts = df_pts[df_pts['target'].isin(transcripts_to_keep)]
    
    
    meta_df = pd.read_csv(Path(flat_files_dir) / 'Quarter_metadata_file.csv') 
    center_y_global = meta_df[meta_df['fov'] == int(fov_id)]['CenterY_global_px'].iloc[0]

    # --- Flip transcript Y-coordinates (Horizontal Mirror) ---
    # df_pts['y_global_px'] = (2 * center_y_global) - df_pts['y_global_px']
    # ----------------------------------------------------------------------
    
    trans_coords = df_pts[['x_global_px', 'y_global_px']].values
    domain.add_points(trans_coords, collection_name="Transcripts")
    
    qTrans = ms.query.query(domain, ('collection',), 'is', 'Transcripts')
    domain.add_labels('Transcript ID', df_pts['target'].values, add_labels_to=qTrans)


    # 3. Add Cell Centroids
    logger.info("Loading Cell Centroids...")
    sdata_table = sdata.tables['table']
    mask = sdata_table.obs['fov'].astype(str) == str(fov_id)
    adata_fov = sdata_table[mask].copy()
    
    # Rather than using obsm['spatial'] which seems corrupted (Y=0), 
    # let's pull directly from the global columns usually present in CosMx Anndata obs
    if 'CenterX_global_px' in adata_fov.obs.columns:
        x_cents = adata_fov.obs['CenterX_global_px'].values
        y_cents = adata_fov.obs['CenterY_global_px'].values
        
        # APPLY Y-INVERSION TO MATCH TRANSCRIPTS
        y_cents = -y_cents 
        cell_coords = np.column_stack((x_cents, y_cents))
    else:
        # Fallback to obsm, but apply the inversion
        cell_coords = adata_fov.obsm['spatial'].copy()
        # cell_coords[:, 1] = -cell_coords[:, 1] # Flip Y to match transcripts

    domain.add_points(cell_coords, collection_name="Cell centroids")
    
    qCells = ms.query.query(domain, ('collection',), 'is', 'Cell centroids')
    cell_ids = adata_fov.obs['cell_ID'].values.astype(str)
    domain.add_labels('Cell ID', cell_ids, add_labels_to=qCells)


    # 4. Map Cell Types to Centroids
    logger.info("Mapping annotated cell types to centroids...")
    cell_id_to_type = dict(zip(external_adata.obs_names.astype(str), external_adata.obs[cluster_labels]))
    
    sdata_indices = adata_fov.obs_names.astype(str)
    cell_types_ordered = [cell_id_to_type.get(idx, "Filtered_in_QC") for idx in sdata_indices]
    domain.add_labels(cluster_labels, cell_types_ordered, add_labels_to=qCells)

    # ==========================================
    # 5. GENERATE POLYGONS FROM IMAGE MASKS
    # ==========================================
    logger.info("Extracting boundaries from CellLabels images...")
    flat_dir = Path(flat_files_dir)
    
    meta_fov = meta_df[meta_df['fov'] == int(fov_id)]
    
    if meta_fov.empty:
        raise ValueError(f"FOV {fov_id} not found in metadata.csv!")
        
    offset_x = meta_fov['CenterX_global_px'].iloc[0] - meta_fov['CenterX_local_px'].iloc[0]
    offset_y = meta_fov['CenterY_global_px'].iloc[0] - meta_fov['CenterY_local_px'].iloc[0]
    
    padded_fov = f"F{int(fov_id):03d}"
    label_dir = flat_dir / "CellLabels"
    mask_files = list(label_dir.glob(f"*{padded_fov}*.tif")) + list(label_dir.glob(f"*{padded_fov}*.png"))
    
    if not mask_files:
        raise FileNotFoundError(f"Could not find segmentation mask for FOV {fov_id}")
        
    mask_img = io.imread(mask_files[0])
    
    # --- GET THE HORIZONTAL CENTER OF THE FOV ---
    fov_center_y = mask_img.shape[0] / 2.0
    # --------------------------------------------

    cache_file = out_dir / f"fov_{fov_id}_geometry_cache.pkl"

    if cache_file.exists():
        logger.info(f"Loading cached boundaries...")
        with open(cache_file, 'rb') as f:
            cache_data = pickle.load(f)
        geometries = cache_data['geometries']
        poly_cell_ids = cache_data['poly_cell_ids']
    else:
        logger.info("Tracing boundaries... (This will be cached)")
        geometries = []
        poly_cell_ids = []
        unique_cells = np.unique(mask_img)

        for cid in unique_cells:
            if cid == 0: continue 
            
            cell_mask = (mask_img == cid)
            contours = measure.find_contours(cell_mask, 0.5)
            if not contours: continue
            
            contour = max(contours, key=len)
            local_y, local_x = contour[:, 0], contour[:, 1]
            
            # --- FLIP POLYGONS ON HORIZONTAL MIRROR FROM CENTER ---
            local_y = (2 * fov_center_y) - local_y
            # ------------------------------------------------------

            poly_coords = np.column_stack((local_x, local_y))

            if not np.array_equal(poly_coords[0], poly_coords[-1]):
                poly_coords = np.vstack((poly_coords, poly_coords[0]))
            
            if len(poly_coords) >= 3:
                geometries.append(poly_coords)
                poly_cell_ids.append(str(cid))

        with open(cache_file, 'wb') as f:
            pickle.dump({'geometries': geometries, 'poly_cell_ids': poly_cell_ids}, f)
    
    # Add to MuSpAn
    domain.add_shapes(shapes=geometries, collection_name="Cell boundaries")
    qBoundaries = ms.query.query(domain, ('collection',), 'is', 'Cell boundaries')

    # Map labels to the newly generated polygons
    raw_id_to_type = dict(zip(cell_ids, cell_types_ordered))
    boundary_types_ordered = [raw_id_to_type.get(cid, "Filtered_in_QC") for cid in poly_cell_ids]
    domain.add_labels(cluster_labels, boundary_types_ordered, add_labels_to=qBoundaries)

    # ==========================================
    # 6. Visualizations
    # ==========================================
    logger.info(f"Visualizing the MuSpAn domain: {domain_name}")
    fig, ax = plt.subplots(figsize=(20, 10), nrows=1, ncols=2)
    
    # Left Plot: Transcripts
    ms.visualise.visualise(
        domain,
        color_by=("label", "Transcript ID"),
        ax=ax[0],
        objects_to_plot=qTrans,
        marker_size=0.5, 
    )
    ax[0].set_title(f"Filtered Transcripts (FOV {fov_id})")

    # Right Plot: Cell Boundaries colored by cluster + black centroids
    # Draw the filled polygons
    ms.visualise.visualise(
        domain,
        color_by=("label", cluster_labels),
        ax=ax[1],
        objects_to_plot=qBoundaries
    )
    # Overlay the centroids as small black dots
    ms.visualise.visualise(
        domain,
        ax=ax[1],
        objects_to_plot=qCells,
        marker_size=2,
        color_by=('constant', '#000000') ,
    )
    ax[1].set_title(f"Cell Boundaries by {cluster_labels} (FOV {fov_id})")

    plt.tight_layout()
    plt.savefig(out_dir / f"muspan_fov_{fov_id}_visualization.png", dpi=300)
    logger.info("Visualization saved.")

    # Save Domain
    ms.io.save_domain(
        domain,
        name_of_file=f"muspan_object_fov_{fov_id}",
        path_to_save=str(out_dir),
    )
    logger.info("Domain saved successfully.")
    
    return domain

# if __name__ == "__main__":
    
#     ZARR_PATH = "Lung13.zarr"
#     FLAT_FILES_DIR = r"C:\Users\bunga\python\Project2\CosMx\Lung13+SMI+Flat+data\Lung13\Lung13-Flat_files_and_images" 
    
#     FOV_TO_LOAD = "8" 
#     ADATA_PATH = "analysis/5_SpatialStat/adata.h5ad" 
#     OUTPUT_DIR = "analysis/6_MuSpan"
    
#     TRANSCRIPTS = [ # List of genes to visualize on tissue
#         "SQSTM1",
#         "CD74",
#         "IGHG1",
#         "COL3A1",
#         "COL4A2",
#     ]
    
#     run_cosmx_muspan(
#         zarr_path=ZARR_PATH,
#         flat_files_dir=FLAT_FILES_DIR,
#         fov_id=FOV_TO_LOAD,
#         adata_path=ADATA_PATH,
#         output_dir=OUTPUT_DIR,
#         cluster_labels="cell_type", 
#         transcripts_to_keep=TRANSCRIPTS
#     )