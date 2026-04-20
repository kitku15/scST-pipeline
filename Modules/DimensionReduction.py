"""Dimension reduction module."""

import scanpy as sc
import squidpy as sq
from pathlib import Path 
from config import settings, get_module


def run_dimension_reduction(data_type, prev_module_dir, module_dir, module_name, n_comps, n_neighbors, resolution, cluster_name):

    """Run dimension reduction on CosMx data."""

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial" 

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Import data
    print("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Perform dimension reduction analysis
    print("Compute PCA...")
    sc.pp.pca(adata, n_comps=n_comps)  # compute principal components
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=50,
        show=False,
        save=f"_{module_name}.png",
    )
    print(f"PCA Variance plot saved to {sc.settings.figdir}")

    print("Compute neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=15)  # compute a neighborhood graph

    print("Create UMAPs and cluster cells..")
    sc.tl.umap(adata)  # calculate umap
    sc.tl.leiden(
        adata,
        resolution=resolution,  # choose resolution for clustering
        key_added=cluster_name,
    )  # name leiden clusters

    # plot UMAP
    print("Plotting UMAPs...")
    sc.pl.umap(
        adata,
        color=[
            "total_counts",
            "n_genes_by_counts",
            cluster_name,
        ],
        wspace=0.4,
        show=False,
        save=f"_{module_name}.png",  # save the figure with the module name
        frameon=False,
    )
    print(f"UMAP plot saved to {sc.settings.figdir}")

    # plot visualization of leiden clusters
    print(f"Plotting {cluster_name} clusters...")
    # Create a plot where each FOV is its own panel

    sq.pl.spatial_scatter(
        adata,
        color=[cluster_name],
        spatial_key=spatial_key, 
        shape=None,
        size=0.05,             # Use a very small size for the full slide
        alpha=0.6,
        frameon=False,
        img=True, 
        save=f"{cluster_name}_full_stitched.png",
    )

    if data_type == "CosMx":
        sq.pl.spatial_scatter(
            adata,
            color=[cluster_name],
            library_key="fov",     # Use the 'fov' column from your adata.obs
            ncols=4,               # Arrange in 4 columns
            shape=None,            # Circles
            size=1,                # Adjust size if dots are too big/small
            img=True,             # Keep False until we confirm coordinates are right
            save=f"{cluster_name}_by_fov.png"
        )
        
    print(f"{cluster_name} spatial scatter plot saved to {module_dir}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    print(f"Data saved to {module_dir / 'adata.h5ad'}")


if __name__ == "__main__":
    data_type = settings['project']['data_type']

    module_1_name, module_1_dir = get_module(1)
    module_2_name, module_2_dir = get_module(2)
   
    n_comps = settings['modules']['DimensionReduction']['n_comps']                    
    n_neighbors = settings['modules']['DimensionReduction']['n_neighbors']                       
    resolution = settings['modules']['DimensionReduction']['resolution']                      
    cluster_name = settings['modules']['DimensionReduction']['cluster_name']   

    run_dimension_reduction(data_type, module_1_dir, module_2_dir, module_2_name, n_comps, n_neighbors, resolution, cluster_name)