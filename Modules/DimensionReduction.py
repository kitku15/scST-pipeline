"""Dimension reduction module."""

import scanpy as sc
import squidpy as sq
from pathlib import Path 


def run_dimension_reduction(prev_module_dir, module_dir, module_name, n_comps, n_neighbors, resolution, cluster_name):

    """Run dimension reduction on CosMx data."""

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Import data
    print("Loading CosMx data...")
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
        spatial_key="global",  # <--- THIS IS THE MAGIC SWITCH
        shape=None,
        size=0.05,             # Use a very small size for the full slide
        alpha=0.6,
        frameon=False,
        img=True, 
        save=f"{cluster_name}_full_stitched.png",
    )
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
    prev_module_dir = Path('analysis/1_QualityControl')
    module_dir = Path('analysis/2_DimensionReduction')
    module_name = "2_dimension_reduction" # Name of the module - will be used in the output directory name
    n_comps = 50                          # number of principal components to compute
    n_neighbors = 15                      # number of neighbors for the neighborhood graph
    resolution = 0.2                    # resolution for leiden clustering
    cluster_name = "leiden"    
   

    run_dimension_reduction(prev_module_dir, module_dir, module_name, n_comps, n_neighbors, resolution, cluster_name)