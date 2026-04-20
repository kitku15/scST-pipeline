"""Image viewing module."""

import scanpy as sc
import squidpy as sq
from pathlib import Path 


def run_view_images(prev_module_dir, module_dir, module_name, gene_list):
    """Run the image viewing module."""
    
    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    sc.settings.figdir = module_dir # set the figures dir to not be figures 


    # Import data
    print("Loading Xenium data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # View plots
    print("Visualize clusters on tissue...")
    sq.pl.spatial_scatter(
        adata,
        spatial_key="global",
        shape=None,
        outline=False,
        color=["leiden", "total_counts"],
        wspace=0.4,
        size=1,
        save="leiden_clusters.png",
        dpi=300,
    )
    print(f"Saved leiden clusters plot to {module_dir / 'leiden_clusters.png'}")

    # View specific gene expression
    print("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        spatial_key="global",
        color=gene_list,
        shape=None,
        size=2,
        img=False,
        save="gene_expression.png",
    )
    print(f"Saved gene expression plot to {module_dir / 'gene_expression.png'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    print(f"Data saved to {module_dir / 'adata.h5ad'}")
    print("Imaging module completed successfully.")


# if __name__ == "__main__":
    # prev_module_dir = Path('analysis/3_Annotate')
    # module_name = "4_ViewImages" 
    # module_dir = Path(f'analysis/{module_name}')
    
    # gene_list = [ # List of genes to visualize on tissue
    #     "SQSTM1",
    #     "CD74",
    #     "IGHG1",
    #     "COL3A1",
    #     "COL4A2",
    # ]

    # run_view_images(prev_module_dir, module_dir, module_name, gene_list)
