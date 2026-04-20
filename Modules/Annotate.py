"""Annotation module."""


import pandas as pd
import scanpy as sc
from pathlib import Path 


def run_annotate(module_dir, module_name, cluster_name, new_clusters, prev_module_dir):
    """Run annotation on CosMx data."""

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    print("Loading CosMx data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Annotate cell clusters
    # Calculate the differentially expressed genes for every cluster,
    # compared to the rest of the cells in our adata
    print("Calculating differentially expressed genes for each cluster...")
    sc.tl.rank_genes_groups(adata, groupby=cluster_name, method="wilcoxon")

    print("Plotting the top differentially expressed genes for each cluster...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_name,
        standard_scale="var",
        n_genes=5,
        show=False,
        save=f"{module_name}.png",
    )
    print(f"Dotplot saved to {sc.settings.figdir}")

    # Plot differentially expressed genes for each cluster
    print("Plot differentially expressed genes for each cluster in elbow plot...")
    sc.pl.rank_genes_groups(
        adata,
        n_genes=10,
        ncols=3,
        legend_fontsize=10,
        show=False,
        save=f"_{module_name}.png",
    )
    print(f"UMAP plot saved to {sc.settings.figdir}")

    # Make a dataframe of marker expression
    print("Save files for differentially expressed genes for each cluster...")
    print("File 1...")
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)]
    markers.to_excel(
        module_dir / "markers.xlsx",
        index=False,
    )
    print(f"Markers saved to {module_dir}")

    print("File 2...")
    # Define the number of clusters
    clusters_list = len(adata.obs[cluster_name].astype(str).unique())

    # Create a list
    list = []
    for cluster_number in range(clusters_list):
        top_genes = adata.uns["rank_genes_groups"]["names"][
            str(cluster_number)
        ]  # Get the names of the top differentially expressed genes
        top_genes = top_genes[:10]  # Get the top 10 genes
        new_row = pd.Series(
            {"Cluster Number": cluster_number, "Top Genes": top_genes}
        )  # Create a new row with the cluster_number and top_genes
        list.append(new_row)

    # Convert list of series to DataFrame
    diff_gene_df = pd.concat(list, axis=1).T
    diff_gene_df.set_index(diff_gene_df.columns[0], inplace=True)
    diff_gene_df.to_csv(
        module_dir / "top_differentially_expressed_genes.csv",
        index=True,
    )
    print(f"Top differentially expressed genes saved to {module_dir}")

    print("File 3...")
    # Create a dictionary to store DataFrames for each cluster
    cluster_dict = {}
    cluster_path = module_dir / "cluster_diff_genes"
    cluster_path.mkdir(exist_ok=True)
    for cluster_number in range(clusters_list):
        # print(cluster_number)
        current_cluster = markers[markers["group"] == str(cluster_number)].sort_values(
            by="logfoldchanges", ascending=False
        )  # make a dataframe of the current cluster
        cluster_dict[f"cluster_{cluster_number}"] = (
            current_cluster  # Store the DataFrame in the dictionary
        )
        # Export the DataFrame to a CSV file
        csv_filename = cluster_path / f"cluster_{cluster_number}_data.csv"

        current_cluster.to_csv(csv_filename, index=False)
        print(f"Exported cluster {cluster_number} data to {csv_filename}")

    # Rename the clusters based on the markers
    print("Renaming clusters based on markers...")
    # Get unique clusters
    unique_clusters = (
        adata.obs[cluster_name].astype(str).unique()
    )  # Get unique cluster names
    cluster_names = {
        cluster: f"Cluster_{cluster}" for cluster in unique_clusters
    }  # Create a mapping of cluster names
    adata.obs[new_clusters] = (
        adata.obs[cluster_name].astype(str).map(cluster_names)
    )  # Map the cluster names to the cell_type column

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    print(f"Data saved to {module_dir / 'adata.h5ad'}")
    print("Annotation module completed successfully.")


if __name__ == "__main__":

    prev_module_dir = Path('analysis/2_DimensionReduction')
    module_name = "3_Annotate" # Name of the module - will be used in the output directory name
    module_dir = Path(f'analysis/{module_name}')
    cluster_name = "leiden"    # name of the cluster column in adata.obs
    new_clusters = "cell_type" # name of the new cluster column in adata.obs

    run_annotate(module_dir, module_name, cluster_name, new_clusters, prev_module_dir)
 