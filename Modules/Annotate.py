"""Annotation module."""

import re
import warnings
from logging import getLogger
from pathlib import Path

import pandas as pd
import scanpy as sc
from CellAnnotation_CellTypist import run_CellTypist
from CellAnnotation_plotting import run_CellType_plotting
from CellAnnotation_ScType import run_ScType
from lists import CellTypist_models, ScType_tissuetypes

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def cluster_DE_analysis(
    adata, cluster_col, module_dir, method=None, celltype_col=None, plot=True
):
    """
    cluster_col: the groupings / clustering column that the DE analysis will be based on
    """
    if method:
        output_dir = Path(f"{module_dir}/{method}/{cluster_col}/DE_analysis")
        sc.settings.figdir = output_dir

    else:
        output_dir = Path(f"{module_dir}/DE_analysis/{cluster_col}")
        sc.settings.figdir = output_dir
        celltype_col = cluster_col

    # Drop unused categories in the metadata column if it is categorical
    if hasattr(adata.obs[celltype_col], "cat"):
        adata.obs[celltype_col] = adata.obs[celltype_col].cat.remove_unused_categories()

    # 2. Remove dendrogram metadata if present
    dendro_key = f"dendrogram_{celltype_col}"
    if dendro_key in adata.uns:
        del adata.uns[dendro_key]

    # Annotate cell clusters

    # Calculate the differentially expressed genes for every cluster,
    # compared to the rest of the cells in our adata
    logger.info("Calculating differentially expressed genes for each cluster...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=celltype_col,
        method="wilcoxon",
        use_raw=True,
    )

    # 1. Plot differentially expressed genes for each cluster
    if plot:
        logger.info(
            "Plotting the top differentially expressed genes for each cluster..."
        )
        sc.pl.rank_genes_groups_dotplot(
            adata,
            groupby=celltype_col,
            standard_scale="var",
            n_genes=5,
            show=False,
            save=f"{celltype_col}.png",
        )
        logger.info(f"Dotplot saved to {sc.settings.figdir}")

        logger.info(
            "Plot differentially expressed genes for each cluster in elbow plot..."
        )
        sc.pl.rank_genes_groups(
            adata,
            n_genes=10,
            ncols=3,
            legend_fontsize=10,
            show=False,
            save=f"_{celltype_col}.png",
        )
        logger.info(f"DE Analysis plots saved to {sc.settings.figdir}")

    # Make a dataframe of marker expression
    logger.info("Save files for differentially expressed genes for each cluster...")
    logger.info("File 1...")
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)]
    markers.to_excel(
        output_dir / f"markers_{celltype_col}.xlsx",
        index=False,
    )
    logger.info(f"Markers saved to {module_dir}")

    logger.info("File 2...")

    # Define the number of clusters
    # clusters_list = len(adata.obs[cluster_col].astype(str).unique())

    # Create a list
    rows_list = []
    # Get the actual group names from the DE results
    result_groups = adata.uns["rank_genes_groups"]["names"].dtype.names

    for group in result_groups:
        top_genes = adata.uns["rank_genes_groups"]["names"][group][:10].tolist()
        new_row = pd.Series({"Cluster Name": group, "Top Genes": top_genes})
        rows_list.append(new_row)

    # Convert list of series to DataFrame
    diff_gene_df = pd.concat(rows_list, axis=1).T
    diff_gene_df.set_index(diff_gene_df.columns[0], inplace=True)
    diff_gene_df.to_csv(
        output_dir / f"top_DEgenes_{celltype_col}.csv",
        index=True,
    )
    logger.info(f"Top differentially expressed genes saved to {module_dir}")

    logger.info("File 3...")
    # Create a dictionary to store DataFrames for each cluster
    cluster_dict = {}
    cluster_path = output_dir / "DEgenes"
    cluster_path.mkdir(exist_ok=True)
    for group_name in adata.obs[celltype_col].unique():
        current_cluster = markers[markers["group"] == str(group_name)].sort_values(
            by="logfoldchanges", ascending=False
        )
        cluster_dict[f"cluster_{group_name}"] = current_cluster

        # Clean the name for the filename (remove spaces/slashes)
        safe_name = re.sub(r"[^\w\s-]", "", str(group_name)).replace(" ", "_")
        csv_filename = cluster_path / f"cluster_{safe_name}_data.csv"

        current_cluster.to_csv(csv_filename, index=False)
        logger.info(f"Exported cluster {group_name} data to {csv_filename}")

    # Rename the clusters only if ScType was falsebcs then cluster names would be integers
    col = adata.obs[celltype_col]
    is_integer_clusters = col.astype(str).str.fullmatch(r"\d+").all()

    if is_integer_clusters:
        logger.info("Integer clusters detected → renaming clusters")

        unique_clusters = col.astype(str).unique()  # Get unique clusters

        cluster_names = {  # Get unique cluster names
            cluster: f"Cluster_{cluster}" for cluster in unique_clusters
        }

        new_clusters_col = f"named_{celltype_col}"
        # Map the cluster names to the new column
        adata.obs[new_clusters_col] = col.astype(str).map(cluster_names)

    else:
        logger.info("Non-integer clusters detected → skipping renaming")


def mode_all(adata):
    umap_keys = [k for k in adata.obsm.keys() if "X_umap_n" in k]

    leiden_keys = [k for k in adata.obs.columns if k.startswith("leiden_n")]

    cluster_cols = []

    for umap in umap_keys:
        match = re.search(r"n\d+", umap)
        if match:
            n_part = match.group()
            for leiden in leiden_keys:
                if n_part in leiden:
                    cluster_cols.append((umap, leiden))

    return cluster_cols


def run_annotate(
    datatype,
    module_dir,
    cluster_name,
    input_adata_path,
    sample_key,
    ScType_anno=False,
    ScType_tissue=None,
    ScType_custom_db=None,
    ScType_mode="All",
    CellTypist_anno=False,
    CellTypist_model=None,
    CellTypist_mode="All",
    CellTypist_custom_model=None,
    CellTypist_train=False,
    CellTypist_train_data=None,
    CellTypist_train_labels=None,
    plot=True,
):
    """Run annotation."""

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading data...")
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # user selects to do ScType automatic cell type annotation
    if ScType_anno:
        # Informative error message if specified tissue type is not in ScType DB
        if ScType_tissue not in ScType_tissuetypes and not ScType_custom_db:
            message = (
                f"Invalid tissue_type: {ScType_tissue}\n"
                "Please select one of the following:\n- "
                + "\n- ".join(ScType_tissuetypes)
            )

            logger.error(message)
            raise ValueError(message)
        else:  # continue with ScType annotation
            if ScType_mode == "All":  # Proceed with ScType Annotation for ALL clusters
                # Build UMAP to leiden clustering tuple
                cluster_cols = mode_all(adata)

                # Run ScType automatic annotation
                for umap_col, cluster_col in cluster_cols:
                    # 1. Run ScType, make predictions
                    adata, sctype_column = run_ScType(
                        adata, cluster_col, ScType_tissue, ScType_custom_db
                    )

                    # 2. Cell Type plotting
                    if plot:
                        run_CellType_plotting(
                            "ScType",
                            datatype,
                            module_dir,
                            adata,
                            sample_key,
                            cluster_col,
                            sctype_column,
                            umap_col,
                        )

                    # 3. DE analysis
                    cluster_DE_analysis(
                        adata,
                        cluster_col,
                        module_dir,
                        method="ScType",
                        celltype_col=sctype_column,
                        plot=plot,
                    )

            else:  # only do ScType Annotation for one clustering
                # extract the resolution block
                parts = cluster_name.split("_")  # ['leiden', 'n10', 'r0.1']
                n_part = parts[1]  # 'n10'
                umap_col = f"X_umap_{n_part}"

                # 1. Run ScType, make predictions
                adata, sctype_column = run_ScType(
                    adata, cluster_name, ScType_tissue, ScType_custom_db
                )

                # 2. Cell Type plotting
                if plot:
                    run_CellType_plotting(
                        "ScType",
                        datatype,
                        module_dir,
                        adata,
                        sample_key,
                        cluster_name,
                        sctype_column,
                        umap_col,
                    )

                # 3. DE analysis
                cluster_DE_analysis(
                    adata,
                    cluster_name,
                    module_dir,
                    method="ScType",
                    celltype_col=sctype_column,
                    plot=plot,
                )

    if CellTypist_anno:
        if (
            not CellTypist_custom_model
            and not CellTypist_train
            and CellTypist_model not in CellTypist_models
        ):
            message = (
                f"Invalid model selection: {CellTypist_model}\n"
                "Please select one of the following:\n- "
                + "\n- ".join(CellTypist_models)
            )

            logger.error(message)
            raise ValueError(message)
        else:  # user selects model thats available, continue with CellTypist annotation
            if (
                CellTypist_mode == "All"
            ):  # Proceed with CellTypist Annotation for ALL clusters
                # Build UMAP to leiden clustering tuple
                cluster_cols = mode_all(adata)

                # Run ScType automatic annotation
                for umap_col, cluster_col in cluster_cols:
                    # 1. Run CellTypist, make predictions
                    adata, indivcellanno_col, majorvotingcellanno_col = run_CellTypist(
                        adata,
                        CellTypist_model,
                        cluster_col,
                        custom_model_path=CellTypist_custom_model,
                        train_model=CellTypist_train,
                        train_data_path=CellTypist_train_data,
                        train_labels_col=CellTypist_train_labels,
                    )
                    # 2. Cell Type plotting
                    if plot:
                        run_CellType_plotting(
                            "CellTypist",
                            datatype,
                            module_dir,
                            adata,
                            sample_key,
                            cluster_col,
                            majorvotingcellanno_col,
                            umap_col,
                            indivcellanno_col,
                        )

                    # 3. DE analysis
                    cluster_DE_analysis(
                        adata,
                        cluster_col,
                        module_dir,
                        method="CellTypist",
                        celltype_col=majorvotingcellanno_col,
                        plot=plot,
                    )

            else:  # only do CellTypist Annotation for one clustering
                # extract the resolution block
                parts = cluster_name.split("_")  # ['leiden', 'n10', 'r0.1']
                n_part = parts[1]  # 'n10'
                umap_col = f"X_umap_{n_part}"

                # 1. Run CellTypist, make predictions
                adata, indivcellanno_col, majorvotingcellanno_col = run_CellTypist(
                    adata,
                    CellTypist_model,
                    cluster_name,
                    custom_model_path=CellTypist_custom_model,
                    train_model=CellTypist_train,
                    train_data_path=CellTypist_train_data,
                    train_labels_col=CellTypist_train_labels,
                )

                # 2. Cell Type plotting
                if plot:
                    run_CellType_plotting(
                        "CellTypist",
                        datatype,
                        module_dir,
                        adata,
                        sample_key,
                        cluster_name,
                        majorvotingcellanno_col,
                        umap_col,
                        indivcellanno_col,
                    )

                # 3. DE analysis
                cluster_DE_analysis(
                    adata,
                    cluster_name,
                    module_dir,
                    method="CellTypist",
                    celltype_col=majorvotingcellanno_col,
                    plot=plot,
                )

    # If using scANVI (or pre-computed labels) and skipping ScType/CellTypist
    if not ScType_anno and not CellTypist_anno:
        logger.info(f"Using pre-computed labels from '{cluster_name}' (e.g., scANVI).")

        # 1. Plotting
        # Assuming UMAP was generated in DimReduc as X_umap_n20, find the first available UMAP
        umap_col = next((k for k in adata.obsm.keys() if "X_umap" in k), "X_umap")
        umap_col = umap_col.removeprefix("X_")  # Remove 'X_' prefix for scanpy plotting

        if plot:
            run_CellType_plotting(
                method="PreAnnotated",
                datatype=datatype,
                module_dir=module_dir,
                adata=adata,
                sample_key=sample_key,
                cluster_col=cluster_name,
                celltype_col=cluster_name,
                umap_col=umap_col,
            )

        # 2. DE analysis
        cluster_DE_analysis(
            adata,
            cluster_name,
            module_dir,
            method="PreAnnotated",
            celltype_col=cluster_name,
            plot=plot,
        )

    # Save anndata object (at the end)
    out_path = module_dir / input_adata_path.name
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
    logger.info("Annotation module completed successfully.")
