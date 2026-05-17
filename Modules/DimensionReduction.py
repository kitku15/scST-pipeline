"""Dimension reduction module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq
import matplotlib
import matplotlib.pyplot as plt
from scvi.external import SCVIVA


matplotlib.use("Agg")


warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_SCVIVA(
    adata,
    scviva_batch_key,
    spatial_key,
    scviva_spatial_knn,
    scviva_epochs,
    scviva_layer,
):
    logger.info("--- Initiating scVIVA Pipeline ---")

    # FIX: Ensure a sample_key exists for scVIVA (it crashes if sample_key is None)
    if scviva_batch_key and scviva_batch_key in adata.obs.columns:
        valid_sample_key = scviva_batch_key
    else:
        valid_sample_key = "_scviva_dummy_sample"
        adata.obs[valid_sample_key] = "sample_1"
        logger.info(
            f"No valid scviva_batch_key provided. Created dummy sample key '{valid_sample_key}' for spatial graph construction."
        )

    # scVIVA requires cell states (labels) and internal expression states (embeddings)
    # We generate a quick preliminary clustering using standard PCA to fulfill this.
    logger.info("Generating preliminary labels for scVIVA environment features...")
    sc.pp.neighbors(adata, n_pcs=30, key_added="pre_scviva_neighbors", n_neighbors=15)
    sc.tl.leiden(
        adata,
        resolution=0.5,
        key_added="pre_scviva_labels",
        neighbors_key="pre_scviva_neighbors",
    )

    setup_kwargs = {
        "labels_key": "pre_scviva_labels",
        "cell_coordinates_key": spatial_key,
        "expression_embedding_key": "X_pca",  # using PCA we generated earlier
        "sample_key": valid_sample_key,
    }

    logger.info("Preprocessing AnnData for scVIVA (computing spatial niche graphs)...")
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=scviva_spatial_knn,
        **setup_kwargs,
    )

    logger.info("Setting up AnnData for scVIVA...")
    SCVIVA.setup_anndata(
        adata,
        layer=scviva_layer,
        batch_key=valid_sample_key,
        **setup_kwargs,
    )

    logger.info("Training scVIVA model...")
    nichevae = SCVIVA(adata)
    nichevae.train(
        max_epochs=scviva_epochs,
        early_stopping=True,
        check_val_every_n_epoch=1,
        batch_size=512,
        plan_kwargs={"lr": 5e-4},
    )

    logger.info("Extracting scVIVA latent representation...")
    adata.obsm["X_scVIVA"] = nichevae.get_latent_representation()
    logger.info("--- scVIVA Pipeline Complete ---")

    return adata


def run_dimension_reduction(
    data_type,
    prev_module_dir,
    module_dir,
    module_name,
    n_comps,
    n_neighbors_list,
    resolution_list,
    cluster_name,
    scviva_layer="counts",  # Layer with raw integer counts (required for scvi)
    scviva_batch_key=None,  # Batch/Sample key if you have multiple FOVs/Samples
    scviva_spatial_knn=20,  # Number of spatial neighbors for scVIVA env graph
    scviva_epochs=400,  # Max epochs for scVIVA training
):
    """Run dimension reduction on data."""

    # Ensure inputs are lists for iteration
    if not isinstance(n_neighbors_list, list):
        n_neighbors_list = [n_neighbors_list]
    if not isinstance(resolution_list, list):
        resolution_list = [resolution_list]

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    cluster_palette_25 = [
        "#be84bf",
        "#ffff34",
        "#e41c1e",
        "#b5df6e",
        "#65c1a4",
        "#d95e01",
        "#b3b3b3",
        "#984da3",
        "#ff8045",
        "#e78ac3",
        "#2d81b9",
        "#050582",
        "#ffd92e",
        "#fb9a74",
        "#92a5cd",
        "#e6aa02",
        "#ff7f00",
        "#fb9998",
        "#f0027f",
        "#ff50a7",
        "#746fb2",
        "#199d76",
        "#8e8e8e",
        "#fc5d5d",
        "#77b975",
        "#bf5c18",
        "#36a230",
        "#4084bb",
        "#989898",
        "#b1df89",
        "#bcb8d9",
        "#a55527",
        "#cd0000",
        "#006300",
        "#ff0f0d",
        "#e5c493",
        "#fb7f71",
        "#7fc97f",
        "#7eeec9",
        "#a6d854",
        "#8dd3c7",
        "#f781bf",
        "#8b8878",
        "#fdb462",
    ]

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir
    sc.set_figure_params(
        facecolor="white", transparent=False, dpi=300, figsize=(12, 12)
    )

    # Import data
    logger.info("Loading data...")
    adata = sc.read_h5ad(prev_module_dir / "adata.h5ad")

    # Ensure layer exists for scvi
    if scviva_layer not in adata.layers:
        logger.warning(
            f"Layer '{scviva_layer}' not found. scVI needs raw counts. Using adata.X instead."
        )
        scviva_layer = None

    # Perform dimension reduction analysis
    logger.info("Compute PCA...")
    sc.pp.pca(adata, n_comps=n_comps)  # compute principal components
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=50,
        show=False,
        save="PCA.png",
    )
    logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    # run SCVIVA
    adata_path = module_dir / "adata.h5ad"

    if not adata_path.exists():
        adata = run_SCVIVA(
            adata,
            scviva_batch_key,
            spatial_key,
            scviva_spatial_knn,
            scviva_epochs,
            scviva_layer,
        )
    else:
        logger.info(
            f"Found existing scVIVA output at {adata_path}. Loading instead of re-running."
        )
        adata = sc.read_h5ad(adata_path)

    for n_neighbors in n_neighbors_list:
        logger.info(f"Compute neighbors for n={n_neighbors}...")

        # 1. Create a unique key for this neighbor graph
        neighbors_key = f"neighbors_n{n_neighbors}"
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            use_rep="X_scVIVA",  # use scVIVA latent space for neighbors
            key_added=neighbors_key,  # Save graph to unique key
        )

        logger.info(f"Create UMAPs and cluster cells for n={n_neighbors}...")
        sc.tl.umap(adata, neighbors_key=neighbors_key)

        # Scanpy saves the UMAP to 'X_umap' by default.
        # We copy it to a unique name so the next loop doesn't overwrite it
        custom_umap_basis = f"umap_n{n_neighbors}"
        adata.obsm[f"X_{custom_umap_basis}"] = adata.obsm["X_umap"].copy()

        for resolution in resolution_list:
            current_cluster_name = f"{cluster_name}_n{n_neighbors}_r{resolution}"

            # combination specific subfolder
            combo_dir = module_dir / f"n{n_neighbors}_r{resolution}"
            combo_dir.mkdir(parents=True, exist_ok=True)

            # Point Scanpy/Squidpy to save figures in this subfolder
            sc.settings.figdir = combo_dir

            # Tell Leiden to use the specific neighbors graph
            logger.info(f"Running Leiden clustering for {current_cluster_name}...")
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=current_cluster_name,
                neighbors_key=neighbors_key,
                flavor="igraph",
            )

            n_clusters = adata.obs[current_cluster_name].nunique()

            # handle palettes when there are > 25 clusters
            if n_clusters <= len(cluster_palette_25):
                adata.uns[f"{current_cluster_name}_colors"] = cluster_palette_25[
                    :n_clusters
                ]
            else:
                # Loop the palette so Scanpy doesn't crash from missing colors
                repeated_palette = cluster_palette_25 * (
                    (n_clusters // len(cluster_palette_25)) + 1
                )
                adata.uns[f"{current_cluster_name}_colors"] = repeated_palette[
                    :n_clusters
                ]

            # plot UMAP
            logger.info(f"Plotting UMAPs for {current_cluster_name}...")
            sc.pl.embedding(
                adata,
                basis=custom_umap_basis,  # Tell plot to use our uniquely saved UMAP
                color=[
                    "total_counts",
                    "n_genes_by_counts",
                    current_cluster_name,
                ],
                wspace=0.4,
                show=False,
                save=f"_{current_cluster_name}.png",
                frameon=False,
            )

            logger.info(f"Plotting Spatial Scatter for {current_cluster_name}...")
            fig, ax = plt.subplots(figsize=(15, 15), facecolor="white")
            ax.set_facecolor("white")

            sq.pl.spatial_scatter(
                adata,
                color=[current_cluster_name],
                spatial_key=spatial_key,
                shape=None,
                facecolor="white",
                size=2,
                frameon=False,
                img=False,
                ax=ax,
                outline=False,
                dpi=300,
            )

            fig.savefig(
                combo_dir / f"{current_cluster_name}_spatial.png",
                dpi=300,
                facecolor="white",
                bbox_inches="tight",
            )
            plt.close(fig)

    # Reset global figdir
    sc.settings.figdir = module_dir

    # Cleanup temporary pre-clustering columns used for scVIVA to keep adata clean
    if "pre_scviva_labels" in adata.obs:
        del adata.obs["pre_scviva_labels"]
    if "_scviva_dummy_sample" in adata.obs:
        del adata.obs["_scviva_dummy_sample"]

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
