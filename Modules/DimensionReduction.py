"""Dimension reduction module."""

import gc
import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi
import squidpy as sq
import torch
from scvi.external import SCVIVA

matplotlib.use("Agg")


warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def run_scanvi_transfer(
    adata_query_full,  # <--- The untouched 10GB original object
    spatial_key,
    ref_path,
    query_batch_key,
    scviva_labels_key,
    scviva_embedding_key,
    ref_batch_key=None,
    query_layer="counts",
    accelerator="cpu",
    devices="auto",
    scanvi_epochs=100,
):
    logger.info("--- Initiating scANVI Reference Mapping (scArches) ---")
    adata_query = adata_query_full.copy()

    # Move raw counts to .X for the query copy
    if query_layer and query_layer in adata_query.layers:
        adata_query.X = adata_query.layers[query_layer].copy()

    # 2. Fetch pre-trained model
    if "/" in ref_path and not Path(ref_path).exists():
        logger.info(f"Downloading/Loading pre-trained model '{ref_path}'...")
        try:
            from scvi.hub import HubModel

            hub_model = HubModel.pull_from_huggingface_hub(ref_path)
            model_path = hub_model.local_dir
        except ImportError:
            logger.error("huggingface_hub is not installed.")
            raise
    else:
        logger.info(f"Using local pre-trained model directory: {ref_path}")
        model_path = ref_path

    # 2.5 HARMONIZE GENE NAMES ON THE COPY: is there a better way to do this?
    logger.info("Loading reference adata to check gene names...")
    try:
        ref_adata = ad.read_h5ad(Path(model_path) / "adata.h5ad", backed="r")
        ref_var_names = ref_adata.var_names.values

        overlap = len(set(adata_query.var_names).intersection(ref_var_names))
        logger.info(f"Initial gene overlap: {overlap} / {len(adata_query.var_names)}")

        if overlap < 10:
            logger.info(
                "Overlap is ~0. Translating CosMx Gene Symbols to Ensembl IDs..."
            )
            symbol_col = None
            for col in ref_adata.var.columns:
                sample_genes = ref_adata.var[col].astype(str).values
                if (
                    "EPCAM" in sample_genes
                    or "KRT8" in sample_genes
                    or "CD8A" in sample_genes
                ):
                    symbol_col = col
                    break

            if symbol_col:
                sym_to_ref = dict(
                    zip(ref_adata.var[symbol_col].astype(str), ref_adata.var_names)
                )
                adata_query.var_names = [
                    sym_to_ref.get(g, g) for g in adata_query.var_names
                ]
                new_overlap = len(
                    set(adata_query.var_names).intersection(ref_var_names)
                )
                logger.info(
                    f"New overlap after translation: {new_overlap} / {len(adata_query.var_names)}"
                )
    except Exception as e:
        logger.warning(f"Could not perform gene translation: {e}")

    # 2.8 MAP QUERY BATCH TO REFERENCE BATCH KEY
    if ref_batch_key:
        logger.info(
            f"Mapping query batch key '{query_batch_key}' to '{ref_batch_key}'..."
        )
        adata_query.obs[ref_batch_key] = adata_query.obs[query_batch_key].astype(str)

    import json

    try:
        with open(Path(model_path) / "_scvi_required_metadata.json", "r") as f:
            meta = json.load(f)
            if "tissue_in_vivo" in str(meta):
                adata_query.obs["tissue_in_vivo"] = "large_intestine"
    except Exception:
        pass

    # 3. Prepare Query Data
    logger.info("Preparing query data to match reference model features...")
    scvi.model.SCANVI.prepare_query_anndata(adata_query, model_path)

    # 4. Load Query Model (scArches architectural surgery)
    logger.info("Loading query data into pre-trained scANVI model...")
    scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, model_path)

    # 5. Train only the query adapters
    logger.info("Training query mapping (scArches)...")
    scanvi_query.train(
        max_epochs=scanvi_epochs,
        plan_kwargs={"weight_decay": 0.0},
        accelerator=accelerator,
        devices=devices,
    )

    logger.info("Extracting scANVI predictions and mapping back to ORIGINAL object...")
    adata_query_full.obs[scviva_labels_key] = scanvi_query.predict()
    adata_query_full.obsm[scviva_embedding_key] = (
        scanvi_query.get_latent_representation()
    )

    del adata_query
    del scanvi_query

    gc.collect()

    logger.info("--- scANVI Label Transfer Complete ---")

    return adata_query_full


def run_scanvi_joint(
    adata_query_full,
    ref_path,
    query_batch_key,
    ref_batch_key,
    ref_label_key,
    scviva_labels_key,
    scviva_embedding_key,
    query_layer="counts",
    accelerator="cpu",
    devices="auto",
    scvi_epochs=400,
    scanvi_epochs=200,
):
    logger.info("--- Initiating De Novo Joint scANVI Integration ---")
    adata_query = adata_query_full.copy()

    # Move raw counts to .X for the query copy
    if query_layer and query_layer in adata_query.layers:
        adata_query.X = adata_query.layers[query_layer].copy()

    # 1. Load reference data
    logger.info(f"Loading reference adata from {ref_path}...")
    ref_adata = ad.read_h5ad(ref_path)

    # 2. Harmonize genes (Find shared genes)
    query_genes = set(adata_query.var_names)
    ref_genes = set(ref_adata.var_names)
    shared_genes = list(query_genes.intersection(ref_genes))

    # Optional case-insensitive check if overlap is surprisingly low
    if len(shared_genes) < len(query_genes) * 0.5:
        logger.info(
            "Initial gene overlap is low. Attempting case-insensitive matching..."
        )
        query_upper = {str(g).upper(): g for g in query_genes}
        ref_upper = {str(g).upper(): g for g in ref_genes}

        shared_upper = set(query_upper.keys()).intersection(set(ref_upper.keys()))

        # Rename query genes to match reference exactly
        rename_dict = {query_upper[g_up]: ref_upper[g_up] for g_up in shared_upper}
        adata_query.var_names = [rename_dict.get(g, g) for g in adata_query.var_names]
        shared_genes = list(set(adata_query.var_names).intersection(ref_genes))

    logger.info(f"Shared genes for joint training: {len(shared_genes)}")

    # Subset both to the shared genes
    adata_query = adata_query[:, shared_genes].copy()
    ref_adata = ref_adata[:, shared_genes].copy()

    # 3. Setup joint metadata
    SCANVI_CELLTYPE_KEY = "scanvi_label"
    adata_query.obs[SCANVI_CELLTYPE_KEY] = "Unknown"

    if ref_label_key not in ref_adata.obs.columns:
        raise ValueError(
            f"reference_label_key '{ref_label_key}' not found in reference."
        )
    ref_adata.obs[SCANVI_CELLTYPE_KEY] = ref_adata.obs[ref_label_key].astype(str)

    # Setup Joint Batch Key
    JOINT_BATCH_KEY = "_joint_batch"
    adata_query.obs[JOINT_BATCH_KEY] = adata_query.obs[query_batch_key].astype(str)

    if ref_batch_key and ref_batch_key in ref_adata.obs.columns:
        ref_adata.obs[JOINT_BATCH_KEY] = ref_adata.obs[ref_batch_key].astype(str)
    else:
        ref_adata.obs[JOINT_BATCH_KEY] = "reference_batch"

    # Distinguish datasets for later extraction
    adata_query.obs["_dataset"] = "query"
    ref_adata.obs["_dataset"] = "reference"

    # 4. Concatenate
    logger.info("Concatenating query and reference datasets...")
    adata_joint = ad.concat([adata_query, ref_adata])

    del adata_query
    del ref_adata

    gc.collect()
    logger.info("Cleared original separate datasets from memory to save RAM.")

    # 5. Train scVI from scratch on the joint object
    logger.info("Setting up scVI for joint object...")
    scvi.model.SCVI.setup_anndata(adata_joint, batch_key=JOINT_BATCH_KEY)
    scvi_model = scvi.model.SCVI(adata_joint, n_layers=2, n_latent=30)

    logger.info(f"Training scVI for {scvi_epochs} epochs...")
    scvi_model.train(
        max_epochs=scvi_epochs,
        accelerator=accelerator,
        devices=devices,
        check_val_every_n_epoch=10,
    )

    # 6. Train scANVI on top of the joint scVI model
    logger.info("Initializing scANVI from joint scVI model...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        labels_key=SCANVI_CELLTYPE_KEY,
    )
    logger.info(f"Training scANVI for {scanvi_epochs} epochs...")
    scanvi_model.train(
        max_epochs=scanvi_epochs,
        accelerator=accelerator,
        devices=devices,
        check_val_every_n_epoch=10,
    )

    # 7. Extract predictions and latent space
    logger.info("Extracting scANVI predictions and mapping back to ORIGINAL object...")
    adata_joint.obs[scviva_labels_key] = scanvi_model.predict(adata_joint)
    adata_joint.obsm[scviva_embedding_key] = scanvi_model.get_latent_representation(
        adata_joint
    )

    # 8. Subset back to just the query cells
    query_mask = adata_joint.obs["_dataset"] == "query"
    adata_query_results = adata_joint[query_mask].copy()

    # Map results back to the untouch original full AnnData
    adata_query_full.obs[scviva_labels_key] = adata_query_results.obs[
        scviva_labels_key
    ].values
    adata_query_full.obsm[scviva_embedding_key] = adata_query_results.obsm[
        scviva_embedding_key
    ]

    # Cleanup memory
    # del adata_query
    # del ref_adata
    del adata_joint
    gc.collect()

    logger.info("--- Joint scANVI Integration Complete ---")
    return adata_query_full


def run_SCVIVA(
    adata,
    scviva_batch_key,
    scviva_sample_key,
    spatial_key,
    scviva_spatial_knn,
    scviva_epochs,
    scviva_layer,
    use_scanvi=False,
    scanvi_mode="scarches",
    ref_path=None,
    ref_batch_key=None,
    ref_label_key=None,
    ref_layer="X",
    scvi_epochs=400,
    scanvi_epochs=200,
    scviva_batch_size=512,
    n_latent=30,
    n_layers=2,
    pre_cluster_res=1.0,
):
    """
    Executes the Spatially-Aware Variational Autoencoder (scVIVA) pipeline.

    This function initializes and trains a spatially-unaware baseline model (scVI)
    to extract batch-corrected embeddings, generates unsupervised biological labels
    using Leiden clustering, and subsequently sets up spatial niche graphs to train
    a spatially-aware scVIVA model. The resulting latent representations are stored
    back into the AnnData object.

    Args:
        adata (anndata.AnnData): The annotated data matrix containing spatial
            transcriptomics data.
        scviva_batch_key (str or None): Column name in `adata.obs` corresponding to
            batch_id / slide_id. If None or missing, a dummy key is automatically created.
        scviva_batch_key (str or None): Column name in `adata.obs` corresponding to
            sample_id. If None or missing, a dummy key is automatically created.
        spatial_key (str): Key in `adata.obsm` containing the spatial coordinates
            (e.g., 'spatial', 'global').
        scviva_spatial_knn (int): Number of nearest neighbors to use when constructing
            the spatial niche graph.
        scviva_epochs (int): Maximum number of training epochs for the scVIVA model.
        scviva_layer (str): Key in `adata.layers` where raw integer counts are stored.

    Returns:
        anndata.AnnData: The updated AnnData object containing the new latent
            representations ('X_scVI' and 'X_scVIVA') in `adata.obsm`.
    """
    logger.info("--- Initiating scVIVA Pipeline ---")

    use_gpu = torch.cuda.is_available()
    accelerator = "gpu" if use_gpu else "cpu"
    devices = 1 if use_gpu else "auto"

    if use_gpu:
        logger.info(
            f"GPU detected. Training will use GPU (accelerator='{accelerator}')."
        )
    else:
        logger.warning(
            f"No GPU detected! Falling back to CPU (accelerator='{accelerator}'). Training will be slower."
        )

    # Ensure a sample_key exists for scVIVA
    if scviva_batch_key and scviva_sample_key in adata.obs.columns:
        logger.info(
            f"Running scVIVA with batch key: {scviva_batch_key}, sample_key: {scviva_sample_key}"
        )
    else:
        scviva_sample_key = "_scviva_dummy_sample"
        adata.obs[scviva_sample_key] = "sample_1"
        logger.info(
            f"No valid scviva_batch_key scviva_sample_key or provided. Created dummy sample and batch key '{scviva_sample_key}' for spatial graph construction."
        )
        scviva_batch_key = scviva_sample_key

    if use_scanvi and ref_path:
        # run scanvi scArches route -> scviva
        scviva_labels_key = "C_scANVI"
        scviva_embedding_key = "X_scANVI"

        if scanvi_mode == "joint":
            # De Novo Joint Integration
            adata = run_scanvi_joint(
                adata_query_full=adata,
                ref_path=ref_path,
                query_batch_key=scviva_batch_key,
                ref_batch_key=ref_batch_key,
                ref_label_key=ref_label_key,
                scviva_labels_key=scviva_labels_key,
                scviva_embedding_key=scviva_embedding_key,
                query_layer=scviva_layer,
                accelerator=accelerator,
                devices=devices,
                scvi_epochs=scvi_epochs,
                scanvi_epochs=scanvi_epochs,
            )
        else:
            # scArches Label Transfer
            adata = run_scanvi_transfer(
                adata_query_full=adata,
                spatial_key=spatial_key,
                ref_path=ref_path,
                query_batch_key=scviva_batch_key,
                scviva_labels_key=scviva_labels_key,
                scviva_embedding_key=scviva_embedding_key,
                ref_batch_key=ref_batch_key,
                query_layer=scviva_layer,
                accelerator=accelerator,
                devices=devices,
                scanvi_epochs=scanvi_epochs,
            )

    else:
        # standard scvi -> scviva
        logger.info("Run Spatially-Unaware Baseline (scVI)")

        scvi.model.SCVI.setup_anndata(
            adata, layer=scviva_layer, batch_key=scviva_batch_key
        )
        scvi_model = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent)
        scvi_model.train(
            max_epochs=scvi_epochs, accelerator=accelerator, devices=devices
        )

        # Extract the clean, batch-corrected baseline embedding from scVI
        adata.obsm["X_scVI"] = scvi_model.get_latent_representation()

        logger.info("Generating unsupervised labels for scVIVA environment features...")
        sc.pp.neighbors(
            adata, use_rep="X_scVI", key_added="pre_scviva_neighbors", n_neighbors=15
        )
        sc.tl.leiden(
            adata,
            resolution=pre_cluster_res,
            key_added="pre_scviva_labels",
            neighbors_key="pre_scviva_neighbors",
        )

        scviva_labels_key = "pre_scviva_labels"
        scviva_embedding_key = "X_scVI"

    logger.info("Run Spatially-Aware scVIVA")

    setup_kwargs = {
        "labels_key": scviva_labels_key,
        "cell_coordinates_key": spatial_key,
        "expression_embedding_key": scviva_embedding_key,
        "sample_key": scviva_sample_key,  # sample key for scviva to ensure spatial correlation stays within samples
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
        batch_key=scviva_batch_key,
        **setup_kwargs,
    )

    logger.info("Training scVIVA model...")
    nichevae = SCVIVA(adata)
    nichevae.train(
        max_epochs=scviva_epochs,
        accelerator=accelerator,
        devices=devices,
        early_stopping=True,
        check_val_every_n_epoch=10,
        batch_size=scviva_batch_size,
        plan_kwargs={"lr": 5e-4},
    )

    logger.info("Extracting scVIVA latent representation...")
    adata.obsm["X_scVIVA"] = nichevae.get_latent_representation()
    logger.info("--- scVIVA Pipeline Complete ---")

    return adata


def run_dimension_reduction(
    data_type,
    input_adata_path,
    module_dir,
    module_name,
    n_neighbors_list,
    resolution_list,
    cluster_name,
    n_comps=None,
    use_scviva=True,
    umap_latent="X_scVI",
    scviva_layer="counts",
    scviva_batch_key=None,
    scviva_sample_key=None,
    scviva_spatial_knn=20,
    scviva_epochs=400,
    run_pca=False,
    use_scanvi=False,
    scanvi_mode="scarches",
    reference_adata_path=None,
    reference_batch_key=None,
    reference_label_key=None,
    reference_layer="X",
    scvi_epochs=400,
    scanvi_epochs=200,
    scviva_batch_size=512,
    n_latent=30,
    n_layers=2,
    pre_cluster_res=1.0,
    dot_size=0.5,
):
    """
    Runs the complete dimension reduction, clustering, and visualization workflow.

    This function orchestrates the loading of spatial transcriptomics data, optional
    PCA computation, execution of the scVIVA pipeline, and iterative Leiden clustering
    and UMAP embeddings across multiple neighbor and resolution parameters. It handles
    dynamically saving output figures and writes the final annotated data matrix to disk.

    Args:
        data_type (str): The platform/type of the spatial data (e.g., "CosMx", "Xenium").
            This determines which default spatial coordinate key to use.
        input_adata_path (str or pathlib.Path): File path to the input AnnData (.h5ad) file.
        module_dir (pathlib.Path): Output directory where generated figures and the
            updated AnnData file will be saved.
        module_name (str): Identifier name for the module.
        n_comps (int): Number of principal components to compute (used only if `run_pca` is True).
        n_neighbors_list (int or list of int): Number of neighbors (or a list of neighbor
            values) to use for computing the neighborhood graph.
        resolution_list (float or list of float): Resolution values (or a list of values)
            to use for Leiden clustering.
        cluster_name (str): Base prefix to use for the output cluster column names
            stored in `adata.obs`.
        scviva_layer (str, optional): Key in `adata.layers` containing raw integer counts,
            which is required by scVI. Defaults to "counts".
        scviva_batch_key (str, optional): Column in `adata.obs` indicating batch or
            sample labels. Defaults to None.
        scviva_spatial_knn (int, optional): Number of spatial neighbors for the scVIVA
            graph construction. Defaults to 20.
        scviva_epochs (int, optional): Maximum training epochs for the scVIVA model.
            Defaults to 400.
        run_pca (bool, optional): If True, computes and plots a PCA variance ratio diagram
            before running scVIVA. Defaults to False.

    Returns:
        None: The function does not return an object; instead, it saves the modified
            AnnData object directly to disk within `module_dir`.
    """

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
    input_adata_path = Path(input_adata_path)
    adata = sc.read_h5ad(input_adata_path)

    # Ensure layer exists for scvi
    if scviva_layer not in adata.layers:
        logger.warning(
            f"Layer '{scviva_layer}' not found. scVI needs raw counts. Using adata.X instead."
        )
        scviva_layer = None

    # --- UPDATED PCA COMPUTATION ---
    if run_pca or (umap_latent == "X_pca"):
        logger.info("Computing PCA...")

        # 1. Find HVGs using the raw integer counts layer
        sc.pp.highly_variable_genes(
            adata, layer="counts", flavor="seurat_v3", n_top_genes=2000
        )

        # 2. Subset the adata (which currently has the log-normalized .X from QC) to just HVGs
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()

        # 3. Scale the log-normalized HVGs (This Centers the data for PCA)
        sc.pp.scale(adata_hvg, max_value=10)

        # 4. Default to 50 if user didn't specify to calculate the curve
        compute_pcs = n_comps if n_comps else 50
        sc.tl.pca(adata_hvg, n_comps=compute_pcs)

        if not n_comps:
            logger.info(
                "No 'n_comps' specified. Running heuristic elbow detection to find optimal PCs..."
            )
            var_ratio = adata_hvg.uns["pca"]["variance_ratio"]
            n_points = len(var_ratio)

            p1 = np.array([0, var_ratio[0]])
            p2 = np.array([n_points - 1, var_ratio[-1]])

            max_dist = -1
            optimal_pc = 10  # Safe fallback

            for i in range(n_points):
                p3 = np.array([i, var_ratio[i]])
                dist = np.abs(np.cross(p2 - p1, p3 - p1)) / np.linalg.norm(p2 - p1)
                if dist > max_dist:
                    max_dist = dist
                    optimal_pc = i + 1

            logger.info(
                f"==> Optimal number of Principal Components determined: {optimal_pc}"
            )
            adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"][:, :optimal_pc].copy()
        else:
            adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

        sc.pl.pca_variance_ratio(
            adata_hvg, log=True, n_pcs=compute_pcs, show=False, save="PCA.png"
        )
        logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    # --- CONDITIONAL SCVIVA BYPASS ---
    adata_path = module_dir / input_adata_path.name

    if use_scviva:
        if not adata_path.exists():
            adata = run_SCVIVA(
                adata=adata,
                scviva_batch_key=scviva_batch_key,
                scviva_sample_key=scviva_sample_key,
                spatial_key=spatial_key,
                scviva_spatial_knn=scviva_spatial_knn,
                scviva_epochs=scviva_epochs,
                scviva_layer=scviva_layer,
                use_scanvi=use_scanvi,
                scanvi_mode=scanvi_mode,
                ref_path=reference_adata_path,
                ref_batch_key=reference_batch_key,
                ref_label_key=reference_label_key,
                ref_layer=reference_layer,
                scvi_epochs=scvi_epochs,
                scanvi_epochs=scanvi_epochs,
                scviva_batch_size=scviva_batch_size,
                n_latent=n_latent,
                n_layers=n_layers,
                pre_cluster_res=pre_cluster_res,
            )
        else:
            logger.info(
                f"Found existing scVIVA output at {adata_path}. Loading instead of re-running."
            )
            adata = sc.read_h5ad(adata_path)
    else:
        logger.info(
            "Bypassing scVI/scVIVA deep learning models (use_scviva=False). Standard Dimension Reduction active."
        )

    # --- DOWNSTREAM CLUSTERING ---
    for n_neighbors in n_neighbors_list:
        logger.info(f"Compute neighbors for n={n_neighbors} using {umap_latent}...")

        # 1. Create a unique key for this neighbor graph
        neighbors_key = f"neighbors_n{n_neighbors}_{umap_latent}"
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            use_rep=umap_latent,
            key_added=neighbors_key,
        )

        logger.info(f"Create UMAPs and cluster cells for n={n_neighbors}...")
        sc.tl.umap(adata, neighbors_key=neighbors_key, min_dist=0.1, spread=0.6)

        # Scanpy saves the UMAP to 'X_umap' by default.
        # We copy it to a unique name so the next loop doesn't overwrite it
        custom_umap_basis = f"umap_n{n_neighbors}_{umap_latent}"
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
                size=dot_size,
                save=f"_{current_cluster_name}.png",
                frameon=False,
            )

            logger.info(f"Plotting Spatial Scatter for {current_cluster_name}...")

            sq.pl.spatial_scatter(
                adata,
                color=[current_cluster_name],
                library_key=scviva_sample_key,  # split spatial plots by sample
                spatial_key=spatial_key,
                shape=None,
                facecolor="white",
                size=dot_size,
                frameon=False,
                img=False,
                outline=False,
                dpi=300,
                figsize=(15, 15),
            )

            plt.savefig(
                combo_dir / f"{current_cluster_name}_spatial.png",
                dpi=300,
                facecolor="white",
                bbox_inches="tight",
            )
            plt.close()

    # Reset global figdir
    sc.settings.figdir = module_dir

    # Cleanup temporary pre-clustering columns used for scVIVA to keep adata clean
    if "pre_scviva_labels" in adata.obs:
        del adata.obs["pre_scviva_labels"]
    if "_scviva_dummy_sample" in adata.obs:
        del adata.obs["_scviva_dummy_sample"]

    # Save anndata object
    out_path = module_dir / input_adata_path.name
    adata.write_h5ad(out_path)
    logger.info(f"Data saved to {out_path}")
