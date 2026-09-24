"""CellTypist Annotation Module."""

import gc
import logging
import os
from contextlib import contextmanager
from typing import Tuple, Optional

import celltypist
import scanpy as sc
from anndata import AnnData
from celltypist import models

logger = logging.getLogger(__name__)


@contextmanager
def temporary_X(adata: AnnData, layer: str = "counts"):
    """
    Context manager to swap .X in-place to save memory, ensuring safe cleanup.
    Prevents deep copying the entire AnnData object.
    """
    logger.debug(f"Temporarily swapping .X with .layers['{layer}']")
    original_X = adata.X.copy()  # Only copies the sparse matrix, NOT the whole AnnData
    try:
        adata.X = adata.layers[layer].copy()
        yield adata
    finally:
        logger.debug("Restoring original .X")
        adata.X = original_X
        del original_X
        gc.collect()


def run_CellTypist(
    adata: AnnData,
    model_name: str,
    cluster_col: str,
    custom_model_path: Optional[str] = None,
    train_model: bool = False,
    train_data_path: Optional[str] = None,
    train_labels_col: Optional[str] = None,
) -> Tuple[AnnData, str, str]:
    """Runs CellTypist annotation with optimal memory management."""

    # 1. Model Loading / Training
    if train_model:
        if custom_model_path and os.path.exists(custom_model_path):
            logger.info(f"Trained model found at {custom_model_path}. Loading...")
            model = models.Model.load(model=custom_model_path)
        else:
            if not train_data_path or not train_labels_col:
                raise ValueError(
                    "Training requires 'train_data_path' and 'train_labels_col'."
                )

            logger.info(f"Loading reference data from {train_data_path} to train...")
            adata_ref = sc.read_h5ad(train_data_path)

            if adata_ref.X.max() > 100:
                logger.info("Normalizing and log1p transforming reference data...")
                sc.pp.normalize_total(adata_ref, target_sum=1e4)
                sc.pp.log1p(adata_ref)

            logger.info(f"Training CellTypist on '{train_labels_col}'...")
            model = celltypist.train(
                adata_ref, labels=train_labels_col, n_jobs=8, feature_selection=True
            )
            del adata_ref
            gc.collect()

            if custom_model_path:
                logger.info(f"Saving newly trained model to {custom_model_path}...")
                os.makedirs(os.path.dirname(custom_model_path), exist_ok=True)
                model.write(custom_model_path)

    elif custom_model_path:
        logger.info(f"Loading custom CellTypist model from {custom_model_path}...")
        model = models.Model.load(model=custom_model_path)
    else:
        logger.info(f"Loading default CellTypist model: {model_name}...")
        models.download_models(force_update=False)
        model = models.Model.load(model=f"{model_name}.pkl")

    # 2. Prediction (Memory-Safe)
    logger.info("Normalizing counts and running predictions...")
    with temporary_X(adata, layer="counts") as temp_adata:
        sc.pp.normalize_total(temp_adata, target_sum=1e4)
        sc.pp.log1p(temp_adata)

        predictions = celltypist.annotate(
            temp_adata,
            model=model,
            majority_voting=True,
            over_clustering=cluster_col,
        )

    # 3. Map predictions back safely
    adata_pred = predictions.to_adata()

    indivcellanno_col = "CellTypist_predictedlabels"
    confidence_col = "CellTypist_confidence"
    majorvotingcellanno_col = f"CellTypist_majorityvoting_{cluster_col}"

    if indivcellanno_col not in adata.obs.columns:
        adata.obs[indivcellanno_col] = adata_pred.obs["predicted_labels"]
    if confidence_col not in adata.obs.columns:
        adata.obs[confidence_col] = adata_pred.obs["conf_score"]

    adata.obs[majorvotingcellanno_col] = adata_pred.obs["majority_voting"]

    prob_key = "CellTypist_probabilities"
    if prob_key not in adata.obsm:
        adata.obsm[prob_key] = predictions.probability_matrix.loc[
            adata.obs_names
        ].to_numpy()
        adata.uns["CellTypist_probability_columns"] = (
            predictions.probability_matrix.columns.tolist()
        )

    return adata, indivcellanno_col, majorvotingcellanno_col
