import os
from logging import getLogger

import celltypist
import scanpy as sc
from celltypist import models

logger = getLogger(__name__)


def run_CellTypist(
    adata,
    model_name,
    cluster_col,
    custom_model_path=None,
    train_model=False,
    train_data_path=None,
    train_labels_col=None,
):
    # Load adata
    adata_ct = adata.copy()

    # restore raw counts and normalize for predictions
    adata_ct.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    # training model
    if train_model:
        # Prevent re-training if it's already saved (e.g., looping through multiple resolutions)
        if custom_model_path and os.path.exists(custom_model_path):
            logger.info(
                f"Trained model already found at {custom_model_path}. Loading it instead of retraining..."
            )
            model = models.Model.load(model=custom_model_path)
        else:
            if not train_data_path or not train_labels_col:
                raise ValueError(
                    "To train a model, 'CellTypist_train_data' and 'CellTypist_train_labels' must be provided in config."
                )

            logger.info(
                f"Loading reference data from {train_data_path} to train custom CellTypist model..."
            )
            adata_ref = sc.read_h5ad(train_data_path)

            # Check if reference data needs normalization (heuristic: max value > 100 usually means raw counts)
            if adata_ref.X.max() > 100:
                logger.info("Normalizing and log1p transforming reference data...")
                sc.pp.normalize_total(adata_ref, target_sum=1e4)
                sc.pp.log1p(adata_ref)

            logger.info(
                f"Training CellTypist model on '{train_labels_col}' (This may take a few minutes)..."
            )
            model = celltypist.train(
                adata_ref, labels=train_labels_col, n_jobs=8, feature_selection=True
            )

            if custom_model_path:
                logger.info(f"Saving newly trained model to {custom_model_path}...")
                os.makedirs(os.path.dirname(custom_model_path), exist_ok=True)
                model.write(custom_model_path)

    elif custom_model_path:
        # Load directly from the provided local file path
        logger.info(f"Loading custom CellTypist model from {custom_model_path}...")
        model = models.Model.load(model=custom_model_path)
    else:
        # Default behavior: download/load from CellTypist registry
        logger.info(f"Loading default CellTypist model: {model_name}...")
        models.download_models(force_update=False)
        model = models.Model.load(model=f"{model_name}.pkl")

    # run prediction
    logger.info("Running CellTypist predictions...")
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,  # This enables the cluster-level consensus
        over_clustering=cluster_col,  # Point this to your cluster column
    )

    # convert predictions to adata
    adata_pred = predictions.to_adata()

    # Cell-level predictions (only save once)
    indivcellanno_col = "CellTypist_predictedlabels"
    confidence_col = "CellTypist_confidence"

    if indivcellanno_col not in adata.obs.columns:
        adata.obs[indivcellanno_col] = adata_pred.obs["predicted_labels"]

    if confidence_col not in adata.obs.columns:
        adata.obs[confidence_col] = adata_pred.obs["conf_score"]

    # Probability matrix (only save once)
    prob_key = "CellTypist_probabilities"
    if prob_key not in adata.obsm.keys():
        adata.obsm[prob_key] = predictions.probability_matrix.loc[
            adata.obs_names
        ].to_numpy()
        adata.uns["CellTypist_probability_columns"] = (
            predictions.probability_matrix.columns.tolist()
        )

    # Majority voting depends on clustering, so always save
    majorvotingcellanno_col = f"CellTypist_majorityvoting_{cluster_col}"
    adata.obs[majorvotingcellanno_col] = adata_pred.obs["majority_voting"]

    return adata, indivcellanno_col, majorvotingcellanno_col
