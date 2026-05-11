import scanpy as sc
import celltypist
from celltypist import models


def run_CellTypist(adata, model_name, cluster_col):
    # Load adata
    adata_ct = adata.copy()

    # restore raw counts
    adata_ct.X = adata.layers["counts"].copy()

    # normalize to 10k
    sc.pp.normalize_total(adata_ct, target_sum=1e4)

    # log transform
    sc.pp.log1p(adata_ct)

    # Enabling `force_update = True` will overwrite existing (old) models.
    models.download_models(force_update=False)

    # Choose the model you want to employ
    model = models.Model.load(model=f"{model_name}.pkl")

    # run prediction
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,  # This enables the cluster-level consensus
        over_clustering=cluster_col,  # Point this to your cluster column
    )

    # add predictions to adata
    adata_pred = predictions.to_adata()

    # save into custom column
    indivcellanno_col = f"CellTypist_predictedlabels_{cluster_col}"
    majorvotingcellanno_col = f"CellTypist_majorityvoting_{cluster_col}"

    adata.obs[indivcellanno_col] = adata_pred.obs["predicted_labels"]
    adata.obs[majorvotingcellanno_col] = adata_pred.obs["majority_voting"]

    return adata, indivcellanno_col, majorvotingcellanno_col
