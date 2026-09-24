"""DE Analysis for specific conditions module using PyDESeq2."""

import re
import warnings
from logging import getLogger
from pathlib import Path
from typing import List, Tuple

import numpy as np
import scanpy as sc
from anndata import AnnData

# Import PyDESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def _run_pydeseq2_comparison(
    adata_sub: AnnData,
    treatment_col: str,
    test_group: str,
    ref_group: str,
    output_dir: Path,
) -> None:
    """Runs PyDESeq2 for a specific pairwise comparison and exports results."""
    comparison_name = f"{test_group}_vs_{ref_group}"

    sc.pp.filter_genes(adata_sub, min_counts=1)
    if adata_sub.n_obs < 2:
        logger.info("      Not enough samples for this comparison.")
        return

    try:
        # Run PyDESeq2
        dds = DeseqDataSet(
            adata=adata_sub, design_factors=treatment_col, refit_cooks=True, n_cpus=4
        )
        dds.deseq2()

        stat_res = DeseqStats(
            dds, contrast=[treatment_col, test_group, ref_group], n_cpus=4
        )
        stat_res.summary()

        # Format Results
        res_df = stat_res.results_df
        markers = res_df.reset_index().rename(
            columns={
                "index": "names",
                "log2FoldChange": "logfoldchanges",
                "pvalue": "pvals",
                "padj": "pvals_adj",
            }
        )

        # Handle NaNs from PyDESeq2
        markers["pvals_adj"] = markers["pvals_adj"].fillna(1.0)
        markers["pvals"] = markers["pvals"].fillna(1.0)
        markers["logfoldchanges"] = markers["logfoldchanges"].fillna(0.0)
        markers = markers.sort_values(by="logfoldchanges", ascending=False)

        # Save outputs
        markers.to_csv(output_dir / f"{comparison_name}_all_genes.csv", index=False)

        sig_markers = markers[
            (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"].abs() > 0.5)
        ]
        sig_markers.to_csv(
            output_dir / f"{comparison_name}_SIGNIFICANT_only.csv", index=False
        )

        logger.info(f"      Found {len(sig_markers)} significant DE genes.")

    except Exception as e:
        logger.error(f"      PyDESeq2 Error: {e}")


def targeted_pairwise_DE(
    pseudobulk_adata_path: str | Path,
    celltype_col: str,
    treatment_col: str,
    comparisons: List[Tuple[str, str]],
    module_dir: str | Path,
) -> None:
    """Performs specific targeted pairwise DE comparisons within each cell type using PyDESeq2."""
    adata = sc.read_h5ad(pseudobulk_adata_path)
    module_dir = Path(module_dir)

    # Ensure raw integer counts for PyDESeq2
    adata.X = np.round(adata.X).astype(int)
    cell_types = adata.obs[celltype_col].dropna().unique()

    for ct in cell_types:
        logger.info(f"\n{'=' * 60}\nRunning PyDESeq2 for cell type: {ct}\n{'=' * 60}")

        adata_ct = adata[adata.obs[celltype_col] == ct].copy()
        safe_ct = re.sub(r"[^\w\s-]", "", str(ct)).replace(" ", "_")
        output_dir = module_dir / safe_ct
        output_dir.mkdir(parents=True, exist_ok=True)

        for test_group, ref_group in comparisons:
            groups_present = adata_ct.obs[treatment_col].unique()
            if test_group not in groups_present or ref_group not in groups_present:
                logger.info(f"  [Skip] {test_group} or {ref_group} missing in {ct}")
                continue

            logger.info(
                f"\n  --> Comparing: {test_group} (Test) vs {ref_group} (Reference)"
            )

            adata_sub = adata_ct[
                adata_ct.obs[treatment_col].isin([test_group, ref_group])
            ].copy()
            _run_pydeseq2_comparison(
                adata_sub, treatment_col, test_group, ref_group, output_dir
            )

    logger.info("\nAll cell types and targeted comparisons processed successfully!")
