"""DE Analysis for specific conditions module using PyDESeq2."""

import re
import warnings
from logging import getLogger
from pathlib import Path

import numpy as np
import scanpy as sc

# Import PyDESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def targeted_pairwise_DE(
    pseudobulk_adata_path, celltype_col, treatment_col, comparisons, module_dir
):
    """
    Performs specific targeted pairwise DE comparisons within each cell type using PyDESeq2.
    """
    adata = sc.read_h5ad(pseudobulk_adata_path)
    module_dir = Path(module_dir)

    # Ensure raw integer counts for PyDESeq2
    adata.X = np.round(adata.X).astype(int)

    cell_types = adata.obs[celltype_col].dropna().unique()

    for ct in cell_types:
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Running PyDESeq2 for cell type: {ct}")
        logger.info(f"{'=' * 60}")

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
            comparison_name = f"{test_group}_vs_{ref_group}"

            # 1. Subset adata to ONLY the two groups we are comparing
            adata_sub = adata_ct[
                adata_ct.obs[treatment_col].isin([test_group, ref_group])
            ].copy()

            # Filter out genes that have 0 counts across all samples in this subset
            sc.pp.filter_genes(adata_sub, min_counts=1)

            if adata_sub.n_obs < 2:
                logger.info("      Not enough samples for this comparison.")
                continue

            try:
                # 2. Run PyDESeq2
                dds = DeseqDataSet(
                    adata=adata_sub,
                    design_factors=treatment_col,
                    refit_cooks=True,
                    n_cpus=4,
                )
                dds.deseq2()

                # Setup statistical test (Contrast: test_group vs ref_group)
                stat_res = DeseqStats(
                    dds, contrast=[treatment_col, test_group, ref_group], n_cpus=4
                )
                stat_res.summary()

                # 3. Extract and Format Results
                res_df = stat_res.results_df

                # Rename PyDESeq2 columns to match Scanpy format for WebVisPrep & UI
                markers = res_df.reset_index().rename(
                    columns={
                        "index": "names",
                        "log2FoldChange": "logfoldchanges",
                        "pvalue": "pvals",
                        "padj": "pvals_adj",
                    }
                )

                # Handle NaNs from PyDESeq2 (genes filtered due to low counts)
                markers["pvals_adj"] = markers["pvals_adj"].fillna(1.0)
                markers["pvals"] = markers["pvals"].fillna(1.0)
                markers["logfoldchanges"] = markers["logfoldchanges"].fillna(0.0)

                # Sort by logfoldchange
                markers = markers.sort_values(by="logfoldchanges", ascending=False)

                # Save ALL genes
                csv_filename = output_dir / f"{comparison_name}_all_genes.csv"
                markers.to_csv(csv_filename, index=False)

                # Save SIGNIFICANT genes (Now using real adjusted p-values!)
                sig_markers = markers[
                    (markers["pvals_adj"] < 0.05)
                    & (markers["logfoldchanges"].abs() > 0.5)
                ]
                sig_filename = output_dir / f"{comparison_name}_SIGNIFICANT_only.csv"
                sig_markers.to_csv(sig_filename, index=False)

                logger.info(f"      Found {len(sig_markers)} significant DE genes.")

            except Exception as e:
                logger.error(f"      PyDESeq2 Error: {e}")
                continue

    logger.info("\nAll cell types and targeted comparisons processed successfully!")
