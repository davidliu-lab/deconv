import logging

import pandas as pd

logger = logging.getLogger(__name__)


def get_genes_not_too_sparse(
    df_rnaseq: pd.DataFrame, sparsity_ceiling: float
) -> set[str]:
    return set(df_rnaseq[(df_rnaseq == 0).mean(axis=1) < sparsity_ceiling].index)


def get_good_genes(
    df_bulkrnaseq_tcga_skcm: pd.DataFrame,
    df_scrnaseq_jerby_arnon: pd.DataFrame,
    sparsity_ceiling_tcga: float,
) -> set[str]:
    genes_not_too_sparse_in_tcga_skcm = get_genes_not_too_sparse(
        df_bulkrnaseq_tcga_skcm, sparsity_ceiling_tcga
    )
    logger.debug(
        "genes not too sparse in tcga_skcm: %s", len(genes_not_too_sparse_in_tcga_skcm)
    )
    genes_in_both_datasets = set(
        df_bulkrnaseq_tcga_skcm.index.intersection(df_scrnaseq_jerby_arnon.index)
    )
    logger.debug("genes in both datasets: %s", len(genes_in_both_datasets))
    good_genes = genes_in_both_datasets & genes_not_too_sparse_in_tcga_skcm
    logger.debug("good genes: %s", len(good_genes))
    return good_genes
