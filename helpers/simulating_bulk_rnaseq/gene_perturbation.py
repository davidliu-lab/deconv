import logging

import numpy as np
import pandas as pd

import helpers

logger = logging.getLogger(__name__)


def determine_genes_to_perturb(
    df_scrnaseq: pd.DataFrame, df_sc_metadata: pd.DataFrame, rng: np.random.Generator
) -> pd.Series:
    is_malignant_cell = df_sc_metadata[helpers.columns.CELL_TYPE] == "Malignant"
    malignant_single_cell_ids = df_sc_metadata[is_malignant_cell][
        helpers.columns.SINGLE_CELL_ID
    ]
    malignant_scrnaseq = df_scrnaseq[malignant_single_cell_ids]
    logger.debug("shape of malignant_scrnaseq: %s", malignant_scrnaseq.shape)
    gene_sparsity_in_malignant_cells = (malignant_scrnaseq == 0).mean(axis="columns")
    genes_top_500 = gene_sparsity_in_malignant_cells.sort_values().iloc[-500:]
    genes_with_sparsity = genes_top_500.sample(100, replace=False, random_state=rng)
    logger.debug("shape of genes_with_sparsity: %s", genes_with_sparsity.shape)
    logger.debug("genes_with_sparsity.head(): %s", genes_with_sparsity.head())
    return genes_with_sparsity


def perturb_scrnaseq_gene_2x(
    df_scrnaseq: pd.DataFrame, genes_to_perturb: pd.Series
) -> pd.DataFrame:
    genes = genes_to_perturb.index
    df_scrnaseq_perturbed = df_scrnaseq.copy()
    logger.debug(
        "Mean expression before: %s",
        df_scrnaseq_perturbed.loc[genes].mean(axis=1),
    )
    assert df_scrnaseq_perturbed.loc[genes].shape[0] == len(genes_to_perturb)
    df_scrnaseq_perturbed.loc[genes] *= 2.0
    logger.debug(
        "Mean expression after: %s",
        df_scrnaseq_perturbed.loc[genes].mean(axis=1),
    )
    assert df_scrnaseq_perturbed.shape == df_scrnaseq.shape, df_scrnaseq_perturbed.shape
    return df_scrnaseq_perturbed
