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
    genes_top_500 = gene_sparsity_in_malignant_cells.nsmallest(500)
    logger.debug("summary of genes_top_500: %s", genes_top_500.describe())
    genes_with_sparsity = genes_top_500.sample(100, replace=False, random_state=rng)
    logger.debug("shape of genes_with_sparsity: %s", genes_with_sparsity.shape)
    logger.debug("genes_with_sparsity.head(): %s", genes_with_sparsity.head())
    return genes_with_sparsity


def perturb_scrnaseq_gene_expression_by_scaling_factor(
    df_scrnaseq: pd.DataFrame, genes_to_perturb: pd.Series, scaling_factor: float
) -> pd.DataFrame:
    genes = genes_to_perturb.index
    df_scrnaseq_perturbed = df_scrnaseq.copy()
    logger.debug(
        "Mean expression before: %s",
        df_scrnaseq_perturbed.loc[genes].mean(axis=1),
    )
    assert df_scrnaseq_perturbed.loc[genes].shape[0] == len(genes_to_perturb)
    df_scrnaseq_perturbed.loc[genes] *= scaling_factor
    logger.debug(
        "Mean expression after: %s",
        df_scrnaseq_perturbed.loc[genes].mean(axis=1),
    )
    assert df_scrnaseq_perturbed.shape == df_scrnaseq.shape, df_scrnaseq_perturbed.shape
    return df_scrnaseq_perturbed


def perturb_scrnaseq_gene_expression_by_scaling_factor_in_cell_type(
    df_scrnaseq: pd.DataFrame,
    df_scrnaseq_metadata: pd.DataFrame,
    cell_type: str,
    genes_to_perturb: pd.Series,
    scaling_factor: float,
) -> pd.DataFrame:
    df = df_scrnaseq.copy()
    cells_of_desired_type = df_scrnaseq_metadata[
        df_scrnaseq_metadata[helpers.columns.CELL_TYPE] == cell_type
    ].index
    assert len(cells_of_desired_type) > 0
    df_cell_type = df[cells_of_desired_type]
    df_cell_type = perturb_scrnaseq_gene_expression_by_scaling_factor(
        df_cell_type, genes_to_perturb, scaling_factor
    )
    df[cells_of_desired_type] = df_cell_type
    assert df.shape == df_scrnaseq.shape, df.shape
    return df
