import logging

import numpy as np
import pandas as pd

import helpers

logger = logging.getLogger(__name__)


def determine_genes_to_perturb(
    df_scrnaseq: pd.DataFrame, df_sc_metadata: pd.DataFrame, rng: np.random.Generator
) -> pd.Index:
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
    genes_to_perturb = genes_with_sparsity.index
    logger.debug("n genes_to_perturb: %s", len(genes_to_perturb))
    logger.debug("genes_to_perturb: %s", genes_to_perturb)
    return genes_to_perturb


def perturb_scrnaseq_gene_expression(
    df_scrnaseq: pd.DataFrame,
    df_scrnaseq_metadata: pd.DataFrame,
    cell_type: str,
    genes: pd.Index,
    scaling_factor: float,
) -> pd.DataFrame:
    df = df_scrnaseq.copy()
    cells = df_scrnaseq_metadata[
        df_scrnaseq_metadata[helpers.columns.CELL_TYPE] == cell_type
    ].index
    assert len(cells) > 0
    df.loc[genes, cells] *= scaling_factor
    return df
