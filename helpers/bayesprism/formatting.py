import pandas as pd

import helpers


def validate_sc_rnaseq(df):
    assert df.index.name == "gene_symbol"
    assert df.columns.name == "single_cell_id"


def validate_sc_rnaseq_metadata(df):
    assert df.index.name == helpers.columns.SINGLE_CELL_ID
    assert helpers.columns.CELL_TYPE in df.columns


def validate_bulk_rnaseq(df):
    assert df.index.name == helpers.columns.GENE_SYMBOL
    assert df.columns.name == helpers.columns.SAMPLE_ID


def format_sc_rnaseq(df: pd.DataFrame, df_annotations: pd.DataFrame) -> pd.DataFrame:
    validate_sc_rnaseq(df)
    validate_sc_rnaseq_metadata(df_annotations)
    # returns df with columns for genes and rows for single cells
    return df.T


def format_sc_cell_types(df_annotations: pd.DataFrame) -> pd.DataFrame:
    # assumes df_annotations has rows for single cells and columns for metadata
    # returns df with rows for single cells and one column with cell type
    raise NotImplemented


def format_sc_cell_states(df_annotations: pd.DataFrame) -> pd.DataFrame:
    # assumes df_annotations has rows for single cells and columns for metadata
    # returns df with rows for single cells and one column with cell state
    raise NotImplemented


def format_bulk_rnaseq(df: pd.DataFrame) -> pd.DataFrame:
    # assumes df has columns for bulk samples and rows for genes
    assert df.index.name == helpers.columns.GENE_SYMBOL
    assert df.columns.name == helpers.columns.SAMPLE_ID
    # returns df with columns for genes and rows for bulk samples
    return df.T
