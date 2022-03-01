from typing import Tuple

import pandas as pd

from helpers import columns
from helpers.cell_type_naming import weird_to_nice


def load_jerby_arnon(n_genes: int = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data

    Args:
        n_genes (int, optional): Max number of genes to read. If None (default), read all genes.

    Returns:
        tuple: two dataframes, one each of sc data and metadata
    """
    sc_rna_seq = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_tpm.csv",
        index_col=0,
        nrows=n_genes,
    )
    sc_rna_seq = sc_rna_seq.rename_axis(index=columns.GENE_SYMBOL, columns=columns.SINGLE_CELL_ID)
    sc_rna_seq = sc_rna_seq.sort_index(axis="columns")
    sc_rna_seq = sc_rna_seq.sort_index(axis="rows")
    metadata = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_cell.annotations.csv",
        na_values={"cell.types": "?"},
    )
    metadata = metadata.rename(
        columns={"cells": columns.SINGLE_CELL_ID, "cell.types": columns.CELL_TYPE}
    )
    metadata = metadata.replace({columns.CELL_TYPE: weird_to_nice})
    metadata = metadata.rename_axis(index=columns.SINGLE_CELL_ID)
    metadata = metadata.set_index(columns.SINGLE_CELL_ID, drop=False)
    metadata = metadata.sort_index()
    return sc_rna_seq, metadata


def load_tcga_skcm(n_genes: int = None) -> pd.DataFrame:
    """Load RNA-seq mixtures for the TCGA SKCM cohort, processed by Derek

    Returns:
        mixtures_tcga_skcm: pandas.DataFrame
    """
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/skcm_rnaseqv2_normalized_clean.txt"
    bulk_rna_seq = pd.read_csv(path, sep="\t", index_col=0, nrows=n_genes)
    bulk_rna_seq = bulk_rna_seq.reindex(sorted(bulk_rna_seq.columns), axis=1)
    bulk_rna_seq = bulk_rna_seq.rename_axis(index=columns.GENE_SYMBOL, columns=columns.SAMPLE_ID)
    bulk_rna_seq = bulk_rna_seq.sort_index()
    return bulk_rna_seq


def load_tcga_skcm_fractions_from_csx() -> pd.DataFrame:
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/CIBERSORTx_Job8_Results.txt"
    fractions = pd.read_csv(path, sep="\t", index_col=0)
    fractions = fractions.drop(fractions.columns[-3:], axis=1)
    fractions = fractions.rename(columns=weird_to_nice)
    fractions = fractions.reindex(sorted(fractions.columns), axis=1)
    fractions = fractions.rename_axis(index=columns.SAMPLE_ID, columns=columns.CELL_TYPE)
    fractions = fractions.sort_index()
    return fractions
