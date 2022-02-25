from typing import Tuple
import pandas as pd

from helpers.cell_type_naming import weird_to_nice


CELL_TYPE_COLUMN_NAME = "cell_type"
GENE_SYMBOL_COLUMN_NAME = "gene_symbol"
SAMPLE_COLUMN_NAME = "sample_id"
SINGLE_CELL_COLUMN_NAME = "single_cell_id"


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
    sc_rna_seq.rename_axis(index=GENE_SYMBOL_COLUMN_NAME, columns=SINGLE_CELL_COLUMN_NAME, inplace=True)
    sc_rna_seq.sort_index(inplace=True)

    metadata = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_cell.annotations.csv",
        na_values={"cell.types": "?"},
    )
    metadata = metadata.replace({"cell.types": weird_to_nice})

    return sc_rna_seq, metadata


def load_tcga_skcm(n_genes: int = None) -> pd.DataFrame:
    """Load RNA-seq mixtures for the TCGA SKCM cohort, processed by Derek

    Returns:
        mixtures_tcga_skcm: pandas.DataFrame
    """
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/skcm_rnaseqv2_normalized_clean.txt"

    mixtures_tcga_skcm = pd.read_csv(path, sep="\t", index_col=0, nrows=n_genes)
    mixtures_tcga_skcm.rename_axis(index=GENE_SYMBOL_COLUMN_NAME, columns=SAMPLE_COLUMN_NAME, inplace=True)
    mixtures_tcga_skcm.sort_index(inplace=True)

    return mixtures_tcga_skcm


def load_tcga_skcm_fractions_from_csx() -> pd.DataFrame:
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/CIBERSORTx_Job8_Results.txt"
    fractions = pd.read_csv(path, sep="\t", index_col=0)
    fractions = fractions.drop(fractions.columns[-3:], axis=1)
    fractions = fractions.rename(columns=weird_to_nice)
    fractions = fractions.reindex(sorted(fractions.columns), axis=1)
    fractions = fractions.rename_axis(index=SAMPLE_COLUMN_NAME, columns=CELL_TYPE_COLUMN_NAME)
    return fractions
