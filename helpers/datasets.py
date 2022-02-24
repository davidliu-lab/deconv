from typing import Tuple
import pandas as pd

from helpers.cell_type_naming import weird_to_nice


GENE_SYMBOL_COLUMN_NAME = "gene_symbol"
SAMPLE_COLUMN_NAME = "sample_id"
SINGLE_CELL_COLUMN_NAME = "single_cell_id"


def load_jerby_arnon(n_genes: int = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon single cell data

    Args:
        n_genes_if_not_all (int, optional): Max number of genes to read. Defaults to None.

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
