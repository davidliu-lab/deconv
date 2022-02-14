import pandas as pd

from helpers.cell_type_naming import weird_to_nice


def load_jerby_arnon(n_genes_if_not_all: int = None):
    """Load Jerby-Arnon single cell data

    Args:
        n_genes_if_not_all (int, optional): Max number of genes to read. Defaults to None.

    Returns:
        tuple: two dataframes, one each of sc data and metadata
    """
    sc_rna_seq = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_tpm.csv",
        index_col=0,
        nrows=n_genes_if_not_all,
    )
    sc_rna_seq.rename_axis(index="GeneSymbol", columns="cells", inplace=True)
    sc_rna_seq.sort_index(inplace=True)

    metadata = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_cell.annotations.csv",
        na_values={"cell.types": "?"},
    )
    metadata = metadata.replace({"cell.types": weird_to_nice})

    return sc_rna_seq, metadata
