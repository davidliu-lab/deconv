import pandas as pd

from helpers.cell_type_naming import weird_to_nice


def load_tirosh():
    sc_rna_seq = pd.read_csv("gs://liulab/ftp/GSE115978/GSE115978_tpm.csv", index_col=0)
    sc_rna_seq = sc_rna_seq.rename_axis(index="GeneSymbol", columns="cells")

    metadata = pd.read_csv("gs://liulab/ftp/GSE115978/GSE115978_cell.annotations.csv")
    metadata = metadata.replace({"cell.types": weird_to_nice})

    return sc_rna_seq, metadata
