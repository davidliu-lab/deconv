import logging
import pathlib
from typing import Union

import pandas as pd
import upath

import helpers

logger = logging.getLogger(__name__)


def create_csx_mixtures_tsv(
    df_bulk_rnaseq: pd.DataFrame, path_target: Union[pathlib.Path, upath.UPath]
):
    df_bulk_rnaseq = df_bulk_rnaseq.rename_axis(index="GeneSymbol", columns=None)
    df_bulk_rnaseq = df_bulk_rnaseq.reset_index()
    df_bulk_rnaseq.to_csv(path_target, index=False, sep="\t")


def create_csx_refsample_tsv(
    df_refsample_sc_rnaseq: pd.DataFrame,
    df_refsample_sc_metadata: pd.DataFrame,
    path_target: Union[pathlib.Path, upath.UPath],
):
    def look_up_cell_type(single_cell_id):
        return df_refsample_sc_metadata.loc[single_cell_id][helpers.columns.CELL_TYPE]

    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename(columns=look_up_cell_type)
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename_axis(
        index="GeneSymbol", columns=None
    )
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.reset_index()
    df_refsample_sc_rnaseq.to_csv(path_target, index=False, sep="\t")


def create_csx_fractions_tsv(
    df_fractions: pd.DataFrame, path_target: Union[pathlib.Path, upath.UPath]
):
    df_fractions = df_fractions.rename_axis(index="Mixture", columns=None)
    df_fractions = df_fractions.reset_index()
    df_fractions.to_csv(path_target, index=False, sep="\t")


def create_refsample_from_jerby_arnon(path_target: Union[pathlib.Path, upath.UPath]):
    logger.debug("loading sc data")
    df_sc_rnaseq, df_metadata = helpers.datasets.load_jerby_arnon_hg19_tpm()
    df_sc_rnaseq *= 1e6 / df_sc_rnaseq.sum()
    logger.debug("shape of scRNA-seq: %s", df_sc_rnaseq.shape)
    nonnull_cell_type_cells = df_metadata[
        ~df_metadata[helpers.columns.CELL_TYPE].isna()
    ].index
    logger.debug("length of nonnull_cell_type_cells: %s", len(nonnull_cell_type_cells))
    logger.debug("%s", nonnull_cell_type_cells)
    df_sc_rnaseq = df_sc_rnaseq[nonnull_cell_type_cells]
    logger.debug("new shape of scRNA-seq: %s", df_sc_rnaseq.shape)
    logger.debug("creating and writing refsample tsv")
    create_csx_refsample_tsv(df_sc_rnaseq, df_metadata, path_target)


def create_csx_mixtures_for_tcga_skcm(path_target: Union[pathlib.Path, upath.UPath]):
    logger.debug("loading bulk rna-seq for tcga skcm")
    df_bulk_rnaseq = pd.read_parquet(
        "gs://liulab/data/pseudobulk_optimization/3_with_tcga_qc/mixtures_real_tcga_skcm/tpm.parquet"
    )
    df_bulk_rnaseq = df_bulk_rnaseq.pivot("gene_symbol", "aliquot_barcode", "tpm")
    logger.debug("creating tsv of bulk rna-seq for tcga skcm")
    create_csx_mixtures_tsv(
        df_bulk_rnaseq,
        path_target,
    )


def create_csx_mixtures_for_pseudobulk(path_target: Union[pathlib.Path, upath.UPath]):
    logger.debug("loading pseudobulk rna-seq")
    df_bulk_rnaseq = pd.read_parquet(
        "gs://liulab/data/pseudobulk_optimization/3_with_tcga_qc/mixtures/n_cells=5/malignant_from_one_sample=True/data.parquet"
    )
    df_bulk_rnaseq = df_bulk_rnaseq.pivot(
        "gene_symbol", "tcga_aliquot_barcode_for_fractions", "tpm"
    )
    logger.debug("creating tsv of pseudobulk rna-seq")
    create_csx_mixtures_tsv(
        df_bulk_rnaseq,
        path_target,
    )


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.datasets").setLevel("DEBUG")

    path_refsample_jerby_arnon = upath.UPath(
        "gs://liulab/data/pseudobulk_evaluation/csx_input_files/refsample_jerby_arnon_hg19_tpm.tsv"
    )
    create_refsample_from_jerby_arnon(path_refsample_jerby_arnon)

    logger.debug("creating csx mixtures for tcga skcm")
    create_csx_mixtures_for_tcga_skcm(
        upath.UPath(
            "gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_tcga_skcm.tsv"
        )
    )

    logger.debug("creating csx mixtures for pseudobulk")
    create_csx_mixtures_for_pseudobulk(
        upath.UPath(
            "gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_pseudobulk.tsv"
        )
    )
