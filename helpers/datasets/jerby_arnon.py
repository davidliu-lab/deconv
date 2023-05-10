import logging

import pandas as pd

from helpers import columns
from helpers.cell_type_naming import weird_to_nice
from helpers.data_io_and_formatting.qa_gene_filtering import get_good_genes
from helpers.datasets.tcga_skcm import load_tcga_skcm_hg19_scaled_estimate_firebrowse

logger = logging.getLogger(__name__)


def load_jerby_arnon(ref_genome="hg19", units="tpm") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data (hg19 tpm) from GEO (GSE115978)

    :return: sc_data, sc_metadata
    """
    if ref_genome == "hg19" and units == "tpm":
        return load_jerby_arnon_hg19_tpm()
    else:
        raise NotImplementedError


def load_jerby_arnon_hg19_tpm(subset=False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data (hg19 tpm) from GEO (GSE115978)

    :return: sc_hg19_tpm, metadata
    """
    logger.debug("loading Jerby-Arnon scRNA-seq data")
    if subset:
        kwargs = {"skiprows": lambda i: i % 25}
    else:
        kwargs = {"engine": "pyarrow"}
    sc_hg19_tpm = (
        pd.read_csv(
            "gs://liulab/ftp/GSE115978/GSE115978_tpm.csv",
            index_col=0,
            **kwargs,
        )
        .rename_axis(index=columns.GENE_SYMBOL, columns=columns.SINGLE_CELL_ID)
        .sort_index(axis="columns")
        .sort_index(axis="rows")
    )
    logger.debug("loading Jerby-Arnon metadata")
    metadata = pd.read_csv(
        "gs://liulab/ftp/GSE115978/GSE115978_cell.annotations.csv",
        na_values={"cell.types": "?"},
    )
    metadata = metadata.rename(
        columns={
            "cells": columns.SINGLE_CELL_ID,
            "cell.types": columns.CELL_TYPE,
            "samples": columns.SAMPLE_ID,
        }
    )
    metadata = metadata.replace({columns.CELL_TYPE: weird_to_nice})
    metadata = metadata.rename_axis(index=columns.SINGLE_CELL_ID)
    metadata = metadata.set_index(columns.SINGLE_CELL_ID, drop=False)
    metadata = metadata.sort_index()
    return sc_hg19_tpm, metadata


def load_scrnaseq_and_filter_genes() -> tuple[pd.DataFrame, pd.DataFrame]:
    logger.debug("loading scRNA-seq data")
    df_scrnaseq, df_sc_metadata = load_jerby_arnon(ref_genome="hg19", units="tpm")
    logger.debug("loading bulk RNA-seq data for filtering")
    df_bulkrnaseq_tcga_skcm = load_tcga_skcm_hg19_scaled_estimate_firebrowse()
    logger.debug("determining good genes")
    good_genes = get_good_genes(df_bulkrnaseq_tcga_skcm, df_scrnaseq, 0.5)
    logger.debug("limiting scRNA-seq data to good genes")
    df_scrnaseq = df_scrnaseq.loc[list(sorted(good_genes))]
    logger.debug("Shape of scRNA-seq data after filtering: %s", df_scrnaseq.shape)
    return df_scrnaseq, df_sc_metadata
