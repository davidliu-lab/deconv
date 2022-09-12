import logging
from typing import Tuple

import pandas as pd

from helpers import columns
from helpers.cell_type_naming import weird_to_nice

logger = logging.getLogger(__name__)


def load_jerby_arnon(
    ref_genome="hg19", units="tpm"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data (hg19 tpm) from GEO (GSE115978)

    :return: sc_data, sc_metadata
    """
    if ref_genome == "hg19" and units == "tpm":
        return load_jerby_arnon_hg19_tpm()
    else:
        raise NotImplementedError


def load_jerby_arnon_hg19_tpm() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data (hg19 tpm) from GEO (GSE115978)

    :return: sc_hg19_tpm, metadata
    """
    logger.debug("loading Jerby-Arnon scRNA-seq data")
    sc_hg19_tpm = (
        pd.read_csv(
            "gs://liulab/ftp/GSE115978/GSE115978_tpm.csv",
            index_col=0,
            engine="pyarrow",
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
