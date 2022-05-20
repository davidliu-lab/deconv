import logging
from typing import Tuple

from google.cloud import bigquery
import pandas as pd

from helpers import columns
from helpers.cell_type_naming import weird_to_nice

logger = logging.getLogger(__name__)


def load_jerby_arnon(ref_genome="hg19", units="tpm") -> Tuple[pd.DataFrame, pd.DataFrame]:
    if ref_genome == "hg19" and units == "tpm":
        return load_jerby_arnon_hg19_tpm()
    else:
        raise NotImplementedError


def load_jerby_arnon_hg19_tpm() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Jerby-Arnon scRNA-seq data (hg19 tpm) from GEO (GSE115978)"""
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


def load_tcga_skcm_hg19_normalized_counts_dereks_file() -> pd.DataFrame:
    """Load RNA-seq (hg19 normalized counts) for TCGA SKCM, processed by Derek"""
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/skcm_rnaseqv2_normalized_clean.txt"
    logger.debug(f"reading {path}")
    bulk_rna_seq = pd.read_csv(path, sep="\t", index_col=0, engine="pyarrow")
    # clean up index (gene symbols)
    bulk_rna_seq = bulk_rna_seq.sort_index()
    # clean up columns (sample IDs)
    bulk_rna_seq = bulk_rna_seq.rename(columns=lambda sample_id: sample_id.replace(".", "-"))
    bulk_rna_seq = bulk_rna_seq.reindex(sorted(bulk_rna_seq.columns), axis=1)
    # clean up everything else
    bulk_rna_seq = bulk_rna_seq.rename_axis(index=columns.GENE_SYMBOL, columns=columns.SAMPLE_ID)
    return bulk_rna_seq


def load_tcga_skcm_hg19_scaled_estimate_firebrowse() -> pd.DataFrame:
    """Load RNA-seq (hg19 scaled_estimate) for TCGA SKCM, from firebrowse"""
    path = "gs://liulab/firebrowse.org/SKCM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"
    logger.debug(f"reading {path}")
    df = pd.read_csv(
        path,
        sep="\t",
        header=[0, 1],
        index_col=0,
    ).iloc[29:]
    df = (
        # df.set_index(pd.MultiIndex.from_tuples(df.index.str.split("|").tolist()))
        df.set_index(df.index.str.replace(r"\|.*", "", regex=True))
        .rename_axis(index=columns.GENE_SYMBOL, columns=[columns.SAMPLE_ID, "values"])
        .loc(axis=1)[:, "scaled_estimate"]
        .groupby(columns.GENE_SYMBOL)
        .sum()
        .droplevel(level="values", axis=1)
        .sort_index()
    )
    return df


def load_tcga_skcm_fractions_from_csx() -> pd.DataFrame:
    logger.debug("loading TCGA SKCM fractions estimated by CIBERSORTx")
    path = "gs://liulab/downloaded_manually/derek_csx_tcga_skcm/CIBERSORTx_Job8_Results.txt"
    fractions = pd.read_csv(path, sep="\t", index_col=0)
    # clean up index (sample IDs)
    fractions = fractions.rename(index=lambda sample_id: sample_id.replace(".", "-"))
    fractions = fractions.sort_index()
    # clean up columns (cell type names)
    fractions = fractions.drop(fractions.columns[-3:], axis=1)
    fractions = fractions.rename(columns=weird_to_nice)
    fractions = fractions.reindex(sorted(fractions.columns), axis=1)
    # clean up everything else
    fractions = fractions.rename_axis(index=columns.SAMPLE_ID, columns=columns.CELL_TYPE)
    return fractions


def load_tcga_skcm_hg38_fpkm_bigquery() -> pd.DataFrame:
    query_string = """
    SELECT
        aliquot_barcode,
        gene_name,
        sum(HTSeq__FPKM) as HTSeq__FPKM_sum
    FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
    WHERE
        project_short_name = "TCGA-SKCM"
        and gene_type = 'protein_coding'
    GROUP BY 1, 2
    """
    client = bigquery.Client()
    logger.debug("loading TCGA SKCM bulk RNA-seq data from BigQuery")
    logger.debug("making query job")
    query_job = client.query(query_string)
    logger.debug("reading data from query to dataframe")
    bulk_rna_seq = query_job.to_dataframe(progress_bar_type="tqdm")
    logger.debug("pivoting results")
    bulk_rna_seq = bulk_rna_seq.pivot(index="gene_name", columns="aliquot_barcode", values="HTSeq__FPKM_sum")
    bulk_rna_seq = bulk_rna_seq.rename_axis(index=columns.GENE_SYMBOL, columns=columns.SAMPLE_ID)
    return bulk_rna_seq
