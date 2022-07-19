import logging

import pandas as pd

import helpers
import helpers.columns

logger = logging.getLogger(__name__)


def create_csx_mixtures_tsv(
    df_bulk_rnaseq: pd.DataFrame, uri_save_destination: str
) -> str:
    # tmp_dir = tmpfile.mkdtemp()
    # uri_csv_file = tmp_dir / "bulk_rnaseq.csv"
    # uri_csv_file = tempfile.mkstemp()
    df_bulk_rnaseq = df_bulk_rnaseq.rename_axis(index="GeneSymbol", columns=None)
    df_bulk_rnaseq = df_bulk_rnaseq.reset_index()
    df_bulk_rnaseq.to_csv(uri_save_destination, index=False, sep="\t")
    return uri_save_destination


def create_csx_refsample_tsv(
    df_refsample_sc_rnaseq: pd.DataFrame,
    df_refsample_sc_metadata: pd.DataFrame,
    uri_save_destination: str,
) -> str:
    # tmp_dir = tempfile.mkdtemp()
    # uri_csv_file = os.path.join(tmp_dir, "refsample_sc_rnaseq.tsv")
    # uri_csv_file = tempfile.mkstemp()

    def look_up_cell_type(single_cell_id):
        return df_refsample_sc_metadata.loc[single_cell_id][helpers.columns.CELL_TYPE]

    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename(columns=look_up_cell_type)
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename_axis(
        index="GeneSymbol", columns=None
    )
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.reset_index()
    df_refsample_sc_rnaseq.to_csv(uri_save_destination, index=False, sep="\t")
    return uri_save_destination


def create_csx_fractions_tsv():
    raise NotImplementedError


def create_refsample_from_jerby_arnon(uri_save_destination):
    logger.debug("loading sc data")
    df_sc_rnaseq, df_metadata = helpers.datasets.load_jerby_arnon_hg19_tpm()
    df_sc_rnaseq *= 1e6 / df_sc_rnaseq.sum()
    logger.debug(f"shape of scRNA-seq: {df_sc_rnaseq.shape}")
    nonnull_cell_type_cells = df_metadata[
        ~df_metadata[helpers.columns.CELL_TYPE].isna()
    ].index
    logger.debug(f"length of nonnull_cell_type_cells: {len(nonnull_cell_type_cells)}")
    logger.debug(f"{nonnull_cell_type_cells}")
    df_sc_rnaseq = df_sc_rnaseq[nonnull_cell_type_cells]
    logger.debug(f"new shape of scRNA-seq: {df_sc_rnaseq.shape}")
    logger.debug("creating and writing refsample tsv")
    create_csx_refsample_tsv(df_sc_rnaseq, df_metadata, uri_save_destination)


def create_csx_mixtures_for_tcga_skcm(uri_save_destination):
    logger.debug("loading bulk rna-seq for tcga skcm")
    df_bulk_rnaseq = pd.read_parquet(
        "gs://liulab/data/pseudobulk_optimization/3_with_tcga_qc/mixtures_real_tcga_skcm/tpm.parquet"
    )
    df_bulk_rnaseq = df_bulk_rnaseq.pivot("gene_symbol", "aliquot_barcode", "tpm")
    logger.debug("creating tsv of bulk rna-seq for tcga skcm")
    create_csx_mixtures_tsv(
        df_bulk_rnaseq,
        uri_save_destination,
    )


def create_csx_mixtures_for_pseudobulk(uri_save_destination):
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
        uri_save_destination,
    )


if __name__ == "__main__":
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(process)d/%(threadName)s %(name)s %(levelname)s\n%(message)s"
        )
    )
    logging.getLogger().addHandler(handler)
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.datasets").setLevel("DEBUG")

    uri_refsample_jerby_arnon = "gs://liulab/data/pseudobulk_evaluation/csx_input_files/refsample_jerby_arnon_hg19_tpm.tsv"
    create_refsample_from_jerby_arnon(uri_refsample_jerby_arnon)

    logger.debug("creating csx mixtures for tcga skcm")
    create_csx_mixtures_for_tcga_skcm(
        "gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_tcga_skcm.tsv"
    )

    logger.debug("creating csx mixtures for pseudobulk")
    create_csx_mixtures_for_pseudobulk(
        "gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_pseudobulk.tsv"
    )
