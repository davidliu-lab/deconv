import logging

import pandas as pd

import helpers.columns


def create_csx_mixture_tsv(df_bulk_rnaseq: pd.DataFrame, uri_file: str) -> str:
    # tmp_dir = tmpfile.mkdtemp()
    # uri_csv_file = tmp_dir / "bulk_rnaseq.csv"
    # uri_csv_file = tempfile.mkstemp()
    df_bulk_rnaseq = df_bulk_rnaseq.rename_axis(index="Gene", columns=None)
    df_bulk_rnaseq *= 1e6 / df_bulk_rnaseq.sum()
    df_bulk_rnaseq = df_bulk_rnaseq.reset_index()
    df_bulk_rnaseq.to_csv(uri_file, index=False, sep="\t")
    return uri_file


def create_csx_refsample_tsv(
    df_refsample_sc_rnaseq: pd.DataFrame,
    df_refsample_sc_metadata: pd.DataFrame,
    uri_file: str,
) -> str:
    # tmp_dir = tempfile.mkdtemp()
    # uri_csv_file = os.path.join(tmp_dir, "refsample_sc_rnaseq.tsv")
    # uri_csv_file = tempfile.mkstemp()

    def look_up_cell_type(single_cell_id):
        return df_refsample_sc_metadata.loc[single_cell_id][helpers.columns.CELL_TYPE]

    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename(columns=look_up_cell_type)
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename_axis(
        index="Gene", columns=None
    )
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.reset_index()
    df_refsample_sc_rnaseq.to_csv(uri_file, index=False, sep="\t")
    return uri_file


if __name__ == "__main__":
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logger.debug("loading sc data")
    (
        df_sc_rnaseq,
        df_metadata,
    ) = helpers.datasets.load_jerby_arnon_hg19_tpm()
    logger.debug(f"length of df_jerby_arnon_sc_rnaseq: {len(df_sc_rnaseq)}")
    nonnull_cell_type_cells = ~df_metadata[helpers.columns.CELL_TYPE].isna()
    df_sc_rnaseq = df_sc_rnaseq.loc[nonnull_cell_type_cells]
    logger.debug(f"new length of df_jerby_arnon_sc_rnaseq: {len(df_sc_rnaseq)}")
    logger.debug("creating and writing refsample tsv")
    uri_file = "gs://liulab/tmp/refsample_from_helpers.tsv"
    create_csx_refsample_tsv(df_sc_rnaseq, df_metadata, uri_file)
    logger.debug("validating refsample tsv")
    df = pd.read_csv(uri_file, sep="\t")
    logger.debug(f"length of df: {len(df)}")
    print("\n".join(df.columns))
