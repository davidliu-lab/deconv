import pandas as pd


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
