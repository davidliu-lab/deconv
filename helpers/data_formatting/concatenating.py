import logging

import cloudpathlib
import pandas as pd

logger = logging.getLogger(__name__)


def load_and_concatenate_bulk_rnaseq(
    cohort_paths: dict[str, cloudpathlib.CloudPath]
) -> pd.DataFrame:
    dataframes = {
        cohort: pd.read_parquet(path / "bulk_rnaseq.parquet")
        for cohort, path in cohort_paths.items()
    }
    df = pd.concat(dataframes, axis="columns")
    column_axis_name = df.columns.names[1]
    df.columns = df.columns.map("/".join)
    df.rename_axis(columns=column_axis_name, inplace=True)
    logger.debug("Shape of concatenated data: %s", df.shape)
    return df


def load_and_concatenate_fractions(
    cohort_paths: dict[str, cloudpathlib.CloudPath]
) -> pd.DataFrame:
    dataframes = {
        cohort: pd.read_parquet(path / "fractions.parquet")
        for cohort, path in cohort_paths.items()
    }
    df = pd.concat(dataframes, axis="rows")
    index_axis_name = df.index.names[1]
    df.index = df.index.map("/".join)
    df.rename_axis(index=index_axis_name, inplace=True)
    logger.debug("Shape of concatenated data: %s", df.shape)
    return df
