import pathlib
from typing import Union

import dask.dataframe as dd
import upath


def read_hires_cell_type_geps(path_pattern: Union[pathlib.Path, upath.UPath]):
    """Read cell type-specific gene expression outputs from CIBERSORTx hires mode.

    Parameters
    ----------
    path_pattern : pathlib.Path or upath.UPath
        path with wildcards to file(s) to read.

    Returns
    -------
    pandas.DataFrame
        The high-resolution cell type gene expression profiles.
    """
    df = dd.read_csv(str(path_pattern), sep="\t", include_path_column=True).compute()
    df[["cell_type"]] = df["path"].str.extract(
        ".*/CIBERSORTxHiRes_NA_(.*)_Window.*\.txt"
    )
    df = df.drop(columns="path")
    df = df.pivot(index=["cell_type", "GeneSymbol"], columns=[])
    df = df.rename_axis(columns="sample_id")
    df = df.stack()
    return df


def read_fractions(path: Union[pathlib.Path, upath.UPath]):
    raise NotImplementedError
