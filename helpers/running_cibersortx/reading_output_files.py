import pathlib
from typing import Union

import cloudpathlib
import dask.dataframe as dd


def read_hires_cell_type_geps(
    path_pattern: Union[pathlib.Path, cloudpathlib.CloudPath]
):
    """Read cell type-specific gene expression outputs from CIBERSORTx hires mode.

    Parameters
    ----------
    path_pattern : pathlib.Path or cloudpathlib.CloudPath
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


def read_fractions(path: Union[pathlib.Path, cloudpathlib.CloudPath]):
    raise NotImplementedError


if __name__ == "__main__":
    path = cloudpathlib.CloudPath(
        "gs://liulab/evaluating_cibersortx/perturbed_gene_expression/2x/2022-09-14_15:17:34"
    )
    path_pattern = path / "outdir" / "CIBERSORTxHiRes_NA_*_Window*txt"
    df = read_hires_cell_type_geps(path_pattern)
    print(df)
