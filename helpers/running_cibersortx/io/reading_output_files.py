import pathlib
from typing import Union

import cloudpathlib
import dask.dataframe as dd


def read_hires_cell_type_geps(path: Union[pathlib.Path, cloudpathlib.CloudPath]):
    """Read the high-resolution cell type gene expression profiles from a file.

    Parameters
    ----------
    path : pathlib.Path or cloudpathlib.CloudPath
        The path to the file to read.

    Returns
    -------
    pandas.DataFrame
        The high-resolution cell type gene expression profiles.
    """
    path_pattern = path / "outdir" / "CIBERSORTxHiRes_NA_*_Window*txt"
    df = dd.read_csv(str(path_pattern), sep="\t", include_path_column=True).compute()
    df[["cell_type"]] = df["path"].str.extract(
        ".*/outdir/CIBERSORTxHiRes_NA_(.*)_Window.*\.txt"
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
    df = read_hires_cell_type_geps(path)
    print(df)
