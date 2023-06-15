import logging
import re
from functools import partial
from typing import Union

import numpy as np
import pandas as pd
import upath
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def get_parquet_paths(path_root: upath.UPath) -> list[upath.UPath]:
    paths = path_root.glob("**/gene_stats_*")
    paths = map(str, paths)
    # paths = filter(re.compile(r".*run_id=0[0-1].*").match, paths)
    paths = sorted(paths)
    paths = list(paths)
    return paths


def extract_from_path(path: str, var_name: str) -> str:
    _ = path.split(var_name + "=")[1]
    return _.split("/")[0]


def test_extract_from_path():
    result = extract_from_path("thing/a=10", "a")
    assert result == "10"


def extract_vars_from_path(path: str) -> list[tuple[str, str]]:
    # extract the variables from the path
    return re.findall(r"(\w+)=(\w+)", path)


def test_extract_vars_from_path():
    test_path = "thing1=0/a=1,b=2/foo=bar/thing.parquet"
    result = extract_vars_from_path(test_path)
    expectation = [("thing1", "0"), ("a", "1"), ("b", "2"), ("foo", "bar")]
    assert result == expectation


def extract_origin_from_path(path: str) -> str:
    # extract the substring between "gene_stats_" and ".parquet"
    return re.findall(r"gene_stats_(.+)\.parquet", path)[0]


def test_extract_origin_from_path():
    test_path = "/foo=bar/gene_stats_thing.parquet"
    result = extract_origin_from_path(test_path)
    assert result == "thing"


ordering_functions = {
    "log2_fc": float,
    "malignant_means": lambda x: "0.715,0.715" if x == "None,None" else x,
}


def load_gene_stats(path_or_paths: Union[upath.UPath, list[upath.UPath]]) -> pd.DataFrame:
    if isinstance(path_or_paths, list):
        logger.debug("Loading from multiple paths: %s", path_or_paths)
        parquet_paths = [
            parquet_path
            for path_root in path_or_paths
            for parquet_path in get_parquet_paths(path_root)
        ]
        path_root = path_or_paths[0]
    else:
        path_root = path_or_paths
        logger.debug("Loading from single path: %s", path_root)
        parquet_paths = get_parquet_paths(path_root)
    logger.debug("Loading %d parquet files", len(parquet_paths))
    # parquet_paths = parquet_paths[::29]
    df = pd.concat(
        {str(path): pd.read_parquet(path) for path in parquet_paths},
        names=["path", "index"],
    )

    # add columns from path (e.g. "/foo=bar")
    for column_name in ["malignant_means", "log2_fc", "run_id"]:
        try:
            s = df.index.get_level_values("path").map(
                partial(extract_from_path, var_name=column_name)
            )
        except IndexError:
            logger.warning("skipping %s because not found", column_name)
            continue
        if column_name in ordering_functions:
            logger.debug("Setting ordering for %s", column_name)
            values = s.unique()
            values = sorted(values, key=ordering_functions[column_name])
            dtype = pd.CategoricalDtype(values, ordered=True)
            s = s.astype(dtype)
        df[column_name] = s

    # add "origin" column
    df["origin"] = df.index.get_level_values("path").map(extract_origin_from_path)

    # add "gene_perturbed" column
    genes_perturbed = pd.read_csv(
        path_root / "genes_perturbed.csv",
        index_col=0,
    )
    gene_in_perturbation_list = df["gene_symbol"].isin(genes_perturbed.index)
    gene_perturbation_was_nonzero = df["log2_fc"].astype(float).abs() > 0
    df["gene_perturbed"] = gene_in_perturbation_list & gene_perturbation_was_nonzero

    df["-log10_pval_signed_directional"] = get_neg_log10_pval_signed_directional(df)

    df = df.set_index(
        [
            "origin",
            "malignant_means",
            "log2_fc",
            "run_id",
            "gene_symbol",
            # "gene_perturbed",  # why should this be in the index?
        ]
    )
    df = df.sort_index()

    # add BH-adjusted p-value (aka q-value) and derived fields
    df_pval_adjusted = compute_pval_adjusted_fields(df)
    df = add_pval_adjusted_fields(df, df_pval_adjusted)
    return df


def get_neg_log10_pval_signed_directional(df: pd.DataFrame) -> pd.Series:
    # add signed p-value fields
    sign_observed_log2_fc = np.sign(df["log2_fold_change"].fillna(0)).replace({0: 1})
    sign_dist_log2_fc = np.sign(df["log2_fc"].astype(float)).replace({0: 1})
    # df["-log10_pval_signed"] already exists
    # df["-log10_pval_signed"] = df["-log10_pval"] * sign_observed_log2_fc
    return df["-log10_pval"] * sign_dist_log2_fc * sign_observed_log2_fc


def check_gene_stats(df_gene_stats: pd.DataFrame):
    # groupby_fields = "origin	malignant_means	log2_fc	run_id".split("\t")
    groupby_fields = list(filter(lambda x: x != "gene_symbol", df_gene_stats.index.names))
    df_gene_stats_groupby = df_gene_stats.groupby(groupby_fields)
    print(df_gene_stats_groupby.size())
    # grouping by origin, count number of unique values of gene_symbol
    print(
        df_gene_stats.groupby("origin").apply(
            lambda df: df.index.get_level_values("gene_symbol").nunique()
        )
    )
    for level in df_gene_stats.index.names:
        # print level and nunique
        print(level, df_gene_stats.index.get_level_values(level).nunique())
    return df_gene_stats


def compute_pval_adjusted_for_group(df: pd.DataFrame) -> pd.DataFrame:
    alpha = 0.1  # false discovery rate
    significance_column = f"significant_bh_fdr={alpha:.2f}"
    reject, pvals_corrected = multipletests(df["pval"], alpha=alpha, method="fdr_bh")[0:2]
    return pd.DataFrame(
        {
            significance_column: reject,
            "pval_adjusted_bh": pvals_corrected,
            "-log10_pval_adjusted_bh": -np.log10(pvals_corrected),
        },
        index=df.index,
    )


def compute_pval_adjusted_fields(
    df_gene_stats: pd.DataFrame, groupby_fields: list[str] = []
) -> pd.DataFrame:
    observed_log2_fc = df_gene_stats["log2_fold_change"].fillna(0)
    sign_observed_log2_fc = np.sign(observed_log2_fc).replace({0: 1})
    if "gene_symbol" in df_gene_stats.index.names:
        dist_log2_fc = df_gene_stats.index.to_frame()["log2_fc"].astype(float)
        # get list of every index level name except "gene_symbol"
        if groupby_fields is None:
            groupby_fields = [x for x in df_gene_stats.index.names if x != "gene_symbol"]
    else:
        dist_log2_fc = df_gene_stats["log2_fc"]
        if groupby_fields is None:
            raise ValueError(
                "groupby_fields must be specified if index doesn't have partition keys"
            )
    sign_dist_log2_fc = np.sign(dist_log2_fc).replace({0: 1})
    # assert only 1 and/or -1 in sign_dist_log2_fc
    assert set(sign_dist_log2_fc.unique()).union({-1, 1}) == {1, -1}
    logger.debug("Grouping by %s", groupby_fields)
    # don't prepend groupby fields to the resulting index
    dfg = df_gene_stats.groupby(groupby_fields, group_keys=False)
    # sanity check - assert size of each group is 5000
    logger.debug("Counts of group sizes: %s", dfg.size().value_counts())
    df_pval_adjusted = dfg.apply(compute_pval_adjusted_for_group)
    assert df_gene_stats.index.equals(df_pval_adjusted.index), "Index of result is different"
    df_pval_adjusted["pval_adjusted_bh_signed"] = (
        df_pval_adjusted["pval_adjusted_bh"] * sign_observed_log2_fc
    )
    df_pval_adjusted["pval_adjusted_bh_signed_directional"] = (
        df_pval_adjusted["pval_adjusted_bh"] * sign_dist_log2_fc * sign_observed_log2_fc
    )
    df_pval_adjusted["-log10_pval_adjusted_bh_signed"] = (
        df_pval_adjusted["-log10_pval_adjusted_bh"] * sign_observed_log2_fc
    )
    df_pval_adjusted["-log10_pval_adjusted_bh_signed_directional"] = (
        df_pval_adjusted["-log10_pval_adjusted_bh"] * sign_dist_log2_fc * sign_observed_log2_fc
    )
    return df_pval_adjusted


def add_pval_adjusted_fields(
    df_gene_stats: pd.DataFrame, df_pval_adjusted: pd.DataFrame
) -> pd.DataFrame:
    assert df_gene_stats.index.equals(df_pval_adjusted.index), "Index mismatch"
    # fields in df_pval_adjusted that are not in df_gene_stats
    fields = filter(lambda x: x not in df_gene_stats.columns, df_pval_adjusted.columns)
    additional_fields = df_pval_adjusted[fields]
    logger.debug("adding fields %s", list(additional_fields.columns))
    return df_gene_stats.join(additional_fields, how="inner", validate="1:1")


def add_more_pval_fields(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add signed_directional and -log10 pval fields
    """
    df = df.copy()
    df["-log10_pval_signed_directional"] = get_neg_log10_pval_signed_directional(df)
    df_pval_adjusted = compute_pval_adjusted_fields(
        df, ["origin", "malignant_means", "log2_fc", "run_id"]
    )
    df = add_pval_adjusted_fields(df, df_pval_adjusted)
    return df
