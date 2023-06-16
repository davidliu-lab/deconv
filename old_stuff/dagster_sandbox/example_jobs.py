import dagster
import pandas as pd
from dagster import (
    AssetIn,
    Definitions,
    Out,
    StaticPartitionsDefinition,
    asset,
    job,
    op,
)
from dagster_gcp_pandas import BigQueryPandasIOManager, bigquery_pandas_io_manager

io_manager = BigQueryPandasIOManager(
    project="keen-dispatch-316219",
    dataset="IRIS",
    timeout=15.0,
)


@asset(code_version="20230523.1")
def iris_data() -> pd.DataFrame:
    return pd.read_csv(
        "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data",
        names=[
            "sepal_length_cm",
            "sepal_width_cm",
            "petal_length_cm",
            "petal_width_cm",
            "species",
        ],
    )


@asset(
    ins={
        "iris_sepal": AssetIn(
            key="iris_data",
            metadata={"columns": ["sepal_length_cm", "sepal_width_cm"]},
        )
    }
)
def sepal_data(iris_sepal: pd.DataFrame) -> pd.DataFrame:
    iris_sepal["sepal_area_cm2"] = iris_sepal["sepal_length_cm"] * iris_sepal["sepal_width_cm"]
    return iris_sepal


@asset
def iris_cleaned(iris_data: pd.DataFrame) -> pd.DataFrame:
    return iris_data.dropna().drop_duplicates()


@asset(
    partitions_def=StaticPartitionsDefinition(["Iris-setosa", "Iris-virginica", "Iris-versicolor"]),
    metadata={"partition_expr": "SPECIES"},
)
def iris_data_partitioned(context) -> pd.DataFrame:
    species = context.asset_partition_key_for_output()

    full_df = pd.read_csv(
        "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data",
        names=[
            "sepal_length_cm",
            "sepal_width_cm",
            "petal_length_cm",
            "petal_width_cm",
            "species",
        ],
    )

    return full_df[full_df["species"] == species]


defs = Definitions(
    assets=[iris_data, sepal_data, iris_cleaned, iris_data_partitioned],
    resources={"io_manager": io_manager},
)
