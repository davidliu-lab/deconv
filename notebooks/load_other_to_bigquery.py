import logging

import duckdb
import numpy as np
import pandas as pd
from google.cloud import bigquery

import helpers

from helpers.logging.configuring_logging import EasternTimeFormatter
from helpers.running_cibersortx.reading_input_files import get_arrow_dataset

logger = logging.getLogger(__name__)


def load_like_cell_type_geps(
    df: pd.DataFrame, client: bigquery.Client, table_id: str
) -> bigquery.LoadJob:
    job_config = bigquery.LoadJobConfig(
        schema=[
            bigquery.SchemaField("malignant_means", bigquery.enums.SqlTypeNames.STRING),
            bigquery.SchemaField("log2_fc", bigquery.enums.SqlTypeNames.FLOAT64),
            bigquery.SchemaField("run_id", bigquery.enums.SqlTypeNames.INTEGER),
            bigquery.SchemaField("name", bigquery.enums.SqlTypeNames.STRING),
        ],
        # write_disposition="WRITE_TRUNCATE",  # overwrite
    )
    logger.debug("creating job to load dataframe to bigquery")
    job = client.load_table_from_dataframe(
        df, table_id, job_config=job_config
    )  # Make an API request.
    return job


if __name__ == "__main__":
    handler = logging.StreamHandler()
    LOGGING_FORMAT = "%(asctime)s %(levelname)s %(name)s %(message)s"
    handler.setFormatter(logging.Formatter(LOGGING_FORMAT, datefmt="%I:%M:%S %p"))
    logging.getLogger().addHandler(handler)
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("google.cloud").setLevel("DEBUG")
    client = bigquery.Client()
    fractions_arrow_dataset = get_arrow_dataset(
        "gs://liulab/differential_composition_and_expression/20230616_03h34m20s/fractions/"
    )
    for i, batch in enumerate(fractions_arrow_dataset.to_batches()):
        df = batch.to_pandas()
        job = load_like_cell_type_geps(df, client, "deconv.fractions")
        logger.info("waiting for job to finish")
        print(job.result())
        logger.info("job finished")
    # cell_type_geps_arrow_dataset = get_arrow_dataset(
    #     "gs://liulab/differential_composition_and_expression/20230616_03h34m20s/cell_type_geps/"
    # )
    # # upload one batch at a time...
    # partitions = cell_type_geps_arrow_dataset.to_batches()
    # for i, partition in enumerate(partitions):
    #     logger.info(f"uploading partition {i}")
    #     df = partition.to_pandas()
    #     df = df.rename(columns={"cell_type": "cell_type_name"})
    #     df["cell_type"] = df["cell_type_name"].apply(
    #         lambda x: helpers.running_cibersortx.reading_input_files.cell_type_name_to_cell_type_id[
    #             x
    #         ]
    #     )
    #     df["run_id"] = 1
    #     df["log2_fc"] = 0
    #     df = df[["cell_type", "cell_type_name", "run_id", "log2_fc"]]
    #     job = load_deg_results(df, client, "deconv.deg_analysis_results")
    #     logger.info("waiting for job to finish")
    #     print(job.result())
    #     logger.info("job finished")
    
    # rows_per_chunk = 16063 * 4296 // 24  # integer division
    # assert rows_per_chunk * 24 == 16063 * 4296
    # start = 0
    # while start < len(df_gene_stats):
    #     end = start + rows_per_chunk
    #     logger.info(f"uploading rows {start} to {end} of {len(df_gene_stats)}")
    #     job = load_deg_results(df_gene_stats.iloc[start:end], client, "deconv.deg_analysis_results")
    #     logger.info("waiting for job to finish")
    #     print(job.result())
    #     logger.info("job finished")
    #     start = end
