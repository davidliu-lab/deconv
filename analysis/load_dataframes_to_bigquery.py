import logging

import duckdb
import numpy as np
import pandas as pd
from google.cloud import bigquery

import helpers
from helpers.deg_analysis import (
    classifier_metrics,
    displaying_tables,
    plotting_curves,
    plotting_volcanos,
)
from helpers.deg_analysis.classifier_metrics_old import (
    calculate_all_curves,
    compute_all_curves_and_metrics,
)
from helpers.deg_analysis.loading_results import (
    get_arrow_dataset_for_deg_analysis_results,
)
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields

logger = logging.getLogger(__name__)


def load_deg_results(df_gene_stats: pd.DataFrame, client: bigquery.Client, table_id: str):
    job_config = bigquery.LoadJobConfig(
        schema=[
            bigquery.SchemaField("log2_fc", bigquery.enums.SqlTypeNames.FLOAT64),
            bigquery.SchemaField("run_id", bigquery.enums.SqlTypeNames.INTEGER),
        ],
        write_disposition="WRITE_TRUNCATE",  # overwrite
    )
    logger.debug("creating job to load dataframe to bigquery")
    job = client.load_table_from_dataframe(
        df_gene_stats, table_id, job_config=job_config
    )  # Make an API request.
    return job


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format=helpers.logging.FORMAT)
    logging.getLogger("helpers").setLevel("DEBUG")
    client = bigquery.Client()
    gene_stats_arrow_dataset = (
        helpers.deg_analysis.loading_results.get_arrow_dataset_for_deg_analysis_results(
            "gs://liulab/differential_composition_and_expression/20230616_03h34m20s/deg_analysis/"
        )
    )
    query_text = """
SELECT
    origin,
    malignant_means,
    log2_fc,
    run_id,
    gene_symbol,
    perturbed AND log2_fc != 0 AS perturbed,
    log2_fold_change,
    "pval",
    "-log10_pval",
    "pval_adjusted_bh",
    -1.0 * log10("pval_adjusted_bh") as "-log10_pval_adjusted_bh",
    "significant_bh_fdr=0.10",
FROM gene_stats_arrow_dataset
--WHERE
--    origin = 'malignant_cibersortx'
--    AND malignant_means in ('0.6,0.8',)
--    AND log2_fc in (1.0,)
--    AND run_id in (0,)
--    AND gene_symbol like '[A-K]%'  -- start with A-K
;
    """
    logger.debug("querying data from arrow dataset")
    df_gene_stats = duckdb.sql(query_text).df()  # query with "duckdb:///:default:"
    logger.debug("adding more pval fields")
    df_gene_stats = add_more_pval_fields(df_gene_stats)
    logger.debug("renaming some columns")
    # rename columns so that "-log10" is replaced with "neg_log10"
    df_gene_stats = df_gene_stats.rename(columns=lambda x: x.replace("-log10", "neg_log10"))
    # rename column "significant_bh_fdr=0.10" to "significant_bh_fdr"
    df_gene_stats = df_gene_stats.rename(columns={"significant_bh_fdr=0.10": "significant_bh_fdr"})
    logger.debug("calling load_deg_results")
    job = load_deg_results(df_gene_stats, client, "deconv.deg_analysis_results")
    logger.info("waiting for job to finish")
    print(job.result())
    logger.info("job finished")
