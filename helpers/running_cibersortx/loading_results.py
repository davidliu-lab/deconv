import logging

import pyarrow.dataset as ds

logger = logging.getLogger(__name__)


def get_arrow_dataset_for_deg_analysis_results(
    base_path: str = "gs://liulab/differential_composition_and_expression/copied/20230505_21h41m44s/deg_analysis/",
) -> ds.Dataset:
    """Get an arrow dataset for the results of a CIBERSORTx run.

    Args:
        base_path: The base path to the CIBERSORTx run.

    Returns:
        An arrow dataset for the results of a CIBERSORTx run.
    """
    logger.debug("loading arrow dataset for results of CIBERSORTx run at %s", base_path)
    return ds.dataset(
        base_path,
        format="parquet",
        partitioning="hive",
    )
