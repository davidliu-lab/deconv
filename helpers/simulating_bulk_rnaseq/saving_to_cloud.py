import logging

import pandas as pd
import upath

logger = logging.getLogger(__name__)


def save_simulated_data(
    df_simulated_bulkrnaseq: pd.DataFrame,
    simulated_cell_type_geps: dict[pd.DataFrame],
    path_root: upath.UPath,
) -> None:
    save_simulated_bulkrnaseq(df_simulated_bulkrnaseq, path_root)
    save_simulated_cell_type_geps(simulated_cell_type_geps, path_root)


def save_simulated_bulkrnaseq(
    df_simulated_bulkrnaseq: pd.DataFrame,
    path_root: upath.UPath,
) -> None:
    uri_simulated_bulkrnaseq = str(path_root / "simulated_bulkrnaseq.parquet")
    logger.debug("saving simulated bulkrnaseq to %s", uri_simulated_bulkrnaseq)
    df_simulated_bulkrnaseq.to_parquet(uri_simulated_bulkrnaseq)


def save_simulated_cell_type_geps(
    simulated_cell_type_geps: dict[pd.DataFrame],
    path_root: upath.UPath,
) -> None:
    df_simulated_cell_type_geps = pd.concat(simulated_cell_type_geps, names=["sample_id"])
    uri_simulated_cell_type_geps = str(path_root / "simulated_cell_type_geps.parquet")
    logger.debug("saving simulated cell type geps to %s", uri_simulated_cell_type_geps)
    df_simulated_cell_type_geps.to_parquet(uri_simulated_cell_type_geps)
