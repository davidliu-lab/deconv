import logging
import re

import pandas as pd

FIELD_MATCHER = re.compile(r"(^[^a-zA-Z])|(\W)")

logging.basicConfig(level="INFO")
logger = logging.getLogger(__name__)


def make_column_name_valid(column):
    original = column
    column = column.lower()
    column = FIELD_MATCHER.sub("_", column)
    logger.info(f"converted {original} -> {column}")
    return column


logger.info(f"reading and loading cell_annotations")
pd.read_csv("/mnt/buckets/liulab/ftp/GSE115978/GSE115978_cell.annotations.csv").rename(
    mapper=make_column_name_valid, axis=1
).to_gbq("gse115978_tirosh.cell_annotations", if_exists="replace")

logger.info(f"reading tpm")
tpm = pd.read_csv(
    "/mnt/buckets/liulab/ftp/GSE115978/GSE115978_tpm.csv", index_col=0
).rename_axis("gene", axis="index")

# doesn't work because of duplicate column names
# logger.info(f"loading tpm")
# tpm.to_gbq("gse115978_tirosh.tpm", if_exists="replace")

logger.info(f"making tpm_unpivoted")
tpm_unpivoted = (
    tpm.rename_axis("single_cell", axis="columns")
    .stack()
    .reset_index(name="expression")[["single_cell", "gene", "expression"]]
)
logger.info(f"loading tpm_unpivoted")
tpm_unpivoted.to_gbq(
    "gse115978_tirosh.tpm_unpivoted", chunksize=1_000_000, if_exists="replace"
)
