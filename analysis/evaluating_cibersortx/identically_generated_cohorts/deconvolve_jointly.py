import logging

import pandas as pd
import upath

import helpers

logger = logging.getLogger(__name__)


def read_simulated_bulkrnaseq(name):
    path_root = upath.UPath(
        "gs://liulab/evaluating_cibersortx/identical_cohorts/2022-08-25_05:53:48/"
    )
    df_simulated_bulkrnaseq = pd.read_parquet(
        path_root / name / "all_simulated_data/simulated_bulkrnaseq.parquet"
    )
    return df_simulated_bulkrnaseq.add_prefix(f"cohort{name}-")


if __name__ == "__main__":
    df_simulated_bulkrnaseq_fused = pd.concat(
        [read_simulated_bulkrnaseq("0"), read_simulated_bulkrnaseq("1")],
        axis="columns",
    )
    path_results = upath.UPath(
        "gs://liulab/evaluating_cibersortx/identical_cohorts/2022-08-25_05:53:48/running_both_jointly"
    )
    path_bulk_rnaseq = path_results / "pseudobulk_rnaseq.tsv"
    helpers.running_cibersortx.creating_input_files.create_csx_mixtures_tsv(
        df_simulated_bulkrnaseq_fused, path_bulk_rnaseq
    )
    helpers.running_cibersortx.hires_with_fractions.run_and_upload(
        str(path_results),
        path_bulk_rnaseq,
        uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )
