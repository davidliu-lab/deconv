import numpy as np
import pandas as pd
from helpers.simulating_bulk_rnaseq.creating_mixtures import make_mixtures


def simulate_data(
    sc_rnaseq: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    fractions: pd.DataFrame,
    rng: np.random.Generator,
    n: int = 50,
    name: str = "control",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_index = pd.Index(
        [f"{name}/sample_{i:03d}" for i in range(n)], name="sample_id"
    )
    sampled_fractions = fractions.sample(n, replace=True, random_state=rng)
    sampled_fractions = sampled_fractions.set_index(sample_index)
    bulk_rnaseq, cell_type_geps = make_mixtures(
        sc_rnaseq, sc_metadata, sampled_fractions, rng=rng
    )
    cell_type_geps = pd.concat(cell_type_geps, names=["sample_id"])
    return sampled_fractions, bulk_rnaseq, cell_type_geps
