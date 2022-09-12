import numpy as np
import pandas as pd


def add_noise_multiplying_uniform(
    df_rnaseq: pd.DataFrame, rng: np.random.Generator = np.random.default_rng()
) -> pd.DataFrame:
    return df_rnaseq * rng.uniform(0.9, 1.1, size=(df_rnaseq.shape))


def normalize_expression(
    geps: pd.DataFrame, normalization_factor: int = 1_000_000
) -> pd.DataFrame:
    """Returns expression matrix with columns that sum to normalization_factor.

    Args:
        geps (pd.DataFrame): index is genes.
        normalization_factor (int, optional): Defaults to 1_000_000.

    Returns:
        pd.DataFrame: normalized gene expression.
    """
    return geps * normalization_factor / geps.sum()
