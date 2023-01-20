import numpy as np
import pandas as pd
from scipy.optimize import root_scalar


def scale_malignant_and_normalize(p: pd.DataFrame, scalar: float):
    p = p.copy()
    p["Malignant"] *= scalar
    p = p.div(p.sum(axis="columns"), axis="index")
    return p


def perturb_malignant_fractions(p: pd.DataFrame, desired_mean: float):
    """Perturb malignant fractions of a dataframe of fractions so that their mean is a desired value."""
    optimand = (
        lambda log2_c: scale_malignant_and_normalize(p, 2**log2_c)["Malignant"].mean()
        - desired_mean
    )
    result = root_scalar(optimand, method="secant", x0=-0.1, x1=0.1)
    log2_c = result.root
    p_perturbed = scale_malignant_and_normalize(p, 2**log2_c)
    np.testing.assert_almost_equal(p_perturbed["Malignant"].mean(), desired_mean)
    np.testing.assert_almost_equal(p_perturbed.sum(axis="columns").to_numpy(), 1.0)
    return p_perturbed


def compute_ratio(p1: pd.DataFrame, p2: pd.DataFrame):
    ratio = p1["Malignant"].mean() / p2["Malignant"].mean()
    return ratio


def compute_ratio_given_perturbation(p1: pd.DataFrame, p2: pd.DataFrame, c: float):
    p1 = scale_malignant_and_normalize(p1, c)
    p2 = scale_malignant_and_normalize(p2, 1 / c)
    ratio = compute_ratio(p1, p2)
    # print("ratio: ", ratio)
    return ratio


def perturb_pair_of_fractions(p1: pd.DataFrame, p2: pd.DataFrame, desired_ratio: float):
    optimand = lambda c: (compute_ratio_given_perturbation(p1, p2, c) - desired_ratio)
    result = root_scalar(optimand, method="secant", x0=1.0, x1=1.1)
    c = result.root
    p1 = scale_malignant_and_normalize(p1, c)
    p2 = scale_malignant_and_normalize(p2, 1 / c)
    df = pd.concat([p1, p2], keys=["higher", "lower"], names=["polarity"])
    return df
