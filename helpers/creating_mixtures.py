import logging

import numpy as np
import pandas as pd

from helpers import columns

logger = logging.getLogger(__name__)


def make_cell_type_geps(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    n_cells_per_gep: int = 5,
    normalization_factor: int = 1_000_000,
    malignant_from_one_sample: bool = True,
    rng: np.random.Generator = np.random.default_rng(),
):
    if malignant_from_one_sample:
        malignant_cells = sc_metadata[columns.CELL_TYPE] == "Malignant"
        random_sample = rng.choice(sc_metadata[malignant_cells][columns.SAMPLE_ID].unique())
        logger.debug(f"randomly chose {random_sample} for malignant cells")
        single_cells_to_include = np.logical_not(np.logical_and(malignant_cells, sc_metadata[columns.SAMPLE_ID] != random_sample))
        sc_metadata = sc_metadata[single_cells_to_include]
    sampled_single_cells_per_type = (
        sc_metadata.groupby(columns.CELL_TYPE)
        .apply(lambda group: list(rng.choice(group[columns.SINGLE_CELL_ID], n_cells_per_gep)))
        .to_dict()
    )
    cell_type_geps = pd.concat(
        {cell_type: sc_data[cells].sum(axis="columns") for cell_type, cells in sampled_single_cells_per_type.items()},
        axis="columns",
    )
    cell_type_geps = normalize_expression(cell_type_geps, normalization_factor)
    return cell_type_geps


def make_fractions_from_dirichlet(n_samples: int, sc_metadata: pd.DataFrame, rng: np.random.Generator):
    samples = pd.Series(
        [f"sample_{j:0{len(str(n_samples))}d}" for j in range(n_samples)],
        name=columns.SAMPLE_ID,
    )
    cell_type_groups = sc_metadata.groupby(columns.CELL_TYPE)
    fraction_values = rng.dirichlet(alpha=(0.5,) * cell_type_groups.ngroups, size=(n_samples,))
    fractions = pd.DataFrame(fraction_values, index=samples, columns=list(cell_type_groups.groups))
    return fractions


def make_mixtures(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    sample_fractions: pd.DataFrame,
    n_cells_per_gep: int = 5,
    normalization_factor: int = 1_000_000,
    malignant_from_one_sample: bool = True,
    rng: np.random.Generator = np.random.default_rng(),
):
    logger.debug(f"using np.random.Generator with BitGenerator state {rng.bit_generator.state['state']}")
    cell_type_geps = {
        sample: make_cell_type_geps(
            sc_data, sc_metadata, n_cells_per_gep, malignant_from_one_sample=malignant_from_one_sample, rng=rng
        )
        for sample in sample_fractions.index
    }
    mixtures = pd.concat(
        {sample: cell_type_geps[sample] @ sample_fractions.loc[sample] for sample in sample_fractions.index},
        axis=1,
    )
    mixtures = add_noise_multiplying_uniform(mixtures, rng)
    mixtures = normalize_expression(mixtures, normalization_factor)
    return mixtures, cell_type_geps


def add_noise_multiplying_uniform(mixtures: pd.DataFrame, rng: np.random.Generator = np.random.default_rng()):
    return mixtures * rng.uniform(0.9, 1.1, size=(mixtures.shape))


def normalize_expression(geps: pd.DataFrame, normalization_factor: int = 1_000_000):
    return geps * normalization_factor / geps.sum()
