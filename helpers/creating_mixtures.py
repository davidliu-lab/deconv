import numpy as np
import pandas as pd


def make_cell_type_geps(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    n_cells_per_gep: int = 5,
    malignant_from_one_sample: bool = True,
    rng: np.random.Generator = np.random.default_rng(),
):
    if malignant_from_one_sample:
        malignant_cells = sc_metadata["cell.types"] == "Malignant"
        random_sample = rng.choice(sc_metadata[malignant_cells]["samples"].unique())
        single_cells_to_exclude = malignant_cells and (
            sc_metadata["samples"] != random_sample
        )
        sc_metadata = sc_metadata[not single_cells_to_exclude]
    sampled_single_cells_per_type = (
        sc_metadata.groupby("cell.types")
        .apply(lambda group: list(rng.choice(group["cells"], n_cells_per_gep)))
        .to_dict()
    )
    cell_type_geps = pd.concat(
        {
            cell_type: sc_data[cells].sum(axis="columns")
            for cell_type, cells in sampled_single_cells_per_type.items()
        },
        axis="columns",
    )
    cell_type_geps = normalize_to_tp100k(cell_type_geps)
    return cell_type_geps


def make_fractions_from_dirichlet(
    n_samples: int, sc_metadata: pd.DataFrame, rng: np.random.Generator
):
    samples = pd.Series(
        [f"sample_{j:0{len(str(n_samples))}d}" for j in range(n_samples)],
        name="Mixture",
    )
    cell_type_groups = sc_metadata.groupby("cell.types")
    fraction_values = rng.dirichlet(
        alpha=(0.5,) * cell_type_groups.ngroups, size=(n_samples,)
    )
    fractions = pd.DataFrame(
        fraction_values, index=samples, columns=list(cell_type_groups.groups)
    )
    return fractions


def make_mixtures(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    sample_fractions: pd.DataFrame = None,
    n_cells_per_gep: int = 5,
    malignant_from_one_sample: bool = True,
    rng: np.random.Generator = np.random.default_rng(),
):
    if sample_fractions is None:
        sample_fractions = make_fractions_from_dirichlet(n_samples, sc_metadata, rng)
    cell_type_geps = {
        sample: make_cell_type_geps(sc_data, sc_metadata, n_cells_per_gep, rng)
        for sample in sample_fractions.index
    }
    mixtures = pd.concat(
        {
            sample: cell_type_geps[sample] @ sample_fractions.loc[sample]
            for sample in sample_fractions.index
        },
        axis=1,
    )
    mixtures = add_noise_multiplying_uniform(mixtures, rng)
    mixtures = normalize_to_tp100k(mixtures)
    return mixtures, cell_type_geps


def add_noise_multiplying_uniform(
    mixtures: pd.DataFrame, rng: np.random.Generator = np.random.default_rng()
):
    return mixtures * rng.uniform(0.9, 1.1, size=(mixtures.shape))


def normalize_to_tp100k(geps: pd.DataFrame):
    return geps * 100_000 / geps.sum()
