import numpy as np
import pandas as pd


def make_a_cell_type_gep_matrix(sc_data, sc_metadata, n_cells_per_gep, rng):
    # for one in silico mixture...
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
    cell_type_geps *= 100_000 / cell_type_geps.sum()
    return cell_type_geps


def make_mixtures(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    n_samples: int,
    n_cells_per_gep: int,
    rng: np.random.RandomState,
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
    cell_type_geps = {
        sample: make_a_cell_type_gep_matrix(sc_data, sc_metadata, n_cells_per_gep, rng)
        for sample in samples
    }
    mixtures = pd.concat(
        {sample: cell_type_geps[sample] @ fractions.loc[sample] for sample in samples},
        axis=1,
    )
    mixtures *= rng.uniform(0.9, 1.1, size=(mixtures.shape))
    mixtures *= 100_000 / mixtures.sum()
    return mixtures, fractions, cell_type_geps
