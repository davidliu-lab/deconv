import numpy as np
import pandas as pd


def make_cohort(
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
    cell_types = list(sorted(sc_metadata["cell.type"].unique()))
    n_cell_types = len(cell_types)
    fraction_values = rng.dirichlet(alpha=(0.5,) * n_cell_types, size=(n_samples,))
    fractions = pd.DataFrame(fraction_values, index=samples, columns=cell_types)

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
    mixtures = pd.concat()

    return mixtures, fractions, cell_type_geps
