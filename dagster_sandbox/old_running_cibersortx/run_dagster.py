"""
a b
a - no gene perturbation, some malignant fraction
b - 100 genes perturbed, some malignant fraction
produce datasets...
- with specified malignant fraction
- with specified log2 fold change
- with specified seed

terminal assets:
for gep_perturbation in (True, False):
    for a_malignant_fraction in (None, 0.5, 0.7, 0.9):
        for b_malignant_fraction in (None, 0.5, 0.7, 0.9):
            data_a = load_data(gep_perturbation, a_malignant_fraction)
            data_b = load_data(gep_perturbation, b_malignant_fraction)
            data_all = pd.concat([data_a.bulk_rnaseq, data_b.bulk_rnaseq])
            fractions_all = pd.concat([data_a.fractions, data_b.fractions], axis=1)
            run_cibersortx(data_all, fractions_all)
"""

"""
would like to know...
- grid of compositional difference and gene expression difference


"""
