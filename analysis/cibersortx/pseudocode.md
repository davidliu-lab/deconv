randomly sample 100 genes to perturb

1. select 50 sample fractions
2. for each sample
    2. create cell type GEPs from scRNA-seq
    3. create bulk rna-seq
2. select 100 genes to perturb
3. for all log2_fc values
    1. for all genes and samples, perturb gene expression
    2. create bulk rna-seq
    3. deconvolve with cibersortx
    4. perform DEG analysis vs control
        1. for bulk GEP
        2. for inferred malignant GEP
