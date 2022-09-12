# Only perturbing gene expression

Does cibersortx identify DEGs and their fold changes, when nothing else is different?

## Experiment 1: 100 DEGs in malignant cells

- Simulate bulk RNA-seq
  - group A: no change
  - group B: 100 genes with 2x expression in malignant cells
- Run cibersortx on this

### Implementation of generating data

1. determine fractions for samples
    - load estimated fractions for TCGA SKCM
    - limit to metastases
    - randomly sample 100 fractions without replacement
1. load scRNA-seq
    - exclude genes sparse in TCGA SKCM
    - exclude gene symbols not also in TCGA SKCM
1. determine genes to perturb
    - in malignant cells, figure out genes 
1. perturb gene expression
    - in simulated cell type GEP of malignant cells, multiply gene expression by 2
1. Generate simulated samples from cell type GEPs, fractions
    - add some noise
    - normalize
