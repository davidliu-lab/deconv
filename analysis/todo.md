# todo / tasks

- improve analyses / incorporate feedback
    - DEG plots...
        -  [x] add line for FDR-adjusted significance
            - [x] in deg code, compute the threshold, store it with results: `signif_bh_adjusted`
            - [x] add line to plots: `fig.add_hline(y=signif_bh_threshold)`
        - [x] On y-axis, plot original p-values (not BH-corrected q-values)
    - only analyze malignant cell DEGs (otherwise too complicated)
        - [x] add plot of DEG limited to malignant cells
            - https://14963bb53b0d1b7a-dot-us-east4.notebooks.googleusercontent.com/http-server/deconv/analysis/_build/html/evaluating_cibersortx/identically_generated_cohorts/samples_like_tcga_skcm_mets.html#malignant-cells-only
- set up experiments:
    - inferring cell type DEGs

    - inferring differential cell type proportions, but with no DEGs
        - [ ]  make script - perturb relative abundance by 2x, 1/2x
    - inferring differential proportions _and_ gene expression
        - [ ]  make script for generating data, deconvolving
    - to do in cohort gen...
        - just do 100 samples total (50 vs 50)
        - apply gene filters
        - normalize to TPM
        - save out all data: (1) final bulk rna-seq, (2) fractions, (3) cell type geps, (4) perturbations


# cibersortx evaluations: 

- how well does cibersortx...
    - infer *cell type proportions* when no genes are diff expressed?
    - infer *DEGs* when cell type proportions are the same?
        - [ ] create & run script
            - [ ] perturb 100 genes in malignant cells by 2x
    - infer *DEGs* and *cell type proportions* when both are different?
        - [ ] make script to simulate data, deconvolve
            - [ ] perturb 100 genes in malignant cells by 2x
            - [ ] perturb cell type fraction of malignant cells by 2x
        - 
        - perturbing 100 genes in malignant cells
        - perturbing malignant cell fractions
