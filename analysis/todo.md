# todo / tasks

- can deconv methods infer cell type-specific differential gene expression? for what fold change magnitudes?
  - experiment: compare (1) unperturbed generated bulk RNA-seq vs (2) simulated bulk RNA-seq with 100 genes 2x perturbed in malignant cells
    - [ ] run CIBERSORTx
      - [x] provide true fractions
        - [x] in helper library, make function for generating fraction file to provide CIBERSORTx
    - generate volano plots for
      - [ ] DGE in bulk RNA-seq
      - [ ] DGE in inferred malignant-specific expression
  - experiment: comparisons with other fold change perturbations (e.g. $(0.25, 0.5, 2, 4)$)
    - [ ] generate simulated data for each fold change
      - [ ] make sure fold change is documented in path or saved in data output
    - [ ] make array of plots for each fold change, comparing bulk DEG vs inferred malignant-specific DEG

- improve DEG plots
  - FDR thresholds...
    - [ ] add thresholds to metadata (e.g. fdr $\in (0.05, 0.1, 0.25, 0.5)$)
    - [ ] adjust y axes of volcano plots to include thresholds
  - [ ] make sure multiple hypothesis stats are computed _after_ any filters (e.g. limiting cell type)

- general feedback...
  - limit variable things in each experiment.
    - e.g., provide true fractions when running cell type-specific gene expression inference
    - want to know the impact of each step on performance
  - no specific suggestion for sampling fractions (maybe random combinations of fractions, dirichlet process, etc.)

- nice refactors
  - [ ] in `helpers.running_cibersortx.copying_to_gcs.copy_local_directory_to_gcs`, use `cloudpathlib.CloudPath` instead of URI `str` for target
  - [ ] in `analysis/evaluating_cibersortx/perturbed_gene_expression/run_cibersortx.py` move `load_and_concatenate` functions to somewhere in `helpers`, because i'm reusing it elsewhere.

# cibersortx evaluations: 

- how well does cibersortx...
    - infer *cell type proportions* when no genes are diff expressed?
    - infer *DEGs* when cell type proportions are the same?
    - infer *DEGs* and *cell type proportions* when both are different?
        - perturbing 100 genes in malignant cells
        - perturbing malignant cell fractions

# previously

## late august

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
        - [x]  make script - perturb relative abundance by 2x, 1/2x
    - inferring differential proportions _and_ gene expression
        - [ ]  make script for generating data, deconvolving
    - to do in cohort gen...
        - just do 100 samples total (50 vs 50)
        - apply gene filters
        - normalize to TPM
        - save out all data: (1) final bulk rna-seq, (2) fractions, (3) cell type geps, (4) perturbations
