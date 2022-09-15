# todo / tasks

- can deconv methods infer cell type-specific differential gene expression?
  - control: compare two unperturbed generated bulk RNA-seq cohorts
    - [x] generate data
      - `gs://liulab/data/simulated/50_samples_no_perturbations/2022-09-13_16:02:18/`
      - `gs://liulab/data/simulated/50_samples_no_perturbations/2022-09-13_21:37:53/`
      - [ ] check if these have identical cell type fractions
    - [x] run cibersortx
      - [x] provide true cell type proportions
    - generate volcano plots
      - [x] at bulk level (Mann-Whitney of simulated bulk RNA-seq)
      - [x] for malignant cell GEPs inferred by CIBERSORTx
  - experiment: compare (1) unperturbed generated bulk RNA-seq vs (2) simulated bulk RNA-seq with 100 genes 2x perturbed in malignant cells
    - [ ] generate perturbed data
      - [x] fixed bug in perturbing cell type GEPs
      - [ ] confirm that cell type fractions are the same as in the unperturbed data
    - [ ] run CIBERSORTx
      - [x] provide true fractions
        - [x] in helper library, make function for generating fraction file to provide CIBERSORTx
    - generate volano plots for
      - [ ] DGE in bulk RNA-seq
      - [ ] DGE in inferred malignant-specific expression
- does accuracy vary by fold change magnitudes?
  - experiment: comparisons with other fold change perturbations (e.g. $(0.25, 0.5, 2, 4)$)
    - [ ] generate simulated data for each fold change
      - [ ] make sure fold change is documented in path or saved in data output
    - [ ] make array of plots for each fold change, comparing bulk DEG vs inferred malignant-specific DEG

- improve DEG volcano plots
  - FDR thresholds...
    - [ ] add thresholds to metadata (e.g. fdr $\in (0.05, 0.1, 0.25, 0.5)$)
    - [ ] adjust y axes of volcano plots to include thresholds
  - [ ] make sure multiple hypothesis stats are computed _after_ any filters (e.g. limiting cell type)
  - make scatter + kde plots
    - see [plotly.figure_factory.create_2d_density](https://plotly.com/python/v3/density-plots/) and [px.density_contour](https://plotly.com/python/2d-histogram-contour/)

- general feedback...
  - limit things that can change in each experiment
    - e.g., provide true fractions when running cell type-specific gene expression inference
    - want to know the impact of each step on performance
  - no specific suggestion for sampling fractions (maybe random combinations of fractions, dirichlet process, etc.)

- refactors
  - [ ] use `bulk_rnaseq` in variable and file names
  - [ ] change timestamps to be more file system friendly (e.g. `2021-01-01_12-34-56`)
  - write dataframes with `cloudpathlib.AnyPath` instead of URI `str` (e.g. `pd.to_csv(path)`)
    - [ ] in `helpers.running_cibersortx.copying_to_gcs.copy_local_directory_to_gcs`
    - [ ] in `helpers.running_cibersortx.creating_input_files`
    - [ ] in `helpers.running_cibersortx.*.run_and_upload`
  - [ ] make `run_and_upload_from_dataframes` for other cibersortx endpoints
  - [ ] move `columns`, `cell_type_naming` to `data_io_and_formatting`
  - [x] in `analysis/evaluating_cibersortx/perturbed_gene_expression/run_cibersortx.py` move `load_and_concatenate` functions to somewhere in `helpers`, because i'm reusing it elsewhere.
  - [ ] look into cloudpathlib alternatives
    - https://github.com/Advestis/transparentpath#usage 
    - https://github.com/fsspec/universal_pathlib

- add back evaluation of simulated data
  - means and stddevs of gene expression
  - PCA of gene expression
    - for first PC, analysis of which genes contribute most to the PC (e.g. sparsity)
  - correlation of gene expression between simulated samples, real TCGA SKCM samples
  - intra-cohort variation
  - parameter tuning

- eventually
  - add runner functions for other deconv methods (BayesPrism, CODEFACS, ...)

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
