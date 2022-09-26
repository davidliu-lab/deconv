# todo / tasks

- can deconv methods infer cell type-specific differential gene expression?
  - negative control: compare two iid generated (unperturbed) bulk RNA-seq cohorts
    - [ ] generate data
      - in `analysis/simulating_bulk_rnaseq/no_perturbations.py`
      - [x] use iid randomly sampled fraction vectors
      - [x] use iid randomly sampled cell type-specific GEPs
    - [ ] run cibersortx
      - [ ] provide true cell type proportions
    - generate volcano plots
      - [ ] at bulk level (Mann-Whitney of simulated bulk RNA-seq)
      - [ ] for malignant cell GEPs inferred by CIBERSORTx
  - experiment: compare with various fold change perturbations of gene expression in malignant cells
    - [x] generate perturbed data (but check results)
      - in `analysis/simulating_bulk_rnaseq/perturbed_malignant_expression.py`
      - [x] do scaling factors $2^n, n \in \{-3, -2, -1, 1, 2, 3\}$)
      - [x] use iid randomly sampled fraction vectors
      - [x] use iid randomly sampled cell type-specific GEPs
      - [x] save fractions
      - [x] save bulk_rnaseq
      - [x] save cell type-specific GEPs
      - [x] save perturbed genes (eg cell type, gene-specific fold changes)
      - [x] log the save paths for each scaling factor
    - [ ] run CIBERSORTx
      - [x] provide true fractions
        - [x] in helper library, make function for generating fraction file to provide CIBERSORTx
    -[ ] make array of volcano plots for each fold change
      - [ ] DGE in bulk RNA-seq
      - [ ] DGE in inferred malignant-specific expression
      - [ ] add malignant fraction to hover data

- improve DEG volcano plots
  - change FDR thresholds...
    - [ ] in stats results, add 0.1, 0.25 FDR thresholds
    - [ ] add lines for these to volcano plot
  - [ ] adjust y axis to be long enough to include FDR thresholds when no DEGs are found
  - [ ] make x axes symmetrical (`fig.update_xaxes(range=[-same, same])`)
  - make scatter + kde plots
    - see [plotly.figure_factory.create_2d_density](https://plotly.com/python/v3/density-plots/) and [px.density_contour](https://plotly.com/python/2d-histogram-contour/)
  - [ ] look into making a `go.ScatterGL` plot instead of `px.scatter`, or using `backend='webgl'`
  - [ ] make sure multiple hypothesis stats are computed _after_ any filters (e.g. limiting cell type)

- general feedback...
  - limit things that can change in each experiment
    - e.g., provide true fractions when running cell type-specific gene expression inference
    - want to know the impact of each step on performance
  - but be clinically plausible
    - sample things iid, not exactly the same (eg fraction composition, cell type GEPs)

- refactors
  - in DEG analysis
    - [ ] compute stats from a `pd.Series` with `"sample_group"` included as an index, not a DataFrame with two columns.
    - [ ] use "baseline" and "other"
  - [ ] use `bulk_rnaseq` in variable and file names
  - write dataframes with `upath.UPath of URI `str` (e.g. `pd.to_csv(path)`)
    - [ ] in `helpers.running_cibersortx.creating_input_files`
    - [ ] in `helpers.running_cibersortx.*.run_and_upload`
  - [ ] make `run_and_upload_from_dataframes` for other cibersortx endpoints
  - [ ] move `columns`, `cell_type_naming` to `data_io_and_formatting`
  - [x] in `analysis/cibersortx/perturbed_gene_expression/run_cibersortx.py` move `load_and_concatenate` functions to somewhere in `helpers`, because i'm reusing it elsewhere.
  - add caching

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
