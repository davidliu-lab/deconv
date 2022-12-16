# project planning

## to-do

experiments
- `gs://liulab/run_everything/20221101_08h45m56s`
  - big
  - `select_100_genes_at_least_somewhat_expressed_in_malignant`
- `gs://liulab/run_everything/20221028_21h18m11s`
  - big
  - `select_100_genes`

others
- `gs://liulab/run_everything/20221021_22h18m15s`
  - big
  - `select_100_genes_at_least_somewhat_expressed_in_malignant`
  - bug!


```
# gs://liulab/run_everything/20221020_17h13m21s/ # big, select_100_genes
# gs://liulab/run_everything/20221020_18h54m44s/ # big, select_100_genes_at_least_somewhat_expressed_in_malignant
# gs://liulab/run_everything/20221021_18h43m00s/ # small
# gs://liulab/run_everything/20221021_20h17m28s/ # ???
# gs://liulab/run_everything/20221021_22h18m15s/ # big, select_100_genes_at_least_somewhat_expressed_in_malignant
# gs://liulab/run_everything/20221028_21h18m11s/ big, select_100_genes
```

evaluating gene expression inference with known cell type composition
- [x] create first draft of one big script that executes all steps
- [x] change which genes to perturb: select 100 at random from all
- [x] change which genes to perturb: select 100 at random from genes expressed in at least 10% of malignant cells
- [ ] rename output subdirectory to `cibersortx_hires_only`
- add other metrics quantifying accuracy (not just plots) for gene expression in malignant cells
  - [ ] pearson correlation
  - [ ] check CIBERSORTx paper for other views, metrics
- [ ] add metadata for genes to plots
  - [ ] sparsity overall
  - [ ] sparsity in malignant cells
  - [ ] differential expression in malignant cells?
- volcano plots for different fold changes
  - [ ] add back FDR lines for plotting volcanos
  - [x] change how plots are generated: pass a dictionary of dataframes, not paths
- [ ] add hires inference using CIBERSORT-estimated fractions

evaluating fraction inference
- [ ] add fraction inference, saving to `cibersortx_fractions`

comparing synthetic data vs TCGA SKCM
- distributions of gene expression
  - compare gene sparsity
    - [ ] distplot/histogram (synthetic and TCGA SKCM)
  - compare mean gene expression
    - [ ] correlation (synthetic vs TCGA SKCM)
    - [ ] scatter plot (synthetic vs TCGA SKCM)
    - [ ] distplot/histogram (synthetic and TCGA SKCM)
  - compare stddev gene expression
    - [ ] correlation
    - [ ] scatter plot (synthetic vs TCGA SKCM)
    - [ ] distplot/histogram (synthetic and TCGA SKCM)
  - [ ] dimensionality of the data, i.e., how much variance explained by k PCs of each dataset
  - PCA
    - [ ] PCA of gene expression, showing TCGA SKCM bulk RNA-seq samples and synthetic samples
    - [ ] on PCA plot, maybe also genes with high loadings
    - [ ] for first PC, analysis of which genes contribute most to the PC (e.g. sparsity)
- distributions of samples
  - [ ] inter-group correlaton
  - [ ] intra-group correlation

other
- [ ] set up new VM
- refactors
  - [ ] use `DataFrame.squeeze` instead of `DataFrame.stack` (but does `dask.dataframe` implement `squeeze`?)
  - make `run_and_upload_from_dataframes` for other cibersortx endpoints
    - [x] `hires_only`
    - `fractions`
      - [x] started
      - [ ] test
    - [ ] `hires_and_fractions`
  - write dataframes with `upath.UPath` of URI (e.g. `pd.to_csv(path)`)
    - [ ] in `helpers.running_cibersortx.creating_input_files`
    - [ ] in `helpers.running_cibersortx.*.run_and_upload`
- add runner functions for other deconv methods
  - [ ] BayesPrism
  - [ ] CODEFACS
- fix table of contents thing on the right
  - [ ] upgrade sphinx book theme? check github issue


## analysis outlines

- can deconv methods infer cell type-specific differential gene expression?
  - negative control: compare two iid generated (unperturbed) bulk RNA-seq cohorts
  - experiment: compare with various fold change perturbations of gene expression in malignant cells

## general feedback

- limit things that can change in each experiment
  - e.g., provide true fractions when running cell type-specific gene expression inference
  - want to know the impact of each step on performance
- be clinically plausible
  - sample things iid, not exactly the same (eg fraction composition, cell type GEPs)

## previous

- generating synthetic data
  - control data
    - [x] generate data
      - in `analysis/simulating_bulk_rnaseq/no_perturbations.py`
      - [x] use iid randomly sampled fraction vectors
      - [x] use iid randomly sampled cell type-specific GEPs
    - [ ] run cibersortx
      - [x] provide true cell type proportions
    - generate volcano plots
      - [ ] at bulk level (Mann-Whitney of simulated bulk RNA-seq)
      - [ ] for malignant cell GEPs inferred by CIBERSORTx
  - experimental data
    - [x] generate perturbed data (but check results)
      - in `analysis/simulating_bulk_rnaseq/perturbed_malignant_expression.py`
      - [ ] do scaling factors $2^n, n \in \{-3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3\}$
    - [ ] run CIBERSORTx
      - [x] provide true fractions
        - [x] in helper library, make function for generating fraction file to provide CIBERSORTx
          - [x] make array of volcano plots for each fold change
      - [ ] DGE in bulk RNA-seq
      - [ ] DGE in inferred malignant-specific expression
      - [ ] add malignant fraction to hover data

- improve DEG volcano plots
  - change FDR thresholds...
    - [x] in stats results, add 0.1, 0.25 FDR thresholds
    - [x] add lines for these to volcano plot
  - [x] adjust y axis to be long enough to include FDR thresholds when no DEGs are found
  - [x] make x axes symmetrical (`fig.update_xaxes(range=[-same, same])`)
  - make scatter + kde plots
    - see [plotly.figure_factory.create_2d_density](https://plotly.com/python/v3/density-plots/) and [px.density_contour](https://plotly.com/python/2d-histogram-contour/)
