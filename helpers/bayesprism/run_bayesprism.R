suppressWarnings(library(BayesPrism))
library("logger")
library("optparse")
library("arrow")
library("tidyverse")

# define the option parser
option_list <- list(
  make_option("--sc_rnaseq_uri"),
  make_option("--sc_rnaseq_cell_types_uri"),
  make_option("--sc_rnaseq_cell_states_uri"),
  make_option("--bulk_rnaseq_uri"),
)
# parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

read_parquet <- function(path_uri) {
  logger::log_info("Loading data from ", path_uri)
  data <- arrow::read_parquet(file = path_uri)
  return(data)
}

sc_rnaseq <- read_parquet(opt$sc_rnaseq_uri)
sc_rnaseq_cell_types <- read_parquet(opt$sc_rnaseq_cell_types_uri)
sc_rnaseq_cell_states <- read_parquet(opt$sc_rnaseq_cell_states_uri) 
bulk_rnaseq <- read_parquet(opt$bulk_rnaseq_uri) %>% 
  column_to_rownames(var = "gene_symbol") %>% 
  t()

for (df in c(sc_rnaseq, sc_rnaseq_cell_types, sc_rnaseq_cell_states, bulk_rnaseq)) {
  logger::log_info("Matrix (hopefully): ")
  typeof(df)
  dim(df)
  head(rownames(df))
  head(colnames(df))
}

suppressWarnings(library("BayesPrism"))

myPrism <- BayesPrism::new.prism(
  reference = sc_rnaseq,
  mixture = bulk_rnaseq,
  input.type = "count.matrix",
  cell.type.labels = sc_rnaseq_cell_types,
  cell.state.labels = sc_rnaseq_cell_states,
  key = "Malignant",
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)

bp.res <- BayesPrism::run.prism(prism = myPrism, n.cores = 45)

# print bp.res to console
bp.res

# extract posterior mean of cell type fraction theta
theta <- BayesPrism::get.fraction(
  bp = bp.res,
  which.theta = "final",
  state.or.type = "type"
)

head(theta)
# write to parquet file
arrow::write_parquet(
  theta,
  path = "theta.parquet"
)

# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)
arrow::write_parquet(
  theta.cv,
  path = "theta_cv.parquet"
)

Z.tumor <- get.exp(
  bp = bp.res,
  state.or.type = "type",
  cell.name = "Malignant"
)

# print first five rows of the transpose of Z.tumor
head(t(Z.tumor)[1:5, ])

arrow::write_parquet(
  t(Z.tumor),
  path = "Z_tumor.parquet"
)
