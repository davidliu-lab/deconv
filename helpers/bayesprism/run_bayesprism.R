library("logger")
library("optparse")

# define the option parser
option_list <- list(
  make_option(
    opt_str = "--reference_sc_rnaseq_uri",
    dest = "reference_sc_rnaseq_uri"
  ),
  make_option(
    opt_str = "--reference_sc_rnaseq_annotations_uri",
    dest = "reference_sc_rnaseq_annotations_uri",
  ),
  make_option(
    opt_str = "--bulk_rnaseq_uri",
    dest = "bulk_rnaseq_uri"
  )
)
# parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

read_data <- function(path_uri) {
  # log message about reading from path_uri with level DEBUG
  log_info("Loading data from ", path_uri)
  library("arrow")
  library("fs")
  library("readr")

  # Check the file extension to determine the file type
  thing <- path(fs::path(path_uri))
  log_info("thing: ", thing)
  file_extension <- tools::file_ext(path_uri)
  if (file_extension == "csv") {
    log_info("Read the data from a CSV file")
    # data <- read.csv()
  } else if (file_extension == "parquet") {
    log_info("Read the data from a Parquet file")
    data <- arrow::read_parquet(path(fs::path(path_uri)))
  } else {
    log_info("Handle unsupported file types")
    stop(paste0("Unsupported file type: ", file_extension))
  }
  return(data)
}
load("./BayesPrism/tutorial.dat/tutorial.gbm.rdata")
ls()

# bk.dat: The sample-by-gene raw count matrix of bulk RNA-seq expression.
# rownames are bulk sample IDs, while colnames are gene names/IDs.
dim(bk.dat)
#> [1]   169 60483
head(rownames(bk.dat))
#> [1] "TCGA-06-2563-01A-01R-1849-01" "TCGA-06-0749-01A-01R-1849-01" "TCGA-06-5418-01A-01R-1849-01" "TCGA-06-0211-01B-01R-1849-01" "TCGA-19-2625-01A-01R-1850-01" "TCGA-19-4065-02A-11R-2005-01"
head(colnames(bk.dat))
#> [1] "ENSG00000000003.13" "ENSG00000000005.5"  "ENSG00000000419.11" "ENSG00000000457.12" "ENSG00000000460.15" "ENSG00000000938.11"

bulk_rnaseq_matrix <- read_data(opt$bulk_rnaseq_uri)
dim(bulk_rnaseq_matrix)
head(rownames(bulk_rnaseq_matrix))
head(colnames(bulk_rnaseq_matrix))

# sc.dat: The cell-by-gene raw count matrix of bulk RNA-seq expression.
# rownames are bulk cell IDs, while colnames are gene names/IDs.
dim(sc.dat)
#> [1] 23793 60294
head(rownames(sc.dat))
#> [1] "PJ016.V3" "PJ016.V4" "PJ016.V5" "PJ016.V6" "PJ016.V7" "PJ016.V8"
head(colnames(sc.dat))

reference_sc_rnaseq_matrix <- read_data(opt$reference_sc_rnaseq_uri)
dim(reference_sc_rnaseq_matrix)
head(rownames(reference_sc_rnaseq_matrix))
head(colnames(reference_sc_rnaseq_matrix))

sort(table(cell.type.labels))

sort(table(cell.state.labels))

table(cbind.data.frame(cell.state.labels, cell.type.labels))

# exit status 0
quit(status = 0)

data_reference_scrnaseq <- read_data(opt$_)
dim(data_reference_scrnaseq)
data_reference_scrnaseq_annotations <- read_data(opt$_)
data_bulkrnaseq <- read_data(opt$bulk_rnaseq_uri)

# print data shape
dim(data_reference_scrnaseq)
dim(data_reference_scrnaseq_annotations)
dim(data_bulkrnaseq)


plot.cor.phi(
  input = sc.dat,
  input.labels = cell.state.labels,
  title = "cell state correlation",
  # specify pdf.prefix if need to output to pdf
  # pdf.prefix="gbm.cor.cs",
  cexRow = 0.2, cexCol = 0.2,
  margins = c(2, 2)
)

plot.cor.phi(
  input = sc.dat,
  input.labels = cell.type.labels,
  title = "cell type correlation",
  # specify pdf.prefix if need to output to pdf
  # pdf.prefix="gbm.cor.ct",
  cexRow = 0.5, cexCol = 0.5,
)

sc.stat <- plot.scRNA.outlier(
  input = sc.dat, # make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels = cell.type.labels,
  species = "hs", # currently only human(hs) and mouse(mm) annotations are supported
  return.raw = TRUE # return the data used for plotting.
  # pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

head(sc.stat)

bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat, # make sure the colnames are gene symbol or ENSMEBL ID
  sc.input = sc.dat, # make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels = cell.type.labels,
  species = "hs", # currently only human(hs) and mouse(mm) annotations are supported
  return.raw = TRUE
  # pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

head(bk.stat)

sc.dat.filtered <- cleanup.genes(
  input = sc.dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)


dim(sc.dat.filtered)


# note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc(
  sc.input = sc.dat.filtered,
  bulk.input = bk.dat
  # pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
#> EMSEMBLE IDs detected.

sc.dat.filtered.pc <- select.gene.type(
  sc.dat.filtered,
  gene.type = "protein_coding"
)

diff.exp.stat <- get.exp.stat(
  sc.dat = sc.dat[, colSums(sc.dat > 0) > 3], # filter genes to reduce memory use
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  psuedo.count = 0.1, # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
  cell.count.cutoff = 50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
  n.cores = 1 # number of threads
)

sc.dat.filtered.pc.sig <- select.marker(
  sc.dat = sc.dat.filtered.pc,
  stat = diff.exp.stat,
  pval.max = 0.01,
  lfc.min = 0.1
)

dim(sc.dat.filtered.pc.sig)

library("BayesPrism")
# print the functions in BayesPrism
ls("package:BayesPrism")

# exit the program with error code 0
quit(status = 0)

myPrism <- new.prism(
  reference = sc.dat.filtered.pc,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = "tumor",
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores = 25)

# print bp.res to console
bp.res

# extract posterior mean of cell type fraction theta
theta <- get.fraction(
  bp = bp.res,
  which.theta = "final",
  state.or.type = "type"
)

head(theta)

# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)

Z.tumor <- get.exp(
  bp = bp.res,
  state.or.type = "type",
  cell.name = "tumor"
)

# print first five rows of the transpose of Z.tumor
head(t(Z.tumor)[1:5, ])
