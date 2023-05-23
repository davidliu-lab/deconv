suppressWarnings(library(BayesPrism))

print("Loading data from ./BayesPrism/tutorial.dat/tutorial.gbm.rdata")
load("./BayesPrism/tutorial.dat/tutorial.gbm.rdata")
ls()

# run BayesPrism
myPrism <- new.prism(
  reference = sc.dat,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = "tumor",
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores = 45)

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
