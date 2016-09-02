#!/usr/bin/env Rscript
"Usage: predictions-smooth.R IN OUT" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Parse options
#
.args <- "../Data/Predictions/predictions.xgboost.chrfold.EGR1.H1-hESC.DNase_Known_KnownWell_DREME_DREMEWell.tsv smoothed-predictions.tsv"
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
in.path <- opts$IN
out.path <- opts$OUT


log.transform <- TRUE
MAX.WIDTH <- 5  # Maximum width of our kernel for smoothing
L <- 200        # Length scale


#
# Load predictions
#
preds <- data.table::fread(
  in.path,
  col.names = c('chrom', 'start', 'end', 'prediction', 'bound'))
head(preds)
sapply(preds, class)
N <- nrow(preds)
message('# predictions: ', N)


#' Distance between two rows
#'
row.similarity <- function(i, j) {
  ifelse(
    preds$chrom[i] == preds$chrom[j],
    exp(-((preds$start[i] - preds$start[j])/L)**2/2),
    0)
}


#
# Create sparse symmetric banded similarity matrix
#
# Calculate non-zero indices and values
i <- as.vector(sapply(1:N, function(i) rep(i, MAX.SEP)))
j <- rep(1:MAX.SEP, N) + i - 1
x <- row.similarity(i, j)
#
# Discard indices outside matrix
keep <- j <= N
sum(! keep)
i <- i[keep]
j <- j[keep]
x <- x[keep]
#
# Create matrix
K <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(N, N), symmetric = TRUE)
object.size(K)
# K[1:12, 1:12]


#
# Normalise smoothing matrix by rows
#
row.sums <- Matrix::rowSums(K)
K.norm <- Matrix::Diagonal(x = 1 / row.sums) %*% K
sums.norm <- Matrix::rowSums(K.norm)
stopifnot(all(abs(1 - sums.norm) < 1e-12))


#
# Smooth predictions
#
logistic <- function(x) 1/(1+exp(-x))
logit <- function(p) log(p/(1-p))
predictions <- preds$prediction
if (log.transform) {
  predictions <- logit(predictions)
}
predictions.smoothed <- K.norm %*% predictions
if (log.transform) {
  predictions <- logistic(predictions)
}


#
# Write predictions
#
preds$prediction <- predictions
data.table::fwrite(preds, out.path, sep = '\t', col.names = FALSE)
