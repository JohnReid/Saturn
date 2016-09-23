#!/usr/bin/env Rscript
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                 # Partition
#SBATCH -J PREDSMOOTH                   # Name of the job
#SBATCH -A MRC-BSU-SL2                  # Which project should be charged
#SBATCH --nodes=1                       # How many whole nodes should be allocated?
#SBATCH --ntasks=1                      # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --cpus-per-task=1               # How many CPUs will be used per task
#SBATCH --mem=7680                      # How many MB each node is allocated
#SBATCH --time=06:00:00                 # How much wallclock time will be required?
#SBATCH -o pred-smooth-%j.out           # stdout
#SBATCH -e pred-smooth-%j.out           # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
                                        # interrupted by node failure or system downtime):
"Usage: predictions-smooth.R [options] IN OUT

Options:
  -l --length-scale=LENGTHSCALE  Use LENGTHSCALE [default: 50]
  --logodds                      Smooth log-odds instead of probabilities [default: FALSE]
  --width=MAXWIDTH               Width of kernel in units of regions [default: 20]" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Parse options
#
# .args <- "--logodds ../Data/Predictions/predictions.xgboost.chrfold.EGR1.H1-hESC.DNase_Known_KnownWell_DREME_DREMEWell.tsv ../slurm/smoothed-predictions.tsv"
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
print(opts)
in.path <- opts$IN
out.path <- opts$OUT
log.transform <- as.logical(opts[['logodds']])
max.width <- as.integer(opts[['width']])
length.scale <- as.numeric(opts[['length-scale']])


#
# Load predictions
#
preds <- data.table::fread(
  in.path,
  col.names = c('chrom', 'start', 'end', 'prediction', 'bound'))
N <- nrow(preds)
message('# predictions: ', N)


#' Distance between two rows
#'
row.similarity <- function(i, j) {
  ifelse(
    preds$chrom[i] == preds$chrom[j],
    exp(-((preds$start[i] - preds$start[j])/length.scale)**2/2),
    0)
}


#
# Create sparse symmetric banded similarity matrix
#
# Calculate non-zero indices and values
i <- as.vector(sapply(1:N, function(i) rep(i, max.width)))
j <- rep(1:max.width, N) + i - 1
x <- row.similarity(i, j)
#
# Discard indices outside matrix
keep <- j <= N
i <- i[keep]
j <- j[keep]
x <- x[keep]
#
# Create matrix
K <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(N, N), symmetric = TRUE)
message('Smoothing matrix size: ', object.size(K))
# K[1:12, 1:12]


#
# Normalise smoothing matrix by rows
#
message('Normalising smoothing matrix')
row.sums <- Matrix::rowSums(K)
K.norm <- Matrix::Diagonal(x = 1 / row.sums) %*% K
sums.norm <- Matrix::rowSums(K.norm)
stopifnot(all(abs(1 - sums.norm) < 1e-12))


#
# Smooth predictions
#
message('Smoothing predictions')
logistic <- function(x) 1/(1+exp(-x))
logit <- function(p) log(p/(1-p))
predictions <- preds$prediction
if (log.transform) {
  predictions <- logit(predictions)
}
predictions.smoothed <- K.norm %*% predictions
if (log.transform) {
  predictions.smoothed <- logistic(predictions.smoothed)
}


#
# Write predictions
#
message('Writing predictions')
preds$prediction <- as.vector(predictions.smoothed)
data.table::fwrite(preds, out.path, sep = '\t', col.names = FALSE)
