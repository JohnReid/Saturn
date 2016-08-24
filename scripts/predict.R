#!/usr/bin/env Rscript
#
# Simple prediction script. Analyses data from other cell types for a TF to
# predict binding in the VALIDATIONCELL type.
#
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                 # Partition
#SBATCH -J PREDICT                      # Name of the job
#SBATCH -A MRC-BSU-SL2                  # Which project should be charged
#SBATCH --nodes=1                       # How many whole nodes should be allocated?
#SBATCH --ntasks=1                      # How many tasks will there be in total? (<= nodes*16)
#SBATCH --cpus-per-task=16              # How many CPUs per task
#SBATCH --mem=61440                     # How many MB each node is allocated
#SBATCH --time=12:00:00                 # How much wallclock time will be required?
#SBATCH -o predict/predict-%j.out       # stdout
#SBATCH -e predict/predict-%j.out       # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
                                        # interrupted by node failure or system downtime):
Sys.setenv(OMP_NUM_THREADS=16)          # Set OpenMP number of threads
"Usage:
predict.R [options] [--features=NAME]... TF VALIDATIONCELL

Options:
  --method=METHOD        Use METHOD [default: glmnet]
  -f --features=NAME     Use NAME features
  --use-zero-dnase       Fit regions with zero DNase levels [default: FALSE]
  --max-boosting=MAXIMUM MAXIMUM number of boost rounds for xgboost method [default: 300]
  -s --sample=PROP       Subsample training regions to with proportion PROP [default: 1]" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load libraries
devtools::load_all()
library(Saturn)


#
# Parse options
#
# .args <- "--motif=Known ARID3A K562"
# .args <- "--motif=Known GATA3 A549"
# .args <- "-f DNase -f KnownMotifs -s .01 GABPA SK-N-SH"
# .args <- "-f DNase -f DREMEWell ATF2 GM12878"
# .args <- "--method=xgboost -f DNase -f DREMEWell ATF2 GM12878"
# Use dummy arguments if they exist otherwise use command line arguments
if (! exists(".args")) .args <- commandArgs(TRUE)
message('Arguments: ', .args)
opts <- docopt::docopt(doc, args = .args)
print(opts)
tf <- factor(opts$TF, tf.levels)
if (is.na(tf)) stop('Unknown TF specified.')
cell.valid <- factor(opts$VALIDATIONCELL, levels = cell.levels)
if (is.na(cell.valid)) stop('Unknown validation cell specified.')
method <- opts$method
feat.names <- opts$features
sample.prop <- as.numeric(opts$sample)
max.boost.rounds <- as.integer(opts[['max-boosting']])
use.zero.dnase <- as.logical(opts[['use-zero-dnase']])
message('Features: ', toString(feat.names))
if (! 'DNase' %in% feat.names) {
  message('WARNING: No DNase feature included in arguments!!!')
}
message('Sample proportion: ', toString(sample.prop))
message('Use zero DNase: ', toString(use.zero.dnase))


#
# Construct output filenames
#
feat.tags <- do.call(stringr::str_c, c(feat.names, list(sep = "_")))
fit.id <- stringr::str_c(method, '.', as.character(tf), '.', as.character(cell.valid), '.', feat.tags)
predictions.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('predictions.', fit.id, '.tsv'))


#
# Load prediction package
#
if ('xgboost' == method) {
  library(xgboost)
  fit.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('fit.', fit.id, '.xgb'))
} else if ('glmnet' == method) {
  library(glmnet)
  fit.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('fit.', fit.id, '.rds'))
} else {
  stop(stringr::str_c('Unknown method: ', method))
}


#
# Check if already run
#
if (file.exists(fit.path) && file.exists(predictions.path)) {
  message('Predictions and fit output files already exist: ', predictions.path)
  quit(save = "no")
}


#
# What sort of validation cell is it?
#
tf.cells <- tfs %>% filter(TF == tf)
valid.split <- (tf.cells %>% filter(cell == cell.valid))$split
message('Validation cell split: ', as.character(valid.split))


#
# Which cells do we have data for?
#
if ('ladder' == valid.split) {
  tf.cells <- tf.cells %>% filter(split != 'submit')
} else if ('submit' == valid.split) {
  tf.cells <- tf.cells %>% filter(split != 'ladder')
} else if ('train' == valid.split) {
  tf.cells <- tf.cells %>% filter(split != 'ladder', split != 'submit')
} else {
  stop('Validation cell split must be one of "train", "ladder" or "submit"')
}
print(tf.cells)
cell.all <- tf.cells$cell
if (! cell.valid %in% cell.all) stop('We have no data for the validation cell type and this TF.')
cell.train <- filter(tf.cells, cell != cell.valid)$cell
if (! length(cell.train)) stop('We have no training data for this cell')


#
# Fix training and validation chromosomes
#
if ('submit' == valid.split) {
  chrs.train <- factor(stringr::str_c('chr', c(2:7, 9:20, 22)), levels = chr.levels)
  chrs.valid <- factor(chr.levels, levels = chr.levels)
} else if ('ladder' == valid.split) {
  chrs.train <- factor(stringr::str_c('chr', c(2:7, 9:20, 22)), levels = chr.levels)
  chrs.valid <- factor(stringr::str_c('chr', c(1, 8, 21)), levels = chr.levels)
} else if ('train' == valid.split) {
  chrs.train <- factor(stringr::str_c('chr', c(3:6, 9:19, 22)), levels = chr.levels)
  chrs.valid <- factor(c('chr2', 'chr7', 'chr20'), levels = chr.levels)
} else {
  stop('Validation cell split must be one of "train", "ladder" or "submit"')
}
message('Chromosomes for training:   ',
        do.call(stringr::str_c, c(as.character(chrs.train), list(sep = ', '))))
message('Chromosomes for validation: ',
        do.call(stringr::str_c, c(as.character(chrs.valid), list(sep = ', '))))


#
# Load the binding data and expand to full test range
#
chip <- load.chip.features(tf)
chip.cols <- colnames(chip)
chip <- do.call(S4Vectors::DataFrame, lapply(chip, train.factor.Rle.to.test))
colnames(chip) <- chip.cols
message('ChIP data size: ', object.size(chip))


#
# Work out which regions to use in each cell, ignoring
# that are not in our training or validation sets.
#
train.idxs <- regions.test$chrom %in% chrs.train
valid.idxs <- regions.test$chrom %in% chrs.valid
#' Get test indices for the particular cell be it a training or a validation cell
#'
regions.for.cell <- function(cell) {
  if (cell %in% cell.train) {
      train.idxs
  } else {
      stopifnot(cell %in% cell.valid)
      valid.idxs
  }
}


#' Get the proportion bound
#'
prop.bound <- function(binding) sum('B' == binding) / length(binding)


#
# Load features
#
load.cell.data <- function(
  cell,
  .use.zero.dnase = FALSE,
  .sample.prop = sample.prop,
  .remove.ambiguous = TRUE)
{
  #
  # Get the binding status
  response <- chip[[as.character(cell)]]
  # Work out which regions to keep
  keep <- regions.for.cell(cell)
  if (! is.null(response) & .remove.ambiguous) {
    # Work out which regions to keep (ignore ambiguously bound regions)
    keep <- keep & ('A' != response)
  }
  response <- response[keep]
  message(cell, ': Proportion bound: ', prop.bound(response))
  #
  # Load and cbind the features
  mat <-
    do.call(
      cbind,
      lapply(feat.names, function(feat.name) load.feature(feat.name, tf, cell)[as.vector(keep), , drop = FALSE]))
  #
  # Remove those with zero DNase levels if we don't want them
  if (! .use.zero.dnase && 'DNase' %in% colnames(mat)) {
    dnase.non.zeros <- 0 != mat[,'DNase']
    message(cell, ': Removing regions with no DNase. ',
            'Reducing regions from ', length(dnase.non.zeros),
            ' to ', sum(dnase.non.zeros))
    mat <- mat[dnase.non.zeros,]
    response <- response[dnase.non.zeros]
    message(cell, ': Proportion bound: ', prop.bound(response))
  }
  #
  # Subsample the regions if requested
  if (.sample.prop < 1) {
    .sample <- sample.int(nrow(mat), size = as.integer(.sample.prop * nrow(mat)))
    mat <- mat[.sample,]
    response <- response[.sample]
  }
  message(cell, ': Has ', nrow(mat), ' regions')
  list(
    features = mat,
    response = factor(as.factor(response), levels = c('U', 'B')))
}
message('Creating training data')
train.data <- lapply(cell.train, load.cell.data)
train.feat <- do.call(rbind, lapply(train.data, function(d) d$features))
train.resp <- do.call(c,     lapply(train.data, function(d) d$response))
message('Training matrix size: ', object.size(train.feat))
message('# regions : ', nrow(train.feat))
message('# features: ', ncol(train.feat))


#
# Do cross-validation on number of boosting rounds
#
xgboost.fit <- function(
  data,
  nround = 300,
  nfold = 5,
  early.stop.round = 50
) {
  #
  # Perform cross-validation to choose number of rounds
  param <- list(silent = 1, objective = 'binary:logistic')
  message('Cross-validating')
  cv.time <- system.time(
    cvresult <- xgb.cv(
      params = param,
      data = data,
      nround = nround,
      nfold = nfold,
      early_stopping_rounds = early.stop.round,
      maximize = TRUE,
      eval_metric = 'map'))
  print(cv.time)
  nround.best <- which.max(cvresult$test.map.mean)
  message('CV best number of rounds: ', nround.best)
  if (nround.best == nround) {
    message('WARNING: Potential underfitting: CV best round is last round: ', nround.best)
  }
  #
  # Fit the boosted tree
  message('Fitting')
  fit.time <- system.time(
    fit <- xgboost(
      data = data,
      params = param,
      nround = nround.best,
      eval_metric = 'map'))
  print(fit.time)
  fit
}


#
# Fit model
#
if ('xgboost' == method) {
  #
  # Fit with xgboost
  #
  message('Fitting model with xgboost')
  dtrain <- xgb.DMatrix(train.feat, label = train.resp - 1)
  fit <- xgboost.fit(data = dtrain)
  xgb.save(fit, fit.path)
} else if ('glmnet' == method) {
  #
  # Fit a logistic regression
  #
  message('Fitting model with glmnet')
  system.time(cvfit <- cv.glmnet(train.feat, train.resp, family = 'binomial'))
  summary(cvfit)
  # plot(cvfit)
  coef(cvfit, s = "lambda.min")
  # Save fit
  message('Saving fit: ', fit.path)
  saveRDS(cvfit, fit.path)
}
rm(train.feat)  # No longer needed


#
# Create validation matrix
#
message('Creating validation data')
valid.feat <- load.cell.data(cell.valid, .use.zero.dnase = TRUE, .sample.prop = 1, .remove.ambiguous = FALSE)$features
message('Validation data size : ', object.size(valid.feat))
message('# validation regions : ', nrow(valid.feat))


#
# Make predictions on validation data
#
message('Making predictions')
if ('xgboost' == method) {
  class(fit)
  system.time(predictions <- predict(fit, valid.feat))
  # system.time(predictions <- predict(fit, as.matrix(valid.feat[1:10000,])))
  # system.time(predictions <- predict(fit, as.matrix(valid.feat)))
  stopifnot(length(predictions) == sum(regions.for.cell(cell.valid)))
} else if ('glmnet' == method) {
  system.time(predictions <- logit.inv(predict(cvfit, valid.feat, s = "lambda.min")[,1]))
}
rm(valid.feat)  # No longer needed


#
# Write predictions in format and order needed for submission
#
valid.keep <- regions.for.cell(cell.valid)
chrom <- regions.test$chrom[valid.keep]
start <- regions.test$start[as.vector(valid.keep)]
misordered.levels <- sort(chr.levels)  # Use alphabetical order that DREAM uses
out <-
  data.frame(
    chrom = Rle(factor(runValue(chrom), levels = misordered.levels), runLength(chrom)),
    start = start,
    end   = start + 200,
    pred  = predictions)
rm(predictions)  # No longer needed
rm(start)  # No longer needed
# Add true binding values to data frame if we know them
if (cell.valid %in% names(chip)) {
  chip.valid <- chip[[as.character(cell.valid)]]
  out$bound <- as.factor(chip.valid[valid.keep])
}
message('Writing predictions to: ', predictions.path)
data.table::fwrite(
  out %>% arrange(chrom, start),
  predictions.path,
  col.names = FALSE,
  sep = '\t')
rm(out)


#
# Session information
#
print(date())
print(Sys.info()[["nodename"]])
devtools::session_info()
largest.objects()
