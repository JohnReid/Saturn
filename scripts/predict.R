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
#SBATCH --time=24:00:00                 # How much wallclock time will be required?
#SBATCH -o predict/predict-%j.out       # stdout
#SBATCH -e predict/predict-%j.out       # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
                                        # interrupted by node failure or system downtime):
# Sys.setenv(OMP_NUM_THREADS=16)          # Set OpenMP number of threads
"Usage:
predict.R [options] [--features=NAME]... TF VALIDATIONCELL

Options:
  --method=METHOD        Use METHOD [default: glmnet]
  --tag=TAG              Add TAG to results name. [default: '']
  -f --features=NAME     Use NAME features
  --max-boosting=MAXIMUM MAXIMUM number of boost rounds for xgboost method [default: 300]
  --expr                 Use expression summary [default: FALSE]
  -r --remove-zero-dnase Don't use regions with zero DNase levels [default: FALSE]
  -d --down-sample       Down-sample training regions (stratified by DNase) [default: False]
  -s --sample=PROP       Subsample training regions to with proportion PROP [default: 1]" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
devtools::load_all()
library(Saturn)


#
# Parse options
#
# .args <- "--tag=test --motif=Known ARID3A K562"
# .args <- "--tag=test --motif=Known GATA3 A549"
# .args <- "--tag=test -f DNase -f KnownMotifs -s .01 --expr GABPA SK-N-SH"
# .args <- "--tag=test --method=xgboost -f DNase -f DREMEWell --expr -s .01 ATF2 GM12878"
# .args <- "--tag=test --method=xgboost --sample=.1 --max-boosting=30 -f DNase -f DREMEWell ATF2 GM12878"
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
tag <- opts$tag
feat.names <- opts$features
sample.prop <- as.numeric(opts$sample)
max.boost.rounds <- as.integer(opts[['max-boosting']])
remove.zero.dnase <- as.logical(opts[['remove-zero-dnase']])
down.sample <- as.logical(opts[['down-sample']])
use.expr <- as.logical(opts[['expr']])
message('Features: ', toString(feat.names))
if (! 'DNase' %in% feat.names) {
  message('WARNING: No DNase feature included in arguments!!!')
}
message('Sample proportion: ', toString(sample.prop))
message('Down-sample: ', toString(down.sample))
message('Remove zero DNase: ', toString(remove.zero.dnase))
message('Use expression: ', toString(use.expr))


#
# Construct output filenames
#
feat.tags <- do.call(stringr::str_c, c(feat.names, list(sep = "_")))
if (use.expr) feat.tags <- stringr::str_c(feat.tags, '_expr')
fit.id <- stringr::str_c(method, '.', tag, '.', as.character(tf), '.', as.character(cell.valid), '.', feat.tags)
message('Fit ID: ', fit.id)
predictions.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('predictions.', fit.id, '.tsv'))


#
# Load prediction package
#
if ('xgboost' == method) {
  # devtools::load_all(stringr::str_c(Sys.getenv('HOME'), '/src/xgboost/R-package'))
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
  .remove.zero.dnase = remove.zero.dnase,
  .sample.prop = sample.prop,
  .down.sample = down.sample,
  .remove.ambiguous = TRUE,
  .use.expr = use.expr)
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
  chrom <- regions.test$chrom[keep]
  message(cell, ': Proportion bound: ', prop.bound(response))
  #
  # Load and cbind the features
  mat <-
    do.call(
      cbind,
      lapply(feat.names, function(feat.name) load.feature(feat.name, tf, cell)[as.vector(keep), , drop = FALSE]))
  #
  # Down sample the unbound examples to match the bound examples (stratified by non-zero DNase if present)
  if (.down.sample) {
    if ('DNase' %in% colnames(mat)) {
      message('Down sampling (stratified by DNase)')
      dnase.non.zeros <- 0 != mat[,'DNase']
      bound <- 'B' == response
      non.zero.bound   <- which(  dnase.non.zeros &   bound)
      zero.bound       <- which(! dnase.non.zeros &   bound)
      non.zero.unbound <- which(  dnase.non.zeros & ! bound)
      zero.unbound     <- which(! dnase.non.zeros & ! bound)
      # Down sample non zero DNase regions
      if (length(non.zero.bound) < length(non.zero.unbound)) {
        non.zero.unbound <- sample(non.zero.unbound, length(non.zero.bound))
      } else {
        non.zero.bound <- sample(non.zero.bound, length(non.zero.unbound))
      }
      # Down sample zero DNase regions
      if (length(zero.bound) < length(zero.unbound)) {
        zero.unbound <- sample(zero.unbound, length(zero.bound))
      } else {
        zero.bound <- sample(zero.bound, length(zero.unbound))
      }
      down.sampled <- sort(c(non.zero.bound, non.zero.unbound, zero.bound, zero.unbound))
    } else {
      message('Down sampling')
      bound <- 'B' == response
      bound.which   <- which(  bound)
      unbound.which <- which(! bound)
      if (sum(bound) < sum(! unbound)) {
        unbound.which <- sample(unbound.which, length(bound.which))
      } else {
        bound.which <- sample(bound.which, length(unbound.which))
      }
      down.sampled <- sort(c(bound.which, unbound.which))
    }
    mat <- mat[down.sampled,]
    response <- response[down.sampled]
    chrom <- chrom[down.sampled]
    message(cell, ': Proportion bound: ', prop.bound(response))
  }
  #
  # Remove those with zero DNase levels if we don't want them
  if (.remove.zero.dnase && 'DNase' %in% colnames(mat)) {
    dnase.non.zeros <- 0 != mat[,'DNase']
    message(cell, ': Removing regions with no DNase. ',
            'Reducing regions from ', length(dnase.non.zeros),
            ' to ', sum(dnase.non.zeros))
    mat <- mat[dnase.non.zeros,]
    response <- response[dnase.non.zeros]
    chrom <- chrom[dnase.non.zeros]
    message(cell, ': Proportion bound: ', prop.bound(response))
  }
  #
  # Subsample the regions if requested
  if (.sample.prop < 1) {
    .sample <- sample.int(nrow(mat), size = as.integer(.sample.prop * nrow(mat)))
    mat <- mat[.sample,]
    response <- response[.sample]
    chrom <- chrom[.sample]
  }
  #
  # Add expression data if requested
  if (.use.expr) {
    cell.expr <- expr.features.for.cell(cell)
    mat <- cbind(
      mat,
      t(matrix(rep(cell.expr, nrow(mat)), ncol = nrow(mat))))
  }
  message(cell, ': Has ', nrow(mat), ' regions')
  list(
    cell = Rle(cell, length(chrom)),
    chrom = chrom,
    features = mat,
    response = factor(as.factor(response), levels = c('U', 'B')))
}
message('Creating training data')
train.data  <- lapply(cell.train, load.cell.data)
train.feat  <- do.call(rbind, lapply(train.data, function(d) d$features))
train.resp  <- do.call(c,     lapply(train.data, function(d) d$response))
train.cell  <- do.call(c,     lapply(train.data, function(d) d$cell))
train.chrom <- do.call(c,     lapply(train.data, function(d) d$chrom))
train.cell.nrows <- lapply(train.data, function(x) nrow(x$features))
rm(train.data)  # No longer needed
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
  early.stop.round = 50,
  folds_test = NULL,
  folds_train = NULL
) {
  #
  # Perform cross-validation to choose number of rounds
  param <- list(silent = 1, objective = 'binary:logistic')
  message('Cross-validating')
  cv.time <- system.time(
    cv.result <- xgb.cv(
      params = param,
      data = data,
      nround = nround,
      nfold = nfold,
      folds_test = folds_test,
      folds_train = folds_train,
      # maximize = TRUE,
      # metrics = list('map'),
      early_stopping_rounds = early.stop.round))
  print(cv.time)
  nround.best <- cv.result$best_ntreelimit
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
  # Decide on CV folds
  #
  if (length(train.cell.nrows) > 1) {
    #
    # we have enough cells to use CV across cell types
    chr.sets <- list(
      chr.levels[seq(1, length(chr.levels), by = 2)],
      chr.levels[seq(2, length(chr.levels), by = 2)])
    .grid <- expand.grid(cell = unique(runValue(train.cell)), chr.set = c(1, 2))
    folds_train <- lapply(
      1:nrow(.grid),
      function(i) {
        test.cell <- .grid$cell[i]
        test.chroms <- chr.sets[[.grid$chr.set[i]]]
        which(train.cell != test.cell & ! train.chrom %in% test.chroms)
      })
    folds_test <- lapply(
      1:nrow(.grid),
      function(i) {
        test.cell <- .grid$cell[i]
        test.chroms <- chr.sets[[.grid$chr.set[i]]]
        which(train.cell == test.cell & train.chrom %in% test.chroms)
      })
  } else {
    #
    # just let xgb.cv use randomly chosen folds
    folds_test <- NULL
    folds_train <- NULL
  }
  #
  # Fit with xgboost
  #
  message('Fitting model with xgboost')
  dtrain <- xgb.DMatrix(train.feat, label = train.resp - 1)
  fit <- xgboost.fit(data = dtrain, nround = max.boost.rounds, folds_test = folds_test, folds_train = folds_train)
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
valid.feat <- load.cell.data(cell.valid, .remove.zero.dnase = FALSE, .sample.prop = 1, .remove.ambiguous = FALSE, .down.sample = FALSE)$features
message('Validation data size : ', object.size(valid.feat))
message('# validation regions : ', nrow(valid.feat))


#'
#' Create a xgb.DMatrix from a dgCMatrix working around bug in xgboost where zero rows at end
#' of dgCMatrix are ignored.
#'
DMatrix.from.dgC <- function(dgC) {
  if (0 == dgC[nrow(dgC), 1]) {
    dgC[nrow(dgC), 1] <- .Machine$double.eps
  }
  xgb.DMatrix(dgC)
}

#
# Make predictions on validation data
#
message('Making predictions')
if ('xgboost' == method) {
  system.time(predictions <- predict(fit, DMatrix.from.dgC(valid.feat)))
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
