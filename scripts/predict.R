#!/usr/bin/env Rscript
#
# Simple prediction script. Analyses data from other cell types for a TF to
# predict binding in the VALIDATIONCELL type.
#
"Usage:
predict.R [--options] [--motifs=TAG ...] TF VALIDATIONCELL

Options:
  -o --output DIR        Specify output directory [default: ./predictions]
  -l --ladder            Use ladder cell types [default: FALSE]
  -s --submit            Use submission cell types [default: FALSE]
  -m --motifs=TAG ...    Use motif features from TAG motifs [default='Known']
  -t --test              Make predictions on all test regions rather
                         than just the validation chromosomes [default: FALSE]" -> doc


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
# .args <- "--motif=Known --motif=Known2 GATA3 A549"
# Use dummy arguments if they exist otherwise use command line arguments
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
tf <- factor(opts$TF, tf.levels)
if (is.na(tf)) stop('Unknown TF specified.')
cell.valid <- factor(opts$VALIDATIONCELL, levels = cell.levels)
if (is.na(cell.valid)) stop('Unknown validation cell specified.')
motifs.tags <- opts$motifs


#
# Which cells do we have data for?
#
tf.cells <- tfs %>% filter(TF == tf)
if (! opts$ladder) tf.cells <- tf.cells %>% filter(split != 'ladder')
if (! opts$submit) tf.cells <- tf.cells %>% filter(split != 'submit')
cell.all <- tf.cells$cell
if (! cell.valid %in% cell.all) stop('We have no data for the validation cell type and this TF.')
cell.train <- filter(tf.cells, cell != cell.valid)$cell
if (! length(cell.train)) stop('We have no training data for this cell')


#
# Fix training and validation chromosomes
#
chrs.train <- factor(stringr::str_c('chr', c(3:6, 9:19, 22)), levels = chr.levels)
if (opts$test) {
  chrs.valid <- factor(chr.levels, levels = chr.levels)
} else {
  chrs.valid <- factor(c('chr2', 'chr7', 'chr20'), levels = chr.levels)
}
message('Chromosomes for training:   ', paste(chrs.train, sep = ', '))
message('Chromosomes for validation: ', paste(chrs.valid, sep = ', '))


#
# Load the binding data and expand to full test range
#
chip <- load.chip.features(tf)
chip.cols <- colnames(chip)
chip <- do.call(S4Vectors::DataFrame, lapply(chip, train.factor.Rle.to.test))
colnames(chip) <- chip.cols


#
# Load DNase features
#
dnase <- load.dnase.features(cell.all)


#
# Load motif features
#
load.motif.features <- function(motifs.tag) {
  motif.feature.dir <- file.path(saturn.data(), 'Features', 'Motifs', motifs.tag)
  message('Loading motifs from: ', motif.feature.dir)
  motifs.meta <- readRDS(file.path(motif.feature.dir, 'motif-names.rds'))
  message('Loaded ', nrow(motifs.meta), ' motifs')
  #
  # Load Rle motif scores for each motif
  motif.features <- lapply(
    1:nrow(motifs.meta),
    function(i) readRDS(file.path(motif.feature.dir, basename(motifs.meta$feature.file[i]))))
  names(motif.features) <- motifs.meta$motif
  filter.motif.features <- function(keep) {
    motif.keep <- lapply(motif.features, function(f) f[keep])
    names(motif.keep) <- names(motif.features)
    do.call(S4Vectors::DataFrame, motif.keep)
  }
}
motif.features <- lapply(motifs.tags, load.motif.features)


#
# Determine which regions DNase levels are non-zero
# Any regions with zero DNase-levels we will assume are
# non-binding.
#
non.zero <- lapply(dnase, function(x) x != 0)
names(non.zero) <- cell.all
num.non.zero <- Reduce('+', lapply(non.zero, sum), 0)
message('# non-zero regions across all cell types: ', num.non.zero)


#
# Create data from the features and response, ignoring those regions without any DNase hits
# and those that are not in our training or validation sets.
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
#' Wrangle the data for the given cell
#'
cell.data <- function(cell) {
  message('Wrangling data for: ', cell)
  # Work out which regions to keep (ignore ambiguously bound regions)
  keep <- regions.for.cell(cell)
  # Get all the motif features
  motif.feats <- do.call(cbind, lapply(motif.features, function(feat) feat(keep)))
  # Return S4Vectors::DataFrame
  cbind(
    S4Vectors::DataFrame(
      # cell  = Rle(cell),
      # chrom = Rle(regions.test$chrom[keep]),
      # start = regions.test$start[as.vector(keep)],
      dnase = Rle(dnase[[as.character(cell)]][keep]),
      bound = chip[[as.character(cell)]][keep]),
    motif.feats)
}
#' Remove ambiguously bound regions from DataFrames
#'
remove.ambiguously.bound <- function(df) df['A' != df$bound,]
#' Convert bound factor to integer
#'
bound.factor.to.numeric <- function(bound) {
  is.unbound <- 'U' == bound
  Rle(ifelse(runValue(is.unbound), 0, 1), runLength(is.unbound))
}
# Build the training data
df.train <- remove.ambiguously.bound(Reduce(rbind, lapply(cell.train, cell.data)))
# Build the validation data
df.valid <- Reduce(rbind, lapply(cell.valid, cell.data))
message('Training data size   : ', object.size(df.train))
message('Validation data size : ', object.size(df.valid))
message('# training regions   : ', nrow(df.train))
message('# validation regions : ', nrow(df.valid))


#
# Fit
#
message('Converting DataFrames to sparse matrices')
df.train.nz <- df.train[0 != df.train$dnase,]
mat.train <- as(subset(df.train.nz, select = -c(bound)), "Matrix")
mat.valid <- as(subset(df.valid   , select = -c(bound)), "Matrix")
message('Training matrix size   : ', object.size(mat.train))
message('Validation matrix size : ', object.size(mat.valid))
# Fit a logistic regression
message('Fitting model')
response <- as.factor(drop.ambiguous.level(df.train.nz$bound))
system.time(cvfit <- glmnet::cv.glmnet(mat.train, response, family = 'binomial'))
summary(cvfit)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")
# Save fit
motifs.id <- do.call(stringr::str_c, c(motifs.tags, list(sep = "_")))
fit.id <- stringr::str_c(as.character(tf), '.', as.character(cell.valid), '.', motifs.id)
fit.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('fit.', fit.id, '.rds'))
message('Saving fit: ', fit.path)
saveRDS(cvfit, fit.path)


#
# Make predictions on validation data
#
lambda.predict <- "lambda.1se"
system.time(predictions <- predict(cvfit, mat.valid, s = lambda.predict)[,1])
# ggplot2::qplot(predictions) + ggplot2::scale_y_log10()

#
# Write predictions
#
# Write predictions in format needed for submission
valid.keep <- regions.for.cell(cell.valid)
start <- regions.test$start[as.vector(valid.keep)]
out <- data.frame(
  chrom = Rle(regions.test$chrom[valid.keep]),
  start = start,
  end   = start + 200,
  pred  = logit.inv(predictions))
# Add true binding values to data frame if we know them
if (cell.valid %in% names(chip)) {
  chip.valid <- chip[[as.character(cell.valid)]]
  out$bound <- chip.valid[valid.keep]
}
predictions.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('predictions.', fit.id, '.tsv'))
message('Writing predictions to: ', predictions.path)
readr::write_tsv(out, predictions.path, col_names = FALSE)


#
# Session information
#
print(date())
print(Sys.info()[["nodename"]])
devtools::session_info()
largest.objects()
