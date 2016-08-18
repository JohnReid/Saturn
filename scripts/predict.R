#!/usr/bin/env Rscript
#
# Simple prediction script. Analyses data from other cell types for a TF to
# predict binding in the VALIDATIONCELL type.
#
"Usage:
predict.R [-t] [--sample=PROPR] [--use-zero-dnase] [--motifs=TAG]... TF VALIDATIONCELL

Options:
  -m --motifs=TAG ...    Use motif features from TAG motifs [default='Known']
  --use-zero-dnase       Fit regions with zero DNase levels [default='False']
  -s --sample=PROP       Down-sample training regions to PROP [default=1]" -> doc


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
# .args <- "--motif=Known --motif=DREME-GABPA -s .1 GABPA SK-N-SH"
# Use dummy arguments if they exist otherwise use command line arguments
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
tf <- factor(opts$TF, tf.levels)
if (is.na(tf)) stop('Unknown TF specified.')
cell.valid <- factor(opts$VALIDATIONCELL, levels = cell.levels)
if (is.na(cell.valid)) stop('Unknown validation cell specified.')
motifs.tags <- opts$motifs
sample.prop <- as.numeric(opts$sample)
use.zero.dnase <- as.logical(opts[['use-zero-dnase']])


#
# Construct output filenames
#
motifs.id <- do.call(stringr::str_c, c(motifs.tags, list(sep = "_")))
fit.id <- stringr::str_c(as.character(tf), '.', as.character(cell.valid), '.', motifs.id)
fit.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('fit.', fit.id, '.rds'))
predictions.path <- file.path(saturn.data(), 'Predictions', stringr::str_c('predictions.', fit.id, '.tsv'))


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
# Load DNase features
#
dnase <- load.dnase.features(cell.all)
message('DNASE data size: ', object.size(dnase))


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
# Create data from the features and response, ignoring those regions
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
#' Wrangle the data for the given cell, remove ambiguous binding if we have
#' binding data
#'
cell.data <- function(cell, remove.ambiguous = TRUE) {
  message('Wrangling data for: ', cell)
  # Get the binding status
  cell.bound <- chip[[as.character(cell)]]
  # Work out which regions to keep
  keep <- regions.for.cell(cell)
  if (is.null(cell.bound)) {
    # If we don't have binding information, use NAs
    cell.bound <- Rle(factor(NA, binding.levels), length(keep))
  } else if (remove.ambiguous) {
    # Work out which regions to keep (ignore ambiguously bound regions)
    keep <- keep & ('A' != cell.bound)
  }
  # Subset bound regions
  cell.bound <- cell.bound[keep]
  # Get all the motif features
  motif.feats <- do.call(cbind, lapply(motif.features, function(feat) feat(keep)))
  # Get the DNase levels
  cell.dnase <- Rle(dnase[[as.character(cell)]][keep])
  # Return S4Vectors::DataFrame
  cbind(
    S4Vectors::DataFrame(
      # cell  = Rle(cell),
      # chrom = Rle(regions.test$chrom[keep]),
      # start = regions.test$start[as.vector(keep)],
      dnase = cell.dnase,
      bound = cell.bound),
    motif.feats)
}
# Build the training data
message('Creating training data')
df.train <- do.call(rbind, lapply(cell.train, cell.data))
message('Training data size   : ', object.size(df.train))
message('# training regions   : ', nrow(df.train))


#
# Remove regions with no DNase if requested
#
if (! use.zero.dnase) {
  message('Removing regions without DNase levels')
  df.train <- df.train[0 != df.train$dnase,]
  message('# training regions: ', nrow(df.train))
}


#
# Create training matrix
#
message('Creating training matrix and response')
mat.train <- as(subset(df.train, select = -c(bound)), "Matrix")
response <- factor(as.factor(df.train$bound), levels = c('U', 'B'))
rm(df.train)  # No longer needed


#
# Down sample regions
#
if (sample.prop < 1) {
  message('Sampling training regions, proportion: ', sample.prop)
  .sample <- sample.int(nrow(mat.train), size = as.integer(sample.prop * nrow(mat.train)))
  mat.train <- mat.train[.sample,]
  response <- response[.sample]
  rm(.sample)
  message('# training regions: ', nrow(mat.train))
}


#
# Fit a logistic regression
#
message('Training matrix size: ', object.size(mat.train))
message('# features: ', ncol(mat.train))
message('Fitting model')
system.time(cvfit <- glmnet::cv.glmnet(mat.train, response, family = 'binomial'))
rm(mat.train)  # No longer needed
summary(cvfit)
# plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")
# Save fit
message('Saving fit: ', fit.path)
saveRDS(cvfit, fit.path)


#
# Create validation matrix
#
message('Creating validation data')
df.valid <- cell.data(cell.valid, remove.ambiguous = FALSE)
message('Validation data size : ', object.size(df.valid))
message('# validation regions : ', nrow(df.valid))
mat.valid <- as(subset(df.valid   , select = -c(bound)), "Matrix")
rm(df.valid)  # No longer needed
message('Validation matrix size : ', object.size(mat.valid))


#
# Make predictions on validation data
#
lambda.predict <- "lambda.min"
system.time(predictions <- predict(cvfit, mat.valid, s = lambda.predict)[,1])
rm(mat.valid)  # No longer needed


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
    pred  = logit.inv(predictions))
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
