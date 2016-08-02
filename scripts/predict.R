#!/usr/bin/env Rscript
#
# Simple prediction script
#
"Usage: predict.R TF VALIDATIONCELL
-o --output DIR  specify output directory [default: ./predictions]
-l --ladder      Use ladder cell types [default: FALSE]
-s --submit      Use submission cell types [default: FALSE]" -> doc


#
# Load libraries
devtools::load_all()
library(Saturn)


#
# Parse options
#
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
tf <- factor(opts$TF, tf.levels)
if (is.na(tf)) stop('Unknown TF specified.')
cell.valid <- factor(opts$VALIDATIONCELL, levels = cell.levels)
if (is.na(cell.valid)) stop('Unknown validation cell specified.')




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

# levels(tfs$cell)
# levels(tf.cells$cell)
# tfs <- tfs %>% mutate(cell = as.character(cell)) %>% mutate(cell = factor(cell, levels = cell.levels))
# devtools::use_data(tfs, overwrite = TRUE)


#
# Fix training and validation chromosomes
#
chrs.train <- factor(stringr::str_c('chr', c(3:6, 9:19, 22)), levels = chr.levels)
chrs.valid <- factor(c('chr2', 'chr7', 'chr20'), levels = chr.levels)


#
# Load the binding data
#
chip <- load.chip.features(tf)


#
# Load DNase features and determine which are non-zero
#
dnase <-
  lapply(
    cell.all,
    function(cell) readRDS(file.path(
      saturn.data(), 'Features', 'DNase',
      stringr::str_c('dnase-summary-', cell, '.rds'))))
names(dnase) <- cell.all
non.zero <- lapply(dnase, function(x) x != 0)
names(non.zero) <- cell.all
num.non.zero <- Reduce('+', lapply(non.zero, sum), 0)
message('# non-zero regions: ', num.non.zero)


#
# Load motif features for regions with non-zero DNase
#
motif.feature.dir <- file.path(saturn.data(), 'Features', 'Motifs', 'Known')
motifs.meta <- readRDS(file.path(motif.feature.dir, 'motif-names.rds'))
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


#
# Create data from the features and response, ignoring those regions without any DNase hits
# and those that are not in our training or validation sets.
#
train.idxs <- regions.test$chrom %in% chrs.train
valid.idxs <- regions.test$chrom %in% chrs.valid
regions.for.cell <- function(cell) {
    if (cell %in% cell.train) {
        train.idxs
    } else {
        stopifnot(cell %in% cell.valid)
        valid.idxs
    }
}
cell.data <- function(cell) {
    message('Wrangling data for: ', cell)
    # Work out which regions to keep
    keep <- non.zero[[cell]] & regions.for.cell(cell)
    # Return S4Vectors::DataFrame
    cbind(
      S4Vectors::DataFrame(
        cell  = Rle(factor(cell, levels = cell.levels)),
        chrom = Rle(regions.test$chrom[keep]),
        start = regions.test$start[as.vector(keep)],
        dnase = Rle(dnase[[cell]][keep]),
        bound = chip[[cell]][keep[regions.train$test.idx]]),
      filter.motif.features(keep))
}
df <- Reduce(rbind, lapply(names(dnase), cell.data))
# sapply(df, class)
message('Data size: ', object.size(df))


#
# Training/validation split
# Match features to regions and split training from validation (ignoring ambiguous training regions)
#
df.train <- df[df$cell %in% cell.train & df$chrom %in% chrs.train & 'A' != df$bound,]
df.train$bound <- Rle(ifelse('U' == df.train$bound, 0, 1))
df.valid <- df[df$cell %in% cell.valid & df$chrom %in% chrs.valid & 'A' != df$bound,]
df.valid$bound <- Rle(ifelse('U' == df.valid$bound, 0, 1))
message('# training regions  : ', nrow(df.train))
message('# validation regions: ', nrow(df.valid))


#
# Fit
#
# Convert the Rle DataFrames to sparse matrices
mat.train <- as(subset(df.train, select = -c(cell, chrom, start, bound)), "Matrix")
mat.valid <- as(subset(df.valid, select = -c(cell, chrom, start, bound)), "Matrix")
# Fit a logistic regression
glmnet.y <- factor(ifelse(df.train$bound, 'B', 'U'), levels=c('U', 'B'))
system.time(cvfit <- glmnet::cv.glmnet(mat.train, glmnet.y, family = 'binomial'))
summary(cvfit)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")


#
# Validate
#
# Make predictions on validation data
predictions <- predict(cvfit, mat.valid, s = "lambda.min")
predictions <- predict(cvfit, mat.valid, s = "lambda.1se")
dim(predictions)
# Assess quality of predictions
labels <- as.vector(df.valid$bound)
scores <- predictions[,1]
length(labels)
length(scores)
prg_curve = prg::create_prg_curve(labels, scores)
auprg = prg::calc_auprg(prg_curve)
# convex_hull = prg::prg_convex_hull(prg_curve)
fig = prg::plot_prg(prg_curve)
# print(prg_curve)
message('AUPRG: ', tf, ', ', cell.valid, ', ', auprg)
# print(fig)


#
# Write predictions
#
# Write predictions in format needed for submission
out <- data.frame(
  chrom = df.valid$chrom,
  start = df.valid$start,
  end   = as.integer(df.valid$start + 200),
  pred  = logit.inv(predictions[,1]))
predictions.file <- stringr::str_c('predictions-', as.character(tf),
                                   '-', as.character(cell.valid), '.tsv')
readr::write_tsv(out, predictions.file, col_names = FALSE)


#
# Session information
#
print(date())
print(Sys.info()[["nodename"]])
devtools::session_info()
largest.objects()
