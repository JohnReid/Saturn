#!/usr/bin/env Rscript
#
# Examine DNase levels against binding
#
"Usage: dnase-binding.R TF" -> doc


#
# Load libraries
devtools::load_all()
library(Saturn)


#
# Parse options
#
# .args <- "E2F1"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
tf <- factor(opts$TF, tf.levels)
if (is.na(tf)) stop('Unknown TF specified.')


#
# Which cells do we have data for?
#
tf.cells <- tfs %>% filter(TF == tf, 'train' == split)
cell.all <- tf.cells$cell


#
# Load the binding data
#
chip <- load.chip.features(tf)


#
# Load DNase features and determine which are non-zero
# Any regions with zero DNase-levels we will assume are
# non-binding
#
dnase <- load.dnase.features(cell.all)
non.zero <- lapply(dnase, function(x) x != 0)
names(non.zero) <- cell.all
num.non.zero <- Reduce('+', lapply(non.zero, sum), 0)
message('# non-zero regions across all cell types: ', num.non.zero)


#
# For each cell type
#
for (cell in cell.all) {
  message(cell)
  cell.df <- S4Vectors::DataFrame(bound = chip[[cell]], dnase = dnase[[cell]][training.region.test.idxs()])
  zero.dnase <- 0 == cell.df$dnase
  is.bound <- 'B' == cell.df$bound
  is.ambig <- 'A' == cell.df$bound
  is.unbnd <- 'U' == cell.df$bound
  by.bound <- function(FUN) c(
    FUN(cell.df$dnase[is.bound]),
    FUN(cell.df$dnase[is.ambig]),
    FUN(cell.df$dnase[is.unbnd]))
  dnase.by.bound <- data.frame(bound = factor(c('B', 'A', 'U'), levels = binding.levels), mean = by.bound(mean), sd = by.bound(sd))
  prop.binding.w.no.dnase <- sum(is.bound & zero.dnase) / sum(is.bound)
  print(dnase.by.bound)
  message('Proportion binding with no DNase: ', prop.binding.w.no.dnase)
}

#
# Create data from the features and response, ignoring those regions without any DNase hits
# and those that are not in our training or validation sets.
#
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
# Make predictions on validation data
#
# predictions <- predict(cvfit, mat.valid, s = "lambda.min")
predictions <- predict(cvfit, mat.valid, s = "lambda.1se")[,1]


#
# Write predictions
#
# Write predictions in format needed for submission
out <- data.frame(
  chrom = df.valid$chrom,
  start = df.valid$start,
  end   = as.integer(df.valid$start + 200),
  pred  = logit.inv(predictions))
dim(out)
# Add truth binding values to data frame if we know them
if (cell.valid %in% names(chip)) {
  chip.valid <- chip[[as.character(cell.valid)]]
  nz.valid <- non.zero[[as.character(cell.valid)]][training.region.test.idxs()]
  idxs.valid <- valid.idxs[training.region.test.idxs()]
  out$bound <- chip.valid['A' != chip.valid & nz.valid & idxs.valid]
}
predictions.file <- stringr::str_c('predictions-', as.character(tf),
                                   '-', as.character(cell.valid), '.tsv')
predictions.path <- file.path(saturn.data(), 'Predictions', predictions.file)
message('Writing predictions to: ', predictions.path)
readr::write_tsv(out, predictions.path, col_names = FALSE)


#
# Session information
#
print(date())
print(Sys.info()[["nodename"]])
devtools::session_info()
largest.objects()
