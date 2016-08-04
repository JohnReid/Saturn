#!/usr/bin/env Rscript
#
# Examine DNase levels against binding
#
"Usage: dnase-binding.R TF" -> doc


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
