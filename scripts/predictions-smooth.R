#!/usr/bin/env Rscript
"Usage: predictions-smooth.R [options] IN OUT

Options:
  -l --length-scale=LENGTHSCALE  Use LENGTHSCALE [default: 50]
  --logodds                      Smooth log-odds instead of probabilities [default: FALSE]
  --width MAXWIDTH               Width of kernel in units of regions [default: 20]" -> doc


#
# Load package
#
devtools::load_all()
library(Saturn)


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
sample_n(preds, 15)
sapply(preds, class)


#
# Smooth predictions
#
predictions.smoothed <- smooth.predictions(preds, length.scale, max.width, log.transform)


#
# Write predictions
#
preds$prediction <- as.vector(predictions.smoothed)
data.table::fwrite(preds, out.path, sep = '\t', col.names = FALSE)
