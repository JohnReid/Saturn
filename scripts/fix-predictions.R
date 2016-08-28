#!/usr/bin/env Rscript
#
# Fix predictions that have floating point end values
#
"Usage: fix-predictions [INPUT ...]" -> doc


#
# Load packages
#
library(docopt)
library(data.table)


#
# Parse options
#
# .args <- "/home/jer15/Dev/Saturn/Data/Predictions/tmp"
# .args <- "/home/jer15/Dev/Saturn/Data/Predictions/predictions.TAF1.HepG2.Known_DREME-TAF1.tsv"
# Use dummy arguments if they exist otherwise use command line arguments
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
prediction.files <- opts$INPUT


#
# Fix files by reading them in and writing them out
#
for (prediction.file in prediction.files) {
  message('Reading: ', prediction.file)
  time.taken <- system.time(predictions <- fread(prediction.file, header = FALSE))
  print(time.taken)
  message('Fixing: ', prediction.file)
  time.taken <- system.time(fwrite(
    predictions,
    prediction.file,
    col.names = FALSE,
    sep = '\t'))
  print(time.taken)
}


# test
out.file <- stringr::str_c(prediction.file, '.test')
time.taken <- system.time(fwrite(
  predictions[, 1:3, with = FALSE],
  out.file,
  col.names = FALSE,
  sep = '\t'))

