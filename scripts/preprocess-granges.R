#!/usr/bin/env Rscript
#
# Load RDS GRanges objects and preprocess with GNCList then save again
#
"Usage: preprocess-granges.R GR..." -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load libraries
#
library(GenomicRanges)


#
# Parse options
#
# .args <- "../Data/Motifs/Known/scan-gr.rds"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
rds.files <- opts[['GR']]


#
# Convert the ranges
#
for (rds.file in rds.files) {
  message('Loading: ', rds.file)
  gr <- readRDS(rds.file)
  if (is(gr, "list")) {
    message('A list: ', rds.file)
    gr.pp <- lapply(gr, GNCList)
  } else if (is(gr, "GNCList")) {
    message('Already a GNCList, skipping: ', rds.file)
    next
  } else {
    message('Converting: ', rds.file)
    gr.pp <- GNCList(gr)
  }
  message('Saving: ', rds.file)
  saveRDS(gr.pp, rds.file, '.tmp')
}
