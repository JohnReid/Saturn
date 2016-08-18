#!/usr/bin/env Rscript
#
# Make features from motif scan
#
"Usage: feature-scan.R SCANDIR SCANTAG" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load libraries
#
devtools::load_all()
library(Saturn)


#
# Parse options
#
# .args <- "/home/john/Dev/DREAM-ENCODE/Data/Motifs/Known/ TestScan"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
scan.dir <- opts[['SCANDIR']]
scan.tag <- opts[['SCANTAG']]


#
# Set up
#
features.file.name <- feature.file.name(scan.tag)
stopifnot(file.exists(scan.dir))
stopifnot(file.exists(dirname(features.file.name)))
if (file.exists(features.file.name)) {
  message('Features already exist, stopping: ', features.file.name)
  quit(save = 'no')
}


#
# Load motifs
#
.scan <- load.motif.dir(scan.dir)


#
# Generate feature for each motif
#
calc.feature <- function(hits.gr) {
  overlaps <- as.data.frame(findOverlaps(ranges.test(), hits.gr, ignore.strand = TRUE))
  overlaps$value <- mcols(hits.gr[overlaps$subjectHits])$logBF
  with(
    aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
    Rle.from.sparse(length(ranges.test()), queryHits, value))
}
features <- lapply(.scan, calc.feature)


#
# Save features
#
message('Saving features: ', features.file.name)
saveRDS(features, features.file.name)
