#!/usr/bin/env Rscript
#
# Make features from motif scan
#
"Usage: feature-scan.R [--tf=TF] [--cell=CELL]... [--wellington] SCANDIR SCANTAG" -> doc


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
.args <- "--cell=A549 --wellington /home/john/Dev/DREAM-ENCODE/Data/Motifs/Known/ TestScan"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
scan.dir <- opts[['SCANDIR']]
scan.tag <- opts[['SCANTAG']]
cells <- opts$cell
tf <- opts$tf
wellington <- opts$wellington


#
# Set up
#
stopifnot(file.exists(scan.dir))
stopifnot(file.exists(file.path(features.dir, scan.tag)))
if (file.exists(features.file.name)) {
  message('Features already exist, stopping: ', features.file.name)
  quit(save = 'no')
}
if (wellington && ! length(cells)) {
  stop('Wellington: no cells specified.')
}


#
# Load motifs
#
.scan <- load.motif.dir(scan.dir)


#' Generate feature for each motif
#'
calc.feature <- function(hits.gr) {
  overlaps <- as.data.frame(findOverlaps(ranges.test(), hits.gr, ignore.strand = TRUE))
  overlaps$value <- mcols(hits.gr[overlaps$subjectHits])$logBF
  with(
    aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
    Rle.from.sparse(length(ranges.test()), queryHits, value))
}


#' Subset ranges by footprints
#'
subsetByFootprints <- function(gr, cell) {
  message('Subsetting by footprints: ', cell)
  footprints <- load.wellington(cell)
  subsetByOverlaps(gr, footprints)
}


#' Save features
#'
save.features <- function(features, features.file.name) {
  message('Saving features: ', features.file.name)
  saveRDS(do.call(S4Vectors::DataFrame, features), features.file.name)
}


#
# Are we subsetting by Wellington footprints?
#
if (wellington) {
  for (cell in cells) {
    message('Calculating features in: ', cell)
    #
    # If a TF is specified name the feature file with it
    #
    if (is.null(tf)) {
      features.file.name <- feature.cell.file.name(scan.tag, cell)
    } else {
      features.file.name <- feature.tf.cell.file.name(scan.tag, tf, cell)
    }
    #
    # Calculate and save features
    features <-
      lapply(
        lapply(.scan, functional::Curry(subsetByFootprints, cell = cell)),
        calc.feature)
    save.features(features, features.file.name)
  }
}

#
# If a TF is specified name the feature file with it
#
if (is.null(tf)) {
  features.file.name <- feature.file.name(scan.tag)
} else {
  features.file.name <- feature.tf.file.name(scan.tag, tf)
}
#
# Calculate and save features
features <- lapply(.scan, calc.feature)
save.features(features, features.file.name)
