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
# .args <- "--cell=A549 --wellington /home/john/Dev/DREAM-ENCODE/Data/Motifs/Known/ KnownMotifs"
# .args <- "--cell=H1-hESC --tf=TEAD4 --wellington ..//Data/Motifs/DREME-TEAD4/genome-scan DREME"
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
stopifnot(file.exists(file.path(features.dir(), scan.tag)))
if (wellington && ! length(cells)) {
  stop('No cells specified for Wellington footprints.')
}


#' Cached function to load motifs
#'
get.scan <- memoise::memoise(function() load.motif.dir(scan.dir))


#' Generate feature for each motif
#'
calc.feature <- function(hits.gr) {
  overlaps <- as.data.frame(findOverlaps(ranges.test(), hits.gr, ignore.strand = TRUE))
  if (! nrow(overlaps)) {
    Rle(0, length(ranges.test()))
  } else {
    overlaps$value <- mcols(hits.gr[overlaps$subjectHits])$logBF
    with(
      aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
      Rle.from.sparse(length(ranges.test()), queryHits, value))
  }
}


#' Subset ranges by footprints
#'
subsetByFootprints <- function(gr, cell) {
  footprints <- load.wellington(cell)
  subsetByOverlaps(gr, footprints)
}


#' Save features as DataFrame
#'
save.features <- function(features, features.file.name) {
  message('Saving features: ', features.file.name)
  dir.create(dirname(features.file.name), showWarnings = FALSE)
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
    well.tag <- stringr::str_c(scan.tag, 'Well')
    if (is.null(tf)) {
      features.file.name <- feature.cell.file.name(well.tag, cell)
    } else {
      features.file.name <- feature.tf.cell.file.name(well.tag, tf, cell)
    }
    if (file.exists(features.file.name)) {
      message('Features already exist, not recreating: ', features.file.name)
    } else {
      #
      # Calculate and save features (subsetted by footprints)
      message('Subsetting by footprints: ', cell)
      features <-
        lapply(
          lapply(get.scan(), functional::Curry(subsetByFootprints, cell = cell)),
          calc.feature)
      names(features) <- stringr::str_c(names(features), '.Well')
      save.features(features, features.file.name)
    }
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
if (file.exists(features.file.name)) {
  message('Features already exist, not recreating: ', features.file.name)
} else {
  #
  # Calculate and save features
  message('Calculating features without footprints: ')
  features <- lapply(get.scan(), calc.feature)
  save.features(features, features.file.name)
}
