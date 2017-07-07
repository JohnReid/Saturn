#' The directory for features
#'
features.dir <- function() file.path(saturn.data(), 'Features')

feature.tf.cell.file.name <- function(feature.name, tf, cell) file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-', tf, '-', cell, '.rds'))
feature.tf.file.name      <- function(feature.name, tf)       file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-', tf,            '.rds'))
feature.cell.file.name    <- function(feature.name, cell)     file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-',          cell, '.rds'))
feature.file.name         <- function(feature.name)           file.path(features.dir(), feature.name, stringr::str_c(feature.name,                     '.rds'))

#' Get file for feature
#'
feature.file <- function(feature.name, tf, cell) {
  tf.cell.file <- feature.tf.cell.file.name(feature.name, tf, cell)
  tf.file      <- feature.tf.file.name     (feature.name, tf      )
  cell.file    <- feature.cell.file.name   (feature.name,     cell)
  .file        <- feature.file.name        (feature.name          )
  if (file.exists(tf.cell.file)) return(tf.cell.file)
  if (file.exists(tf.file))      return(tf.file)
  if (file.exists(cell.file))    return(cell.file)
  if (file.exists(.file))        return(.file)
  stop(stringr::str_c('No features for this combination: ', feature.name, ', ', tf, ', ', cell))
}


#' Load a feature
#'
load.feature <- function(feature.name, tf, cell) load.feature.from.file(feature.file(feature.name, tf, cell))


#' Read a feature from a file (cached)
#'
read.feature.from.file <- memoise::memoise(function(file.name) {
  message('Loading feature: ', file.name)
  readRDS(file.name)
})


#' Load a feature from a file and convert to a Matrix
#'
load.feature.from.file <- function(file.name, keep.idx = NULL) {
  feat <- read.feature.from.file(file.name)
  if (! is.null(keep.idx)) {
    feat <- feat[keep.idx, ]
  }
  as(feat, "Matrix")
}


#' Make a feature by calculating the width of its intersect with each test range
#'
feature.from.ranges <- function(feat.gr) {
  feat.reduced <- reduce(feat.gr)
  #
  # Calculate the overlaps between the test ranges and the feature locations
  system.time(overlaps <- findOverlaps(ranges.test.gnclist(), reduce(feat.reduced)))
  #
  # Calculate the widths of the overlaps per test range
  overlapwidths <-
      as.data.frame(overlaps) %>%
      mutate(width = width(pintersect(ranges.test()[queryHits], feat.reduced[subjectHits]))) %>%
      group_by(queryHits) %>%
      summarise(width = sum(width))
  #
  # Convert into a run-length-encoded vector feature
  with(overlapwidths, Rle.from.sparse(length(ranges.test()), queryHits, width))
}
