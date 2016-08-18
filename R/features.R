#' The directory for features
#'
features.dir <- function() file.path(saturn.data(), 'Features')

feature.tf.cell.file <- function(feature.name, tf, cell) file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-', tf, '-', cell, '.rds'))
feature.tf.file      <- function(feature.name, tf)       file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-', tf,            '.rds'))
feature.cell.file    <- function(feature.name, cell)     file.path(features.dir(), feature.name, stringr::str_c(feature.name, '-',          cell, '.rds'))

#' Get file for feature
#'
feature.file <- function(feature.name, tf, cell) {
  tf.cell.file <- feature.tf.cell.file(feature.name, tf, cell)
  tf.file      <- feature.tf.file     (feature.name, tf      )
  cell.file    <- feature.cell.file   (feature.name,     cell)
  if (file.exists(tf.cell.file)) return(tf.cell.file)
  if (file.exists(tf.file))      return(tf.file)
  if (file.exists(cell.file))    return(cell.file)
  stop(stringr::str_c('No features for this combination: ', feature.name, ', ', tf, ', ', cell))
}


#' Load a feature
#'
load.feature <- function(feature.name, tf, cell) load.feature.from.file(feature.file(feature.name, tf, cell))


#' Load a feature from a file (cached)
#'
load.feature.from.file <- memoise::memoise(function(feature.file.name) {
  message('Loading feature: ', feature.file.name)
  as(readRDS(feature.file.name), "Matrix")
})
