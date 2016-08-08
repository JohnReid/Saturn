#' The directory for ChIP features
#'
chip.features.dir <- function() file.path(saturn.data(), 'Features', 'ChIP')


#' The file for the TF-cell ChIP features
#'
chip.feature.file <- function(tf, cell) file.path(chip.features.dir(), stringr::str_c('labels-', tf, '-', cell, '.rds'))


#' Load the ChIP features for the TF-cell combo
#'
load.chip.feature <- function(tf, cell) {
  message('Loading ChIP features for ', tf, ' in cell ', cell)
  readRDS(chip.feature.file(tf, cell))
}


#' Load all the ChIP features for the TF and return as a DataFrame
#'
load.chip.features <- memoise::memoise(function(tf) {
  message('Loading ChIP features for ', tf)
  #
  # Get the TF-cell combinations we should have features for
  tf.cells <- tfs %>% filter(TF == tf, 'train' == split)
  .args <- lapply(tf.cells$cell, function(cell) load.chip.feature(tf, cell))
  res <- do.call(S4Vectors::DataFrame, .args)
  names(res) <- as.character(tf.cells$cell)
  res
})


#' The file with the TF's binding labels
#'
tf.chip.labels.file <- function(tf)
  file.path(saturn.data(), 'ChIPseq', 'labels', stringr::str_c(tf, '.train.labels.tsv.gz'))

#' Read the ChIP binding labels for the TF into a S4Vectors::DataFrame
#'
read.chip.labels <- memoise::memoise(function(tf) {
  path <- tf.chip.labels.file(tf)
  message('Loading: ', path)
  # Read ChIP labels into data.table
  labels.dt <- data.table::fread(stringr::str_c('zcat ', path), sep = '\t', drop = 'stop', verbose = FALSE) %>% rename(chrom = chr)
  # Make chrom into a factor with correct levels
  labels.dt$chrom <- factor(labels.dt$chrom, levels = chr.levels)
  # Sort by chrom then start
  data.table::setkey(labels.dt, chrom, start)
  # Convert into S4Vectors::Dataframe using Rle
  # First convert binding factors
  label.cells <- colnames(labels.dt)[3:ncol(labels.dt)]
  binding.rles <- lapply(label.cells, function(l) S4Vectors::Rle(factor(labels.dt[[l]], levels=binding.levels)))
  names(binding.rles) <- label.cells
  # Now convert chromosomes
  df.args <- c(list(chrom=S4Vectors::Rle(labels.dt$chrom), start=labels.dt$start, check.names = FALSE), binding.rles)
  # Check the labels are ordered exactly the same as the test regions
  stopifnot(all(df.args$chrom == regions.train$chrom))
  stopifnot(all(df.args$start == regions.train$start))
  # Build DataFrame
  do.call(S4Vectors::DataFrame, df.args)
})


#' Melt ChIP data into long format
#'
chip.melt <- function(chip) {
  # Melt the binding vector
  binding.m <- do.call(c, lapply(colnames(chip)[3:ncol(chip)], function(col) chip[[col]]))
  # Melt the cell names
  cells.m <- do.call(
      c,
      lapply(colnames(chip)[3:ncol(chip)],
              function(col) S4Vectors::Rle(factor(col, levels = cell.levels), nrow(chip))))
  # Create DataFrame
  DataFrame(
    cell  = cells.m,
    chrom = rep(chip$chrom, ncol(chip) - 2),
    start = rep(chip$start, ncol(chip) - 2),
    bound = binding.m)
}

#' Convert melted ChIP data into data.table
#'
chip.data.table <- function(chip.m) data.table::setkey(
  data.table(
    cell  = as.factor(chip.m$cell),
    chrom = as.factor(chip.m$chrom),
    start = as.vector(chip.m$start),
    bound = as.factor(chip.m$bound)),
  cell, chrom, start)


#' Read ChIP-seq labels into data frame
read.chip.labels.old <- function(tf) {
  path <- file.path(saturn.data(), 'ChIPseq', 'labels', stringr::str_c(tf, '.train.labels.tsv.gz'))
  message('Loading: ', path)
  res <- readr::read_tsv(
    path, progress = FALSE,
    col_types = cols(
      chr=col_factor(seqnames(hg19)),
      start=col_integer(),
      stop=col_integer(),
      .default=col_factor(c('U', 'A', 'B'))))
  names(res) <- rify(names(res))
  res
}


#' Load ChIP training labels into GRanges
load.chip.labels <- memoise::memoise(function(tf) labels.granges(read.chip.labels(tf)))


#' Load ChIP peaks
load.chip.peaks <- memoise::memoise(function(cell, tf, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'ChIPseq', 'peaks', type,
              stringr::str_c('ChIPseq.', cell, '.', tf, '.', .type, '.narrowPeak.gz'))))
})


#' Drop the ambiguous level from the factor-Rle
#'
drop.ambiguous.level <- function(bound) Rle(factor(runValue(bound), levels = c('U', 'B')), runLength(bound))


