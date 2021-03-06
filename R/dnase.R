#' The directory for DNase features
#'
dnase.features.dir <- function() file.path(saturn.data(), 'Features', 'DNase')


#' The file for the TF-cell DNase features
#'
dnase.feature.file <- function(cell) file.path(
  dnase.features.dir(), stringr::str_c('dnase-summary-', cell, '.rds'))


#' Load the DNase features for the cell
#'
load.dnase.feature <- memoise::memoise(function(cell) {
  message('Loading DNase features for cell ', cell)
  readRDS(dnase.feature.file(cell))
})


#' Load DNase features
#'
load.dnase.features <- function(cells) {
  l <- lapply(cells, load.dnase.feature)
  names(l) <- cells
  res <- do.call(S4Vectors::DataFrame, l)
  names(res) <- cells
  res
}


#' Combine ChIP and DNAse data
#'
combine.chip.dnase <- function(chip.labels, dnase) {
  # Find the overlaps between the DNAse data and the ChIP labels
  overlaps <- as.data.frame(findOverlaps(dnase, chip.labels, ignore.strand = TRUE))
  # Add the p-values to the overlaps
  overlaps$dnase <- dnase[overlaps$queryHits,]$pValue
  # Summarise the overlaps by the maximum p-value for each ChIP label
  label.dnase <- overlaps %>%
    group_by(subjectHits) %>%
    summarise(dnase=max(dnase))
  dnase <- rep(0, length(chip.labels))
  dnase[label.dnase$subjectHits] <- label.dnase$dnase
  chip.labels$dnase <- dnase
  chip.labels
}



#' Load DNase peaks
#'
load.dnase.peaks <- memoise::memoise(function(cell, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'DNASE', 'peaks', type,
              stringr::str_c('DNASE.', cell, '.', type, '.narrowPeak.gz'))))
})


#' File path for DNase signal
#'
dnase.signal.file <- function(cell) file.path(
  saturn.data(), 'DNASE', 'fold_coverage_wiggles',
  stringr::str_c('DNASE.', cell, '.fc.signal.bigwig'))


#' Load DNase signal
#'
load.dnase.signal <- memoise::memoise(function(cell) {
  rtracklayer::import(dnase.signal.file)
})


#' Summarise DNase peaks by ChIP-regions
#'
summarise.dnase <- function(cell, type='conservative', aggregation.fn=max) {
  # Load the peaks
  dnase <- load.dnase.peaks(cell, type)
  # Aggregate
  agg.by.region(dnase, ranges.test(), 'pValue', aggregation.fn)
}


#' DNase footprint file
#'
wellington.footprint.filename <- function(cell, type = 'relaxed') file.path(
    saturn.data(), 'DNASE', 'bams', 'footprints', cell, type,
    stringr::str_c('DNASE.', cell, '.merged.bam.DNASE.', cell,
                   '.relaxed.bed.WellingtonFootprints.FDR.0.01.bed'))


#' Load Wellington DGF predictions
#'
load.wellington <- memoise::memoise(function(cell, type = 'relaxed') {
  dgf.file <- wellington.footprint.filename(cell = cell, type = type)
  message('Loading Wellington DGF predictions: ', dgf.file)
  rtracklayer::import(dgf.file)
})
