.NARROWPEAK.COLS <- c(
  'chrom',
  'chromStart',
  'chromEnd',
  'name',
  'score',
  'strand',
  'signalValue',
  'pValue',
  'qValue',
  'peak')

chrs.levels <- c(str_c('chr', 1:22), 'chrX')

#' The root directory of the Saturn data.
saturn.data <- function() system.file('Data', package='Saturn')


#' Load expression data.
saturn.expr <- function(cell, biorep) {
  readr::read_tsv(file.path(saturn.data(),
                     'RNAseq',
                     sprintf('gene_expression.%s.biorep%d.tsv', cell, biorep)))
}


#' Load ChIP training labels
load.chip.labels <- function(tf) {
  readr::read_tsv(file.path(saturn.data(),
                            'ChIPseq',
                            'labels',
                            str_c(tf, '.train.labels.tsv.gz')))
}


#' Load a narrowPeak file.
load.narrowpeak <- function(path) readr::read_tsv(path, col_names = .NARROWPEAK.COLS)


#' Load ChIP peaks
load.chip.peaks <- function(cell, tf, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  load.narrowpeak(file.path(saturn.data(),
                            'ChIPseq',
                            'peaks',
                            type,
                            str_c('ChIPseq.', cell, '.', tf, '.', .type, '.narrowPeak.gz')))
}


#' Load DNase peaks
load.dnase.peaks <- function(cell, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  load.narrowpeak(file.path(saturn.data(),
                            'DNASE',
                            'peaks',
                            type,
                            str_c('DNASE.', cell, '.', type, '.narrowPeak.gz')))
}
