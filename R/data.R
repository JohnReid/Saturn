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
saturn.data <- function() getOption('saturn.data',
                                    system.file('Data', package='Saturn'))


#' Load a narrowPeak file.
load.narrowpeak <- function(path) readr::read_tsv(path, col_names = .NARROWPEAK.COLS)


#' Convert a narrowPeak data frame to GRanges.
narrowpeak.granges <- function(narrowpeak) with(narrowpeak,
  GRanges(
    seqnames = Rle(chrom),
    ranges = IRanges(start = chromStart, end = chromEnd),
    signal = signalValue,
    pValue = pValue,
    qValue = qValue,
    peak = peak))


#' Convert binding factor to numeric
binding.as.numeric <- function(binding) ifelse('B' == binding, 1, ifelse('A' == binding, .5, 0))


#' Convert a ChIP-seq labels data frame to GRanges.
labels.granges <- function(labels) with(labels,
  GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(start = start, end = stop),
    mcols = labels %>% dplyr::select(-chr, -start, -stop)))


#' Load ChIP training labels
load.chip.labels <- memoise::memoise(function(tf) {
  labels.granges(readr::read_tsv(
    file.path(saturn.data(),
    'ChIPseq',
    'labels',
    str_c(tf, '.train.labels.tsv.gz'))))
})


#' Load expression data.
saturn.expr <- memoise::memoise(function(cell, biorep) {
  readr::read_tsv(file.path(saturn.data(),
                     'RNAseq',
                     sprintf('gene_expression.%s.biorep%d.tsv', cell, biorep)))
})


#' Combine ChIP and DNAse data
combine.chip.dnase <- function(chip.labels, dnase) {
  # Find the overlaps between the DNAse data and the ChIP labels
  overlaps <- as.data.frame(findOverlaps(dnase, chip.labels))
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



#' Load ChIP peaks
load.chip.peaks <- memoise::memoise(function(cell, tf, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'ChIPseq', 'peaks', type,
              str_c('ChIPseq.', cell, '.', tf, '.', .type, '.narrowPeak.gz'))))
})


#' Load DNase peaks
load.dnase.peaks <- function(cell, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'DNASE', 'peaks', type,
              str_c('DNASE.', cell, '.', type, '.narrowPeak.gz'))))
}

#' Load motif scan results
load.motif.scan <- memoise::memoise(function(results.path, seqs.path) {
  motif.seqs <- readr::read_csv(seqs.path)
  with(
    readr::read_csv(
      results.path, skip = 1, col_types = 'cciicnnn', progress = FALSE,
      col_names = c('motif', 'w.mer', 'seq', 'position', 'strand',
                    'Z', 'score', 'p.value')) %>%
    mutate(chr = motif.seqs$ID[seq+1],
           end = position + str_length(w.mer)),
    GRanges(
      seqnames = Rle(chr),
      ranges = IRanges(start = position, end = end),
      strand = strand,
      motif = motif,
      Z = Z,
      score = score,
      p.value = p.value))
})
