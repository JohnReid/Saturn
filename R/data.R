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

hg19 <- BSgenome.Hsapiens.UCSC.hg19

#' Make names suitable for R.
#'
#' E.g. substitute '_' for '-'
#'
rify <- function(x) str_replace_all(x, '[-{}]', '_')


#' The root directory of the Saturn data.
saturn.data <- function() getOption('saturn.data',
                                    system.file('Data', package='Saturn'))


#' Load a narrowPeak file.
load.narrowpeak <- function(path) {
  message('Loading: ', path)
  readr::read_tsv(
    path, progress = FALSE,
    col_names = .NARROWPEAK.COLS,
    col_types = cols(
      chrom = col_factor(seqnames(hg19)),
      chromStart = col_integer(),
      chromEnd = col_integer(),
      name = col_character(),
      score = col_integer(),
      strand = col_character(),
      signalValue = col_double(),
      pValue = col_double(),
      qValue = col_double(),
      peak = col_integer()))
}

#' Convert a narrowPeak data frame to GRanges.
narrowpeak.granges <- function(narrowpeak) with(narrowpeak,
  GRanges(
    seqnames = Rle(chrom),
    ranges = IRanges(start = chromStart+1, end = chromEnd),
    seqinfo = seqinfo(hg19),
    signal = signalValue,
    pValue = pValue,
    qValue = qValue,
    peak = peak))


#' Convert binding factor to numeric
binding.as.numeric <- function(binding) ifelse('B' == binding, 1, ifelse('A' == binding, .5, 0))


#' Convert a ChIP-seq labels data frame to GRanges.
labels.granges <- function(labels) {
  mcol.names <- names(labels)[4:ncol(labels)]
  .args = lapply(labels[,mcol.names], function(x) Rle(as.integer(x)))
  .args$seqnames = Rle(labels$chr)
  .args$ranges = IRanges(start = labels$start+1, end = labels$stop)
  .args$seqinfo = seqinfo(hg19)
  do.call(GRanges, args = .args)
}


#' Read ChIP-seq labels into data frame
read.chip.labels <- function(tf) {
  path <- file.path(saturn.data(), 'ChIPseq', 'labels', str_c(tf, '.train.labels.tsv.gz'))
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


#' Load expression data.
saturn.expr <- memoise::memoise(function(cell, biorep) {
  path <- file.path(
    saturn.data(), 'RNAseq',
    sprintf('gene_expression.%s.biorep%d.tsv', cell, biorep))
  message('Loading: ', path)
  readr::read_tsv(
    path, skip = 1, progress = FALSE,
    col_names = c(
      'gene_id',
      'transcript_ids',
      'length',
      'effective_length',
      'expected_count',
      'TPM',
      'FPKM'),
    col_types = cols(
      gene_id = col_character(),
      transcript_ids = col_character(),
      .default = col_double()))
})

#' Combine ChIP and DNAse data
#'
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

#' Summarise DNase peaks by ChIP-regions
summarise.dnase <- memoise::memoise(function(cell, type='conservative') {
  dnase <- load.dnase.peaks(cell, type)
  # Find the overlaps between the DNAse data and the ChIP regions
  overlaps <- as.data.frame(findOverlaps(dnase, regions))
  # Add the p-values to the overlaps
  overlaps$dnase <- mcols(dnase)[overlaps$queryHits, 'pValue']
  # Summarise the overlaps by the maximum p-value for each ChIP label
  label.dnase <- overlaps %>%
    group_by(subjectHits) %>%
    summarise(dnase=max(dnase))
  dnase <- Rle(0, length(regions))
  dnase[label.dnase$subjectHits] <- label.dnase$dnase
  dnase
})

#' Binding for TF/cell type combination
#'
binding.tf.cell <- memoise::memoise(function(tf, cell) mcols(load.chip.labels(tf))[,rify(cell)])

#' Load motif scan results
#'
load.motif.scan <- memoise::memoise(function(
  results.path,
  seqs.path,
  prior.log.odds = .PRIOR.LOG.ODDS,
  maximum.BF = 10)
{
  motif.seqs <- readr::read_csv(
    seqs.path, skip = 1, progress = FALSE,
    col_names = c('length', 'ID'),
    col_types = cols(length=col_integer(), ID = col_character()))
  message('Loading: ', results.path)
  readr::read_csv(
    results.path, skip = 1, progress = FALSE,
    col_names = c('motif', 'w.mer', 'seq', 'position', 'strand',
                  'Z', 'score', 'p.value'),
    col_types = cols(
      motif = col_character(),
      w.mer = col_character(),
      seq = col_integer(),
      position = col_integer(),
      strand = col_character(),
      Z = col_double(),
      score = col_integer(),
      p.value = col_double())) %>%
    mutate(motif = rify(motif),
           chr = motif.seqs$ID[seq+1],
           end = position + str_length(w.mer),
           logBF = Z.to.log.BF(Z, prior.log.odds, maximum = maximum.BF),
           neg.log.p = -log10(p.value)) %>%
  group_by(motif) %>%
  do(gr = with(.,
    GRanges(
      seqnames = Rle(chr),
      ranges = IRanges(start = position+1, end = end),
      strand = strand,
      seqinfo = seqinfo(hg19),
      Z = Z,
      logBF = logBF,
      neg.log.p = neg.log.p)))
})




.PRIOR <- 1e-3
.PRIOR.LOG.ODDS <- log10(.PRIOR) - log10(1 - .PRIOR)
#' Motif Z to log Bayes factor
#'
Z.to.log.BF <- function(Z, prior.log.odds = .PRIOR.LOG.ODDS, maximum = 10) {
  ifelse(1 == Z, maximum, log10(Z) - log10(1 - Z) - prior.log.odds)
}


#' Construct a run-length encoded vector with zeros everywhere else than specified
#'
#' DEPRECATED: Rle.from.sparse is much more efficient.
#'
sparse.to.rle <- function(len, idxs, x) {
  warning("This function (sparse.to.rle) has been deprecated. Rle.from.sparse is much more efficient.")
  res <- Rle(0, len)
  res[idxs] <- x
  res
}


#' Construct a run-length encoded vector with zeros everywhere else than specified
#'
Rle.from.sparse <- function(len, idxs, x) {
  # The number of non-zero entries
  N <- length(idxs)
  # Check the indexes are sorted
  stopifnot(N < 2 || all(idxs[1:(N-1)] < idxs[2:N]))
  # Will store the values (including the zeros)
  .vals <- rep(0, 2*N+1)
  # Will store the lengths (for the zeros and non-zeros)
  .lens <- rep(1, 2*N+1)
  # Extend the non-zero indices to beginning and end of range
  .idxs.ext <- c(0, idxs, len+1)
  # Calculate the length of the zero runs
  zero.lens <- .idxs.ext[2:(N+2)] - .idxs.ext[1:(N+1)] - 1
  # Where are the zeros in the .vals and .lens vector?
  zero.idxs <- seq.int(1, 2*N+1, 2)
  # Set the lengths of the zeros
  .lens[zero.idxs] <- zero.lens
  # Where are the non-zeros in the .vals and .lens vector?
  non.zero.idxs <- seq.int(2, 2*N,   2)
  # Set the non-zero values
  .vals[non.zero.idxs] <- x
  # Which lengths are positive?
  pos.lens <- .lens > 0
  # Some of the lengths will be zero, ignore them and their values
  Rle(.vals[pos.lens], .lens[pos.lens])
}


#' Summarise motif scan hits on GRanges
#'
summarise.scan.hits <- function(gr, scan.hits) {
  # Get the scores as Rle by range
  with(
    # Find the overlaps between the GRanges and the scan hits
    as.data.frame(findOverlaps(gr, scan.hits, ignore.strand=TRUE)) %>%
      # Attach the scores to the overlaps
      mutate(
        logBF     = scan.hits$logBF[subjectHits],
        neg.log.p = scan.hits$neg.log.p[subjectHits]) %>%
      # Group by range
      group_by(queryHits) %>%
      # Summarise by maximum scores
      summarise(logBF = max(logBF), neg.log.p = max(neg.log.p)),
    # Return a list of scores
    list(
      logBF     = Rle.from.sparse(length(gr), queryHits, logBF),
      neg.log.p = Rle.from.sparse(length(gr), queryHits, neg.log.p)))
}
