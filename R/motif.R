.PRIOR <- 1e-3
.PRIOR.LOG.ODDS <- log10(.PRIOR) - log10(1 - .PRIOR)
.MAX.LOG.BF <- 10


#' Motif Z to log Bayes factor
#'
Z.to.log.BF <- function(Z, prior.log.odds = .PRIOR.LOG.ODDS, maximum = .MAX.LOG.BF) {
  ifelse(1 == Z, maximum, log10(Z) - log10(1 - Z) - prior.log.odds)
}


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
           end = position + stringr::str_length(w.mer),
           logBF = Z.to.log.BF(Z, prior.log.odds, maximum = maximum.BF),
           neg.log.p = -log10(p.value)) %>%
  group_by(motif) %>%
  do(gr = with(.,
    GRanges(
      seqnames = S4Vectors::Rle(chr),
      ranges = IRanges(start = position+1, end = end),
      strand = strand,
      seqinfo = seqinfo(hg19),
      Z = Z,
      logBF = logBF,
      neg.log.p = neg.log.p)))
})


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
