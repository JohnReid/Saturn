.PRIOR <- 1e-3
.PRIOR.LOG.ODDS <- log10(.PRIOR) - log10(1 - .PRIOR)
.MAX.LOG.BF <- 10


#' Motif Z to log Bayes factor
#'
Z.to.log.BF <- function(Z, prior.log.odds = .PRIOR.LOG.ODDS, maximum = .MAX.LOG.BF) {
  ifelse(1 == Z, maximum, log10(Z) - log10(1 - Z) - prior.log.odds)
}


#' Load motif scan results from directory
#'
load.motif.dir <- function(
  motifs.dir,
  prior.log.odds = .PRIOR.LOG.ODDS,
  maximum.BF = 10)
{
  load.motif.scan(
    results.path = file.path(motifs.dir, 'steme-pwm-scan.out'),
    seqs.path = file.path(motifs.dir, 'steme-pwm-scan.seqs'),
    prior.log.odds,
    maximum.BF)
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
    col_types = readr::cols(length = readr::col_integer(), ID = readr::col_character()))
  message('Loading: ', results.path)
  motif.scan <-
    data.table::fread(
      results.path,
      sep = ',',
      showProgress = FALSE,
      verbose = FALSE,
      col.names = c('motif', 'w.mer', 'seq', 'position', 'strand',
                    'Z', 'score', 'p.value'),
      colClasses = c(
        'character',
        'character',
        'integer',
        'integer',
        'integer',
        'numeric',
        'integer',
        'numeric')) %>%
    mutate(
      motif = factor(rify(motif)),
      chr = motif.seqs$ID[seq+1],
      end = position + stringr::str_length(w.mer),
      logBF = Z.to.log.BF(Z, prior.log.odds, maximum = maximum.BF),
      neg.log.p = -log10(p.value))
  motif.names <- levels(motif.scan$motif)
  scans <- lapply(
    motif.names,
    function(m) with(
      filter(motif.scan, motif == m),
      sort(GRanges(
        seqnames = S4Vectors::Rle(chr),
        ranges = IRanges(start = position + 1, end = end),
        strand = strand,
        seqinfo = seqinfo(hg19),
        Z = Z,
        logBF = logBF,
        neg.log.p = neg.log.p))))
  names(scans) <- motif.names
  scans
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
