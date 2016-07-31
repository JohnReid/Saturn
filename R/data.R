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

cell.levels <- c(
    'A549',
    'GM12878',
    'H1-hESC',
    'HCT116',
    'HeLa-S3',
    'HepG2',
    'IMR90',
    'induced_pluripotent_stem_cell',
    'K562',
    'liver',
    'MCF-7',
    'Panc1',
    'PC-3',
    'SK-N-SH')
chrs.levels <- c(str_c('chr', 1:22), 'chrX')
binding.levels <- c('U', 'A', 'B')

hg19 <- BSgenome.Hsapiens.UCSC.hg19


#' Convert from Rle to one column matrix
#'
setAs("Rle", "Matrix", function(from) {
    rv <- runValue(from)
    nz <- rv != 0
    i <- as.integer(ranges(from)[nz])
    x <- rep(rv[nz], runLength(from)[nz])
    sparseMatrix(i=i, p=c(0L, length(x)), x=x, dims=c(length(from), 1))
})


#' Convert from DataFrame of Rle to sparse Matrix
#'
#' Conversion from DataFrame of Rle to sparse Matrix
#' from: https://support.bioconductor.org/p/66586/#85609
#'
setAs("DataFrame", "Matrix", function(from) {
  mat = do.call(cbind, lapply(from, as, "Matrix"))
  colnames(mat) <- colnames(from)
  rownames(mat) <- rownames(from)
  mat
})


#' Retrieve N largest objects in global environment.
#'
largest.objects <- function(N = 10) {
  z <- sapply(ls(envir = globalenv()), function(x) object.size(get(x, envir = globalenv())))
  as.matrix(rev(sort(z))[1:N])
}


#' Make names suitable for R.
#'
#' E.g. substitute '_' for '-'
#'
rify <- function(x) str_replace_all(x, '[-{}]', '_')


#' Which string matches
str.match <- function(strs, x) which(strs == x)

#' Convert back from R names to original cell names
#'
#' E.g. substitute '-' for '_'
#'
unrify <- function(x) {
  idxs <- vapply(x, function(c.r) str.match(cells$cell.r, c.r), FUN.VALUE=0)
  cells$cell[idxs]
}


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


#' Convert a ChIP-seq labels data frame to GRanges.
labels.granges <- function(labels) {
  mcol.names <- names(labels)[4:ncol(labels)]
  .args = lapply(labels[,mcol.names], function(x) Rle(as.integer(x)))
  .args$seqnames = Rle(labels$chr)
  .args$ranges = IRanges(start = labels$start+1, end = labels$stop)
  .args$seqinfo = seqinfo(hg19)
  do.call(GRanges, args = .args)
}


#' The file with the TF's binding labels
#'
tf.chip.labels.file <- function(tf)
  file.path(saturn.data(), 'ChIPseq', 'labels', str_c(tf, '.train.labels.tsv.gz'))

#' Read the ChIP binding labels for the TF into a S4Vectors::DataFrame
#'
# read.chip.labels <- memoise::memoise(function(tf) {
read.chip.labels <- (function(tf) {
  path <- tf.chip.labels.file(tf)
  message('Loading: ', path)
  # Read ChIP labels into data.table
  labels.dt <- fread(str_c('zcat ', path), sep = '\t', drop = 'stop', verbose = FALSE) %>% rename(chrom = chr)
  # Make chrom into a factor with correct levels
  labels.dt$chrom <- factor(labels.dt$chrom, levels = chrs.levels)
  # Sort by chrom then start
  setkey(labels.dt, chrom, start)
  # Convert into S4Vectors::Dataframe using Rle
  # First convert binding factors
  label.cells <- colnames(labels.dt)[3:ncol(labels.dt)]
  binding.rles <- lapply(label.cells, function(l) Rle(factor(labels.dt[[l]], levels=binding.levels)))
  names(binding.rles) <- label.cells
  # Now convert chromosomes
  df.args <- c(list(chrom=Rle(labels.dt$chrom), start=labels.dt$start, check.names = FALSE), binding.rles)
  # Check the labels are ordered exactly the same as the test regions
  stopifnot(all(df.args$chrom == regions.train$chrom))
  stopifnot(all(df.args$start == regions.train$start))
  # Build DataFrame
  do.call(DataFrame, df.args)
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
              function(col) Rle(factor(col, levels = cell.levels), nrow(chip))))
  # Create DataFrame
  DataFrame(
    cell  = cells.m,
    chrom = rep(chip$chrom, ncol(chip) - 2),
    start = rep(chip$start, ncol(chip) - 2),
    bound = binding.m)
}

#' Convert melted ChIP data into data.table
#'
chip.data.table <- function(chip.m) setkey(
  data.table(
    cell  = as.factor(chip.m$cell),
    chrom = as.factor(chip.m$chrom),
    start = as.vector(chip.m$start),
    bound = as.factor(chip.m$bound)),
  cell, chrom, start)


#' Read ChIP-seq labels into data frame
read.chip.labels.old <- function(tf) {
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



#' Load ChIP peaks
load.chip.peaks <- memoise::memoise(function(cell, tf, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'ChIPseq', 'peaks', type,
              str_c('ChIPseq.', cell, '.', tf, '.', .type, '.narrowPeak.gz'))))
})


#' Load DNase peaks
load.dnase.peaks <- memoise::memoise(function(cell, type='conservative') {
  if ('conservative' == type) .type <- 'conservative.train'
  else .type <- type
  narrowpeak.granges(load.narrowpeak(
    file.path(saturn.data(), 'DNASE', 'peaks', type,
              str_c('DNASE.', cell, '.', type, '.narrowPeak.gz'))))
})


#' Load regions from annotations
#'
load.regions <- function(name) {
    file <- file.path(saturn.data(), 'annotations', str_c(name, '_regions.blacklistfiltered.bed.gz'))
    # sort(import(file, seqinfo=seqinfo(hg19)))
    df <-
       readr::read_tsv(file, col_names = c('chrom', 'start', 'end'), progress = FALSE) %>%
       select(-end) %>%
       mutate(chrom = factor(chrom, levels = chrs.levels)) %>%
       arrange(chrom, start)
    S4Vectors::DataFrame(chrom = Rle(df$chrom), start = df$start)
}


#' Convert regions to ranges
#'
regions.to.ranges <- function(regions) GRanges(
  GRanges(
    seqnames = regions$chrom,
    ranges = IRanges(start = regions$start, width = 200),
    strand = '*',
    seqinfo = seqinfo(hg19)))


#' The ladderboard ranges
#'
ranges.ladder <- memoise::memoise(function() regions.to.ranges(regions.ladder))


#' The test ranges
#'
ranges.test <- memoise::memoise(function() regions.to.ranges(regions.test))


#' The training ranges
#'
ranges.train <- memoise::memoise(function() regions.to.ranges(regions.train))


#' Summarise DNase peaks by ChIP-regions
#'
summarise.dnase <- function(cell, type='conservative', aggregation.fn=max) {
  # Load the peaks
  dnase <- load.dnase.peaks(cell, type)
  # Aggregate
  agg.by.region(dnase, ranges.test(), 'pValue', aggregation.fn)
}

#' Aggregate the named values in gr by region using the aggregation function.
#'
agg.by.region <- function(gr, regions, value.col, aggregation.fn=max) {
  # Find the overlaps between the data and the regions
  overlaps <- as.data.frame(findOverlaps(gr, regions, ignore.strand = TRUE))
  # Add the values to the overlaps
  overlaps$value <- mcols(gr)[overlaps$queryHits, value.col]
  # Summarise the overlaps by the aggregated value for each region
  region.values <- overlaps %>%
    group_by(subjectHits) %>%
    summarise(value=aggregation.fn(value))
  # Return the result as a Rle vector
  Rle.from.sparse(length(regions), region.values$subjectHits, region.values$value)
}

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
