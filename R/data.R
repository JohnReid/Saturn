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

tf.levels <- c(
  "ARID3A", "ATF2",   "ATF3",   "ATF7",   "CEBPB",  "CREB1",  "CTCF",   "E2F1",
  "E2F6",   "EGR1",   "EP300",  "FOXA1",  "FOXA2",  "GABPA",  "GATA3",  "HNF4A",
  "JUND",   "MAFK",   "MAX",    "MYC",    "NANOG",  "REST",   "RFX5",   "SPI1",
  "SRF",    "STAT3",  "TAF1",   "TCF12",  "TCF7L2", "TEAD4",  "YY1",    "ZNF143")

cell.levels <- c(
    'A549',
    'GM12878',
    'H1-hESC',
    'HCT116',
    'HeLa-S3',
    'HepG2',
    'IMR-90',
    'induced_pluripotent_stem_cell',
    'K562',
    'liver',
    'MCF-7',
    'Panc1',
    'PC-3',
    'SK-N-SH')

chr.levels <- c(stringr::str_c('chr', 1:22), 'chrX')

binding.levels <- c('U', 'A', 'B')

hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

#
# Regexes to extract TF and cell from filename
#
tf.regex <- stringr::str_c(tf.levels, collapse = '|')
cell.regex <- stringr::str_c(cell.levels, collapse = '|')


# Convert from Rle to one column matrix
#
setAs("Rle", "Matrix", function(from) {
    rv <- runValue(from)
    nz <- rv != 0
    i <- as.integer(ranges(from)[nz])
    x <- rep(rv[nz], runLength(from)[nz])
    Matrix::sparseMatrix(i=i, p=c(0L, length(x)), x=x, dims=c(length(from), 1))
})


# Convert from DataFrame of Rle to sparse Matrix
#
# Conversion from DataFrame of Rle to sparse Matrix
# from: https://support.bioconductor.org/p/66586/#85609
#
setAs("DataFrame", "Matrix", function(from) {
  mat = do.call(cbind, lapply(from, as, "Matrix"))
  colnames(mat) <- colnames(from)
  rownames(mat) <- rownames(from)
  mat
})


#' Logit function
#'
logit <- function(p) log(p) - log(1-p)


#' Inverse of logit function
#'
logit.inv <- function(l) {
  e <- exp(l)
  e / (e + 1)
}


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
rify <- function(x) stringr::str_replace_all(x, '[-{}]', '_')


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
    col_types = readr::cols(
      chrom = readr::col_factor(seqnames(hg19)),
      chromStart = readr::col_integer(),
      chromEnd = readr::col_integer(),
      name = readr::col_character(),
      score = readr::col_integer(),
      strand = readr::col_character(),
      signalValue = readr::col_double(),
      pValue = readr::col_double(),
      qValue = readr::col_double(),
      peak = readr::col_integer()))
}

#' Convert a narrowPeak data frame to GRanges.
narrowpeak.granges <- function(narrowpeak) with(narrowpeak,
  GRanges(
    seqnames = S4Vectors::Rle(chrom),
    ranges = IRanges(start = chromStart+1, end = chromEnd),
    seqinfo = seqinfo(hg19),
    signal = signalValue,
    pValue = pValue,
    qValue = qValue,
    peak = peak))


#' Convert binding factor to numeric
binding.as.numeric <- function(binding) ifelse('B' == binding, 1, ifelse('A' == binding, .5, 0))


#' Convert a labels data frame to GRanges.
labels.granges <- function(labels) {
  mcol.names <- names(labels)[4:ncol(labels)]
  .args = lapply(labels[,mcol.names], function(x) S4Vectors::Rle(as.integer(x)))
  .args$seqnames = S4Vectors::Rle(labels$chr)
  .args$ranges = IRanges(start = labels$start+1, end = labels$stop)
  .args$seqinfo = seqinfo(hg19)
  do.call(GRanges, args = .args)
}


#' Load regions from annotations
#'
load.regions <- function(name) {
    file <- file.path(saturn.data(), 'annotations', stringr::str_c(name, '_regions.blacklistfiltered.bed.gz'))
    # sort(import(file, seqinfo=seqinfo(hg19)))
    df <-
       readr::read_tsv(file, col_names = c('chrom', 'start', 'end'), progress = FALSE) %>%
       select(-end) %>%
       mutate(chrom = factor(chrom, levels = chr.levels)) %>%
       arrange(chrom, start)
    S4Vectors::DataFrame(chrom = S4Vectors::Rle(df$chrom), start = df$start)
}


#' Assume all columns that are not named cell, chrom, start, stop, end are cell names
#'
cell.names <- function(df) colnames(df)[! colnames(df) %in% c('cell', 'chrom', 'start', 'end', 'stop')]


#' The index of the region in the test regions, given its index in the training regions
#'
train.to.test.idx <- function(idx) regions.train$test.idx[idx]


#' The index of the region in the test regions, given its index in the ladder regions
#'
ladder.to.test.idx <- function(idx) regions.ladder$test.idx[idx]


#' The index of the region in the training regions, given its index in the test regions
#'
test.to.train.idx <- function(idx) regions.test$train.idx[idx]


#' The index of the region in the ladder regions, given its index in the test regions
#'
test.to.ladder.idx <- function(idx) regions.test$ladder.idx[idx]


#' The test indexes of the training regions
#'
training.region.test.idxs <- memoise::memoise(function() ! Rle(is.na(regions.test$train.idx)))


#' The test indexes of the ladder regions
#'
ladder.region.test.idxs <- function() ! training.region.test.idxs()


#' Convert a factor-Rle that is indexed by the training ranges to a factor-Rle
#' that is indexed by the test ranges, with NA filling in the missing values
#'
train.factor.Rle.to.test <- function(x) {
  res <- Rle(factor(NA, levels = levels(x)), nrow(regions.test))
  res[regions.train$test.idx] <- x
  res
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


#' Construct a run-length encoded vector with zeros everywhere else than specified
#'
#' DEPRECATED: Rle.from.sparse is much more efficient.
#'
sparse.to.rle <- function(len, idxs, x) {
  warning("This function (sparse.to.rle) has been deprecated. Rle.from.sparse is much more efficient.")
  res <- S4Vectors::Rle(0, len)
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
  S4Vectors::Rle(.vals[pos.lens], .lens[pos.lens])
}


#' Parse predictions file name
#'
parse.predictions.file <- function(file.name) {
  fields <- stringr::str_split_fixed(basename(file.name), stringr::fixed('.'), 5)
  list(
    tf = fields[2],
    cell = fields[3],
    motif.tags = stringr::str_split(fields[4], stringr::fixed('_')))
}
