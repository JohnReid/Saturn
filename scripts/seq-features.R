#!/usr/bin/env Rscript
#
# Creates features from the sequence
#
"Usage:
seq-feat.R [options]

Options:
  --max-K=MAXK              Longest k-mers we should inspect [default: 2]" -> doc


#
# Set warnings as errors
#
options(warn = 2)




# Load packages
#
devtools::load_all()
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19


#
# Parse options
#
.args <- "--max-K=2"
# Use dummy arguments if they exist otherwise use command line arguments
if (! exists(".args")) .args <- commandArgs(TRUE)
message('Arguments: ', do.call(paste, as.list(.args)))
opts <- docopt::docopt(doc, args = .args)
print(opts)
max.k <- as.integer(opts[['max-K']])
message('Maximum K: ', toString(max.k))


k.stats <- function(seq, k) {
  counts <- oligonucleotideFrequency(seq, width = k, as.array = TRUE, simplify.as = 'collapsed')
  dim(counts)
}

#
# Calculate statistics for a region
#
region.stats <- function(region) {
  seq <- getSeq(genome, region)
  for (k in 1:max.k) {
  }

}


stats.for.ranges <- function(ranges, k) {
  #
  # Allocate frequencies matrix
  freqs <- matrix(nrow = length(ranges), ncol = 4**k)
  #
  # Calculate frequencies by chromosome
  for (chr in runValue(chrom(ranges))) {
    print(chr)
    idxs <- chr == chrom(ranges)
    freqs[which(idxs),] <- oligonucleotideFrequency(getSeq(genome, ranges[idxs]), width = k)
  }
  #
  # Make each row's (range's) counts sum to same value and convert to frequencies
  sums <- rowSums(freqs)
  max.sum <- max(sums)
  zero.rows <- sums == 0
  freqs[zero.rows,] = max.sum / ncol(freqs)
  sums[zero.rows] = max.sum
  freqs <- freqs / sums
  .sums <- rowSums(freqs)
  range(.sums)
  freqs[sample(nrow(freqs), 20),]
  #
  # Calculate quantile for each column (i.e. k-mer)
  # and only keep those in first or last 5%
  quantise.k.mer <- function(k.mer.freqs) {
    quant <- quantile(k.mer.freqs, probs = c(.05, .95))
    idxs <- k.mer.freqs < quant[1] | k.mer.freqs > quant[2]
    Matrix::sparseMatrix(
      i = which(idxs),
      j = rep(1, sum(idxs)),
      x = k.mer.freqs[idxs],
      dims = c(nrow(freqs), 1))
  }
  #
  # Bind columns (k-mers) together
  feats.sparse <- do.call(Matrix::cBind, apply(freqs, 2, quantise.k.mer))
}
# ranges <- ranges.test()[chrom(ranges.test()) == 'chr19']
ranges <- ranges.test()
feats <- stats.for.ranges(ranges, k = 1)
