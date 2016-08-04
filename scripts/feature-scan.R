#!/usr/bin/env Rscript
#
# Make features from motif scan
#
"Usage: feature-scan.R SCANDIR SCANTAG" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load libraries
#
devtools::load_all()
library(Saturn)


#
# Parse options
#
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
scan.dir <- opts[['SCANDIR']]
scan.tag <- opts[['SCANTAG']]


#
# Set up file paths
#
scan.out   <- file.path(scan.dir, 'steme-pwm-scan.out')
scan.seqs  <- file.path(scan.dir, 'steme-pwm-scan.seqs')
features.dir <- file.path(saturn.data(), 'Features', 'Motifs', scan.tag)
stopifnot(dir.exists(scan.dir))
stopifnot(file.exists(scan.out))
stopifnot(file.exists(scan.seqs))
stopifnot(! dir.exists(features.dir))
dir.create(features.dir)


#
# Read hits
#
# Read sequence names
message('Read sequence names: ', scan.seqs)
seqs <-
  data.table::fread(scan.seqs, col.names = c('length', 'chrom'), header = TRUE) %>%
  mutate(chrom = factor(chrom, chr.levels))
sapply(seqs, class)
# Read hits
message('Read scan hits: ', scan.out)
hits <-
  data.table::fread(
    scan.out,
    header = TRUE,
    col.names = c('motif', 'w.mer', 'seq', 'start', 'strand', 'Z', 'score', 'p.value'),
    colClasses = c("character", "character", "integer", "integer",
                   "character", "numeric", "integer", "numeric")) %>%
  mutate(
    motif = factor(motif),
    chrom = seqs$chrom[seq + 1],  # Use 1-based indexing
    start = start + 1,            # Change to 1-based indexing
    strand = factor(strand)) %>%
  select(-seq)
sapply(hits, class)
object.size(hits)


#
# Generate feature for each motif
#
message('Generate features for motifs')
motifs <- levels(hits$motif)
features <- lapply(motifs,
  function(m) {
    message('Motif: ', m)
    # Separate into motifs
    motif.hits <- hits %>% filter(motif == m) %>% select(-motif)
    # Convert to GRanges
    hits.gr <- with(motif.hits,
      GRanges(
        seqnames = chrom,
        ranges = IRanges(start = start, width = stringr::str_length(w.mer)),
        strand = strand,
        seqinfo = seqinfo(hg19)))
    # Find overlaps with test ranges
    overlaps <- as.data.frame(findOverlaps(ranges.test(), hits.gr, ignore.strand = TRUE))
    overlaps$value <- Z.to.log.BF(motif.hits$Z[overlaps$subjectHits])
    with(
      aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
      Rle.from.sparse(length(ranges.test()), queryHits, value))
  })
names(features) <- motifs

#
# Save features and motif names
#
message('Saving features: ', features.dir)
feature.file <- function(m) stringr::str_c('motif-', m, '.rds')
feature.path <- function(m) file.path(features.dir, feature.file(m))
motifs.meta <- data.frame(motif = motifs, motif.r = rify(motifs)) %>%
  mutate(feature.file = feature.file(motif))
saveRDS(motifs.meta, file.path(features.dir, 'motif-names.rds'))
lapply(motifs, function(m) saveRDS(features[[m]], feature.path(m)))
