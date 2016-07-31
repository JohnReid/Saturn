#!/usr/bin/env Rscript
#
# Make features from motif scores
#

devtools::load_all()
library(Saturn)

#
# Set up file paths
#
motifs.dir   <- file.path(saturn.data(), 'motifs')
motifs.out   <- file.path(motifs.dir, 'steme-pwm-scan.out')
motifs.seqs  <- file.path(motifs.dir, 'steme-pwm-scan.seqs')
features.dir <- file.path(saturn.data(), 'Features', 'Motifs', 'Known')
stopifnot(dir.exists(motifs.dir))
stopifnot(dir.exists(features.dir))

#
# Read hits
#
# Read sequence names
seqs <-
  fread(motifs.seqs, col.names = c('length', 'chrom'), header = TRUE) %>%
  mutate(chrom = factor(chrom, chrs.levels))
sapply(seqs, class)
# Read hits
hits <-
  fread(
    motifs.out,
    header = TRUE,
    col.names = c('motif', 'w.mer', 'seq', 'start', 'strand', 'Z', 'score', 'p.value'),
    colClasses=c("character", "character", "integer", "integer",
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
        ranges = IRanges(start = start, width = str_length(w.mer)),
        strand = strand,
        seqinfo = seqinfo(hg19)))
    # Find overlaps with test ranges
    overlaps <- as.data.frame(findOverlaps(ranges.test(), hits.gr, ignore.strand = TRUE))
    overlaps$value <- motif.hits$Z[overlaps$subjectHits]
    with(
      aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
      Rle.from.sparse(length(ranges.test()), queryHits, value))
  })
names(features) <- motifs

#
# Save features and motif names
#
feature.file <- function(m) file.path(features.dir, str_c('motif-', m, '.rds'))
motifs.meta <- data.frame(motif = motifs, motif.r = rify(motifs)) %>%
  mutate(feature.file = feature.file(motif))
saveRDS(motifs.meta, file.path(features.dir, 'motif-names.rds'))
lapply(motifs, function(m) saveRDS(features[[m]], feature.file(m)))
