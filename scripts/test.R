library(ggplot2)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
options(dplyr.width = Inf)

devtools::load_all('.')
print(saturn.data())

#
# DNase data
liver.dnase <- load.dnase.peaks('liver', 'conservative')
dnase <- as.data.frame(liver.dnase)
sample_n(dnase, 8)
dnase <- dnase %>%
  mutate(len=end-start,
         peakBias=peak/len)
sample_n(dnase, 8)
ggplot(dnase, aes(x=peakBias)) + geom_histogram()
ggplot(dnase, aes(x=len)) + geom_histogram()
ggplot(dnase, aes(x=len)) + geom_histogram() + scale_x_log10()


#
# ChIP data
myc.labels <- load.chip.labels('MYC')
sample_n(as.data.frame(myc.labels), 8)
ctcf.h1.hesc.peaks <- load.chip.peaks('H1-hESC', 'CTCF', 'conservative')
sample_n(as.data.frame(ctcf.h1.hesc.peaks), 8)


#
# Combine ChIP and DNAse data
K562.dnase <- load.dnase.peaks('K562', 'conservative')
myc.K562.dnase <- combine.chip.dnase(myc.labels, K562.dnase)
binding <- binding.as.numeric(myc.labels$mcols.K562)
with(myc.K562.dnase, cor(binding, dnase))
data.frame(binding = myc.K562.dnase$mcols.K562,
           dnase = myc.K562.dnase$dnase) %>%
  group_by(binding) %>%
  summarise(dnase.mean = mean(dnase),
            dnase.sd = sd(dnase))


#
# Expression data
K562.expr.1 <- saturn.expr('K562', 1)
names(K562.expr.1)
dim(K562.expr.1)
sample_n(K562.expr.1, 2)
ggplot(K562.expr.1, aes(x=TPM, y=FPKM)) + geom_point(alpha=.1) +
  scale_x_log10() +
  scale_y_log10()
ggplot(K562.expr.1, aes(x=FPKM)) +
  geom_histogram() +
  geom_rug(alpha=.1) +
  scale_x_log10()


#
# Examine available data
#
# For example, ARID3A in K562 is in the ladderboard split
arid3a.labels <- load.chip.labels('ARID3A')
arid3a.labels
#
# So there is only binding data for HepG2 not K562.
#
# Let's check if we have chromosomes 1, 21 and 8
'chr1' %in% seqlevels(arid3a.labels)
'chr2' %in% seqlevels(arid3a.labels)
'chr8' %in% seqlevels(arid3a.labels)
'chr21' %in% seqlevels(arid3a.labels)
# No


#
# Load motif scan results
results.path <- file.path(saturn.data(), 'motifs', 'steme-pwm-scan.out')
seqs.path <- file.path(saturn.data(), 'motifs', 'steme-pwm-scan.seqs')
motif.seqs <- readr::read_csv(seqs.path)
prior.log.odds <- .PRIOR.LOG.ODDS
hits <- readr::read_csv(
    results.path, skip = 1, col_types = 'cciicnnn', progress = FALSE,
    col_names = c('motif', 'w.mer', 'seq', 'position', 'strand',
                  'Z', 'score', 'p.value')) %>%
  mutate(chr = motif.seqs$ID[seq+1],
         end = position + str_length(w.mer),
         logBF = Z.to.log.BF(Z, prior.log.odds),
         neg.log.p = -log10(p.value))
hits.gr <- hits %>%
  group_by(motif) %>%
  do(ranges = with(.,
    GRanges(
      seqnames = Rle(chr),
      ranges = IRanges(start = position+1, end = end),
      strand = strand,
      seqinfo = seqinfo(hg19),
      motif = motif,
      Z = Z,
      logBF = logBF,
      neg.log.p = neg.log.p)))
hits.gr[hits.gr$motif=='REST.p3',2][[1]]

motif.scan <- load.motif.scan(
  file.path(saturn.data(), 'motifs', 'steme-pwm-scan.out'),
  file.path(saturn.data(), 'motifs', 'steme-pwm-scan.seqs'))
motif.scan

#
# Look at how REST hits overlap with ChIP-seq ranges
#
# Just get the REST hits
rest.scan <- (motif.scan %>% filter('REST.p3' == motif))$gr[[1]]
length(rest.scan)
# Merge them
rest.labels <- load.chip.labels('REST')
system.time(rest.labels <- merge.scan.hits(rest.labels, rest.scan))

#
# Check overlaps
# Calculate which REST hits overlap other REST hits
rest.scan.reduced <- reduce(rest.scan, ignore.strand=TRUE)
rest.scan.overlaps <- rest.scan.reduced[width(rest.scan.reduced) > 22,]
overlapping <- sort(rest.scan[queryHits(findOverlaps(rest.scan, rest.scan.overlaps, ignore.strand=TRUE))])
example <- overlapping[205:207,]
regions <- unique(queryHits(findOverlaps(myc.labels, example)))
rest.scan[ov[queryHits(ov) %in% regions]@subjectHits]
myc.labels[regions]
myc.logBF[regions]


#
# Test BSgenome.Hsapiens.UCSC.hg19
hg19 <- BSgenome.Hsapiens.UCSC.hg19
seqinfo(hg19)
starts <- c(
  40715041,
  5014919,
  24662264)
gr <- GRanges(
  seqnames = c('chr17', 'chr17', 'chr14'),
  ranges = IRanges(start = starts+1, width = 14),
  seqinfo = seqinfo(hg19))
seqs <- getSeq(hg19, gr)
seqs[[1]]
seqs[[2]]
seqs[[3]][1:4]

