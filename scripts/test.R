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
rest.labels <- load.chip.labels('REST')
object.size(rest.labels)
sample_n(as.data.frame(rest.labels), 8)
ctcf.h1.hesc.peaks <- load.chip.peaks('H1-hESC', 'CTCF', 'conservative')
sample_n(as.data.frame(ctcf.h1.hesc.peaks), 8)


#
# Combine ChIP and DNAse data
h1.hesc.dnase <- summarise.dnase('H1-hESC', 'relaxed')
binding <- binding.tf.cell('CTCF', 'H1_hESC')
cor(as.vector(h1.hesc.dnase), as.integer(binding))
data.frame(binding = binding,
           dnase = h1.hesc.dnase) %>%
  group_by(binding) %>%
  summarise(dnase.mean = mean(dnase),
            dnase.sd = sd(dnase))

#
# Check whether there is any binding where DNase is absent
#
# First check at region level with summaries of DNase
dnase.absent <- h1.hesc.dnase != 0
has.binding <- binding != 1
dnase.absent & has.binding
#
# Now check at base-pair level with DNase bigwig file
dnase.path <- file.path(saturn.data(), 'DNASE', 'fold_coverage_wiggles', 'DNASE.H1-hESC.fc.signal.bigwig')
dnase.bw <- import(dnase.path)
dnase.bw.absent = dnase.bw[dnase.bw$score == 0]
regions.binding <- regions[binding != 1]


#
# Expression data
h1.hesc.expr.1 <- saturn.expr('H1-hESC', 1)
names(h1.hesc.expr.1)
dim(h1.hesc.expr.1)
sample_n(h1.hesc.expr.1, 2)
ggplot(h1.hesc.expr.1, aes(x=TPM, y=FPKM)) + geom_point(alpha=.1) +
  scale_x_log10() +
  scale_y_log10()
ggplot(h1.hesc.expr.1, aes(x=FPKM)) +
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
# Check if labels have same ranges for distinct TFs
all(arid3a.labels == rest.labels)

#
# Load motif scan results
results.path <- file.path(saturn.data(), 'motifs', 'steme-pwm-scan.out')
seqs.path <- file.path(saturn.data(), 'motifs', 'steme-pwm-scan.seqs')
motif.scan <- load.motif.scan(results.path, seqs.path)
motif.scan

#
# Look at how REST hits overlap with ChIP-seq ranges
#
# Just get the REST hits
rest.scan <- (motif.scan %>% filter('REST.p3' == motif))$gr[[1]]
length(rest.scan)
# Merge them
rest.labels <- load.chip.labels('REST')
system.time(rest.labels <- summarise.scan.hits(rest.labels, rest.scan))

#
# Check overlaps
# Calculate which REST hits overlap other REST hits
rest.scan.reduced <- reduce(rest.scan, ignore.strand=TRUE)
rest.scan.overlaps <- rest.scan.reduced[width(rest.scan.reduced) > 22,]
overlapping <- sort(rest.scan[queryHits(findOverlaps(rest.scan, rest.scan.overlaps, ignore.strand=TRUE))])
example <- overlapping[205:207,]
ov <- findOverlaps(rest.labels, example)
regions <- unique(queryHits(ov))
rest.scan[subjectHits(ov[queryHits(ov) %in% regions])]
rest.labels[regions]


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


#
# Test plotting
library(ggbio)
rest.hepg2.labels <- rest.labels[mcols(rest.labels)[,'HepG2'] != 'U']
autoplot(seqinfo(rest.hepg2.labels)[chr.levels]) +
  layout_karyogram(rest.hepg2.labels,
                   aes_string(fill='HepG2', color='HepG2'))
region <- reduce(rest.labels[regions,])
sites <- subsetByOverlaps(rest.scan, region, ignore.strand = TRUE)
autoplot(sites, geom = "bar", alpha = .2, aes(y = logBF))
