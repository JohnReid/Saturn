library(ggplot2)
library(dplyr)
library(GenomicRanges)
options(dplyr.width = Inf)

devtools::load_all('.')
print(saturn.data())

#
# DNase data
liver.dnase <- load.dnase.peaks('liver', 'conservative')
sample_n(liver.dnase, 8)
liver.dnase <- liver.dnase %>%
  mutate(len=chromEnd-chromStart,
         peakBias=peak/len)
sample_n(liver.dnase, 8)
ggplot(liver.dnase, aes(x=peakBias)) + geom_histogram()
ggplot(liver.dnase, aes(x=len)) + geom_histogram()
ggplot(liver.dnase, aes(x=len)) + geom_histogram() + scale_x_log10()
gr <- narrowpeak.granges(liver.dnase)
gr

#
# ChIP data
myc.labels <- load.chip.labels('MYC')
class(myc.labels)
sapply(myc.labels, class)
dim(myc.labels)
sample_n(myc.labels, 8)
ctcf.h1.hesc.peaks <- load.chip.peaks('H1-hESC', 'CTCF', 'conservative')
sample_n(ctcf.h1.hesc.peaks, 8)


#
# Combine ChIP and DNAse data
K562.dnase <- load.dnase.peaks('K562', 'conservative')
K562.dnase.gr <- narrowpeak.granges(K562.dnase)
myc.gr <- labels.granges(myc.labels)
overlaps <- as.data.frame(findOverlaps(K562.dnase.gr, myc.gr))
overlaps$dnase <- K562.dnase.gr[overlaps$queryHits,]$pValue
overlaps %>%
  group_by(subjectHits) %>%
  summarise(dnase=max(dnase))
dnase <- rep(0, length(myc.gr))
dnase[overlaps$subjectHits] <- overlaps$dnase
qplot(dnase)
myc.gr$k562.dnase <- dnase
binding <- binding.as.numeric(myc.gr$mcols.K562)
with(myc.gr, cor(binding, k562.dnase))
with(myc.gr, qplot(x=mcols.K562, y=myc.gr$k562.dnase) + geom_boxplot())


#
# Expression data
K562.1 <- saturn.expr('K562', 1)
names(K562.1)
dim(K562.1)
sample_n(K562.1, 2)

ggplot(K562.1, aes(x=TPM, y=FPKM)) + geom_point(alpha=.1) +
  scale_x_log10() +
  scale_y_log10()

ggplot(K562.1, aes(x=FPKM)) +
  geom_histogram() +
  geom_rug(alpha=.1) +
  scale_x_log10()

#
# Examine available data
#
# For example, ARID3A in K562 is in the ladderboard split
arid3a.labels <- load.chip.labels('ARID3A')
names(arid3a.labels)
#
# So there is only binding data for HepG2 not K562.
#
# Let's check if we have chromosomes 1, 21 and 8
'chr1' %in% arid3a.labels$chr
'chr2' %in% arid3a.labels$chr
'chr8' %in% arid3a.labels$chr
'chr21' %in% arid3a.labels$chr
#
# No
