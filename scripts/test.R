library(ggplot2)
library(dplyr)
library(stringr)
library(GenomicRanges)
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
