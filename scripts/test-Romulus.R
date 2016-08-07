library(Romulus)

# Clip the DNase-seq data for NRSF at 99.9% quantile
thresh <- max(quantile(NRSF.cuts, 0.999), 1L)
cuts <- pmin(NRSF.cuts, thresh)

# Forward strand cuts only upstream and within the candidate binding site
cuts1 <- cuts[, 1:(ncol(cuts) / 2 - NRSF.margin)]
# Reverse strand cuts only within the candidate binding site and downstream
cuts2 <- cuts[, ncol(cuts) / 2 + (NRSF.margin + 1):(ncol(cuts) / 2)]

# Take 20 bp bins outside the candidate binding site and 1 bp bins within it
NRSF.width <- (ncol(cuts) - 4 * NRSF.margin) / 2
bins1 <- c(rep(1:10, each = 20), 10 + 1:NRSF.width)
bins2 <- c(1:NRSF.width, NRSF.width + rep(1:10, each = 20))

# Fit the Romulus model
r.fit <- fitRomulus(cuts1, cuts2, NRSF.anno, list(c("score")), bins1, bins2)
