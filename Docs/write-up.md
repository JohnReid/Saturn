---
title: "Saturn team write-up"
output: pdf_document
bibliography: ENCODE.bib
---


We used `xgboost` [@chen_xgboost:_2016] to fit a gradient tree boosting model
to a set of features derived from the called DNase peaks and a genome wide scan
for known motifs and *de novo* motifs. We used the DGF method `Wellington` to
create new features of just those binding sites overlapping a footprint. We
found that using a kernel smoother on our predictions increased our
cross-validation performance.


## Features

The most predictive feature we used were the p-values of the relaxed DNase peaks.

We compiled a list of motifs that represented the binding specificities of all
the TFs in the challenge from the following public motif databases:

- Human and Mouse high-throughput SELEX motifs from [@jolma_multiplexed_2010]
- Human and Mouse HT-SELEX motifs from [@jolma_dna-binding_2013]
- The Swiss Regulon human and mouse motifs [@pachkov_swissregulon_2013]
- From [@zhao_quantitative_2011]
- JASPAR CORE [@mathelier_jaspar_2016]
- Direct and inferred motifs for *Homo sapiens* from [@weirauch_determination_2014]

We then perfomed a genome-wide scan of all these motifs using the `STEME`
[@reid_steme:_2014] software to generate a set of putative binding sites. These
were summarised as region-level features by the maximum of their log-odds ratio
(Bayes factor?).

We used a discriminative motif finder, `DREME` [@bailey_dreme:_2011], to find
motifs that discriminated between the bound sequences and those that were
unbound for each TF. In exactly the same way as for the known motifs, we
perfomed a genome-wide scan and summarised the putative binding sites as
region-level features on a per-TF basis.

We used `Wellington` [@piper_wellington:_2013], a DNase footprint detection
algorithm to determine TF binding footprints in the DNase peaks. We used these
to filter both the known motif binding sites and those the *de novo* sites. We
included both the filtered binding sites and the unfiltered sites as features.


## Predictions

We applied a kernel smoothing method to the predictions. We adjusted the log
odds ratio for each region using a Gaussian kernel. We used TF-specific
length-scales between 0 and 200 base pairs that we chose using
cross-validation.


## References
