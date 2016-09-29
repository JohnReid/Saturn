# Saturn team method write-up

We used `xgboost` to fit a gradient tree boosting model to hand-crafted features.

## Features

### Known motifs

We compiled a list of motifs that represented the binding specificities of all the TFs in the
challenge from the following public motif databases:

- blah
- blah

We then perfomed a genome-wide scan of all these motifs using the STEME software to generate
a set of putative binding sites. These were summarised as region-level features by the maximum
of their log-odds ratio (Bayes factor?).


### *De novo* motifs

We used a discriminative motif finder, DREME, to find motifs that discriminated between
the bound sequences and those that were unbound for each TF. In exactly the same way as for
the known motifs, we perfomed a genome-wide scan and summarised the putative binding sites as
region-level features on a per-TF basis.


### DNase data

### Combining DNase with motifs

We used Wellington, a DNase footprint detection algorithm to determine TF binding footprints
in the DNase peaks. We used these

## Smoothing

We applied a kernel smoothing method to the predictions.
