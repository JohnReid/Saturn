
# Data pre-processing

## DNase

### BAMs

Summarise the DNase `bigwig` files per region:

    submit-dnase-summarise

Download the BAM files from synapse:

    download-bams

Index the downloaded BAM files:

    bams-list | submit-index

Merge the technical replicates:

    submit-dnase-merge

Run Wellington on the merged BAMs:

    submit-wellington

Link expression of TF to footprint quality (see Gusmao et al. 2016 "Association of TF expression with footprint quality")


## Annotations

### Genome

- strip
- build suffix array


### Regions

Define/retrieve positive/negative sequences

    submit-get-fasta


### Motifs

- Gather motifs from 3rd party libraries

Run DREME on positive and negative sequences:

    submit-dreme

Scan genome with DREME motifs:

    submit-scan-genome

Make motif hits into per-region features:

    submit-featurise-dreme-scan
