#!/usr/bin/env Rscript
TF <- commandArgs()[2]
devtools::load_all()
rmarkdown::render(
    "get-TF-fasta.Rmd",
    output_file = str_c('FASTA-', TF, '.html'),
    params = list(TF = TF))
