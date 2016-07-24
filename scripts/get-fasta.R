#!/usr/bin/env Rscript
devtools::load_all()
for (TF in levels(tfs$TF)) {
    rmarkdown::render(
        "get-TF-fasta.Rmd",
        output_file = str_c('FASTA-', TF, '.html'),
        params = list(TF = TF))
}
