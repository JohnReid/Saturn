#!/usr/bin/env Rscript
devtools::load_all()
for (TF in levels(tfs$TF)) {
    rmarkdown::render(
        "get-TF-fasta.Rmd",
        output_file = str_c('FASTA-', TF, '.pdf'),
        params = list(TF = TF))
}
