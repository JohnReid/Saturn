#!/usr/bin/env Rscript
#
# Make DNase summary features
#


#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
devtools::load_all()
library(Saturn)


#
# For each cell type
#
for (cell in levels(cells$cell)) {
  if ("" != cell) {
    message(cell)
    # Load data
    system.time(dnase.summary  <- summarise.dnase(cell, type = 'relaxed'))
    message(cell, ' DNase summary: ', object.size(dnase.summary), ' bytes')
    # Save data
    feature.file <- file.path(saturn.data(), 'Features', 'DNase', str_c('dnase-summary-', cell, '.rds'))
    saveRDS(dnase.summary, feature.file)
  }
}
