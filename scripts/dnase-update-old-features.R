#!/usr/bin/env Rscript

devtools::load_all()

for (cell in cell.levels) {
  message(cell)
  dnase.levels <- load.dnase.feature(cell)
  saveRDS(S4Vectors::DataFrame(DNase = dnase.levels), feature.cell.file.name('DNase', cell))
}
