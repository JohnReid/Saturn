#!/usr/bin/env Rscript

devtools::load_all()

#
# Convert the known motifs
#
known.feature.dir <- file.path(features.dir(), 'Motifs', 'Known')
saveRDS(
  load.motif.features(known.feature.dir),
  feature.file.name('KnownMotifs'))

#
# Convert the DREME motifs
#
for (tf in tfs$TF) {
  dreme.feature.dir <- file.path(features.dir(), 'Motifs', stringr::str_c('DREME-', tf))
  saveRDS(
    load.motif.features(dreme.feature.dir),
    feature.tf.file.name('DREME', tf))
}
