#!/usr/bin/env Rscript

devtools::load_all()

#
# Convert the known motifs
#
known.feature.dir <- file.path(features.dir(), 'Motifs', 'Known')
new.feature.file <- feature.file.name('KnownMotifs')
if (! file.exists(new.feature.file)) {
  saveRDS(
    load.motif.features(known.feature.dir),
    new.feature.file)
}

#
# Convert the DREME motifs
#
for (tf in tf.levels) {
  if ('E2F1' == tf) next
  dreme.feature.dir <- file.path(features.dir(), 'Motifs', stringr::str_c('DREME-', tf))
  new.feature.file <- feature.tf.file.name('DREME', tf)
  if (! file.exists(new.feature.file)) {
    saveRDS(
      load.motif.features(dreme.feature.dir),
      new.feature.file)
  }
}
