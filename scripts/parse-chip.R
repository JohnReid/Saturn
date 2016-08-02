#!/usr/bin/Rscript
#
# Parse ChIP label files
#

devtools::load_all()
library(Saturn)


#
# For each TF
#
for (tf in tf.levels) {
  tf.cells <-
    tfs %>%
    filter(TF == tf, 'train' == split) %>%
    mutate(feature.file = chip.feature.file(tf, cell))
  #
  # If all ChIP features exist then don't recreate them
  #
  if (all(file.exists(tf.cells$feature.file))) {
    message('Already have label features for ', tf)
    next
  }
  #
  # Load the ChIP data
  #
  chip <- read.chip.labels(tf)
  # Check we have the correct columns, i.e. they match the rows in tf.cells
  stopifnot(all(sort(colnames(chip)[3:ncol(chip)]) == sort(as.character(tf.cells$cell))))
  #
  # Write each feature separately
  #
  message('Writing ', tf, ' features')
  invisible(tf.cells %>%
    rowwise() %>%
    do(written = saveRDS(chip[.$cell], .$feature.file)))
}
