#!/usr/bin/env Rscript
"Usage: scores-analyse.R SCORESTSV" -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
devtools::load_all()
library(ggplot2)
library(ggthemes)
library(stringr)


#' Parse TF and cell from filename
#'
# parse.filename <- function(file.name) str_split_fixed(basename(file.name), fixed('.'), 5)[,2:4]
# Some old filenames didn't include the method (glmnet) so parse these differently
parse.filename <- function(file.name) {
  split <- str_split(basename(file.name), fixed('.'))
  .list <- lapply(
    split,
    function(s) {
      if (6 == length(s)) {
        c(s[[3]], s[[4]], s[[5]], s[[2]])
      } else {
        c(s[[2]], s[[3]], s[[4]], 'glmnet')
      }
    })
  do.call(rbind, .list)
}



#
# Functions
#
plots.filename <- function(plot.name) file.path(plots.dir, str_c(plots.tag, '-', plot.name, '.pdf'))
save.plot <- function(plot.name) ggsave(plots.filename(plot.name), width = 18, height = 12, units = 'in')


#
# Parse options
#
# .args <- "../slurm/scores-with-Well.tsv"
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
scorestsv <- opts$SCORESTSV
plots.dir <- '../Plots'
plots.tag <- 'scores'


#
# Load scores
#
message('Loading scores: ', scorestsv)
scores <- readr::read_tsv(scorestsv)
parsed <- parse.filename(scores$file)
scores$TF <- factor(parsed[,1], levels = tf.levels)
scores$cell <- factor(parsed[,2], levels = cell.levels)
scores$motif.tags <- factor(stringr::str_replace(parsed[,3], 'DREME-.*', 'DREME'))
scores$method <- factor(parsed[,4])
sapply(scores, class)


#
# Only use scores with all predictions
#
num.preds <-
  scores %>%
  group_by(TF, cell) %>%
  summarise(num.preds = n())
most.preds = max(num.preds$num.preds)
scores.filtered <- scores %>% left_join(num.preds) %>% filter(most.preds == num.preds)

#
# Plot scores
#
# AUROC vs. AUPRC
ggplot(scores, aes(x = AUROC, y = AUPRC, label = TF, colour = interaction(motif.tags, method))) +
  geom_label() +
  scale_colour_few() +
  theme_few()
save.plot('AUROC-AUPRC')
# AUPRC by TF
ggplot(scores.filtered, aes(x = reorder(TF, AUPRC, FUN = median), y = AUPRC, colour = interaction(motif.tags, method))) +
  geom_boxplot() +
  # geom_jitter(height = 0) +
  labs(x = 'TF') +
  scale_colour_few() +
  theme_few()
save.plot('AUPRC-by-TF')
# AUPRC by cell
ggplot(scores.filtered, aes(x = reorder(cell, AUPRC, FUN = median), y = AUPRC, colour = interaction(motif.tags, method))) +
  geom_boxplot() +
  # geom_jitter(height = 0) +
  labs(x = 'cell') +
  scale_colour_few() +
  theme_few()
save.plot('AUPRC-by-cell')
# recall at 10% FDR by TF
ggplot(scores.filtered, aes(x = reorder(TF, recall_10, FUN = median), y = recall_10, colour = interaction(motif.tags, method))) +
  geom_boxplot() +
  # geom_jitter(height = 0) +
  labs(x = 'TF') +
  scale_colour_few() +
  theme_few()
save.plot('recall-10-by-TF')
# recall at 50% FDR by TF
ggplot(scores.filtered, aes(x = reorder(TF, recall_50, FUN = median), y = recall_50, colour = interaction(motif.tags, method))) +
  geom_boxplot() +
  # geom_jitter(height = 0) +
  labs(x = 'TF') +
  scale_colour_few() +
  theme_few()
save.plot('recall-50-by-TF')
