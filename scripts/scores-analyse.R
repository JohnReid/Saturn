#!/usr/bin/env Rscript
"Usage: scores-analyse.R SCORESTSV" -> doc


#
# Load libraries
#
library(ggplot2)
library(stringr)


#' Parse TF and cell from filename
#'
parse.filename <- function(file.name) str_split_fixed(basename(file.name), fixed('.'), 4)[,2:3]


#
# Parse options
#
# .args <- "my-scores.tsv"
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
scorestsv <- opts$SCORESTSV


#
# Load scores
#
message('Loading scores: ', scorestsv)
scores <- readr::read_tsv(scorestsv)
parsed <- parse.filename(scores$file)
scores$TF <- factor(parsed[,1])
scores$cell <- factor(parsed[,2])


#
# Plot scores
#
# AUROC vs. AUPRC
ggplot(scores, aes(x=AUROC, y=AUPRC, label=TF)) + geom_label()
# AUPRC by TF
ggplot(scores, aes(x=reorder(TF, AUPRC, FUN=median), y=AUPRC)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'TF')
# AUPRC by cell
ggplot(scores, aes(x=reorder(cell, AUPRC, FUN=median), y=AUPRC)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'cell')
# recall 10 by TF
ggplot(scores, aes(x=reorder(TF, recall_10, FUN=median), y=recall_10)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'TF')
# recall 10 by cell
ggplot(scores, aes(x=reorder(cell, recall_10, FUN=median), y=recall_10)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'cell')
# recall 50 by TF
ggplot(scores, aes(x=reorder(TF, recall_50, FUN=median), y=recall_50)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'TF')
# recall 50 by cell
ggplot(scores, aes(x=reorder(cell, recall_50, FUN=median), y=recall_50)) + geom_boxplot() + geom_jitter(height = 0) + labs(x = 'cell')
