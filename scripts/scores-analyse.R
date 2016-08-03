#!/usr/bin/env Rscript
"Usage: scores-analyse.R SCORESTSV" -> doc


#
# Load libraries
#
library(ggplot2)
library(stringr)


parse.filename <- function(file.name) {
  str_split_fixed(basename(file.name), fixed('.'), 4)[,2:3]
}


#
# Parse options
#
# .args <- "../slurm/my-scores.tsv"
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
scores


#
# Plot scores
#
ggplot(scores, aes(x=AUROC, y=AUPRC, colour=TF)) + geom_point()
