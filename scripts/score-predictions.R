#!/usr/bin/env Rscript
"Usage: score-predictions.R PREDICTIONSFILE" -> doc


#
# Load libraries
#


#
# Parse options
#
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
predictions.path <- opts$PREDICTIONSFILE


#
# Load predictions
#
predictions <- readr::read_tsv(predictions.path, col_names = c('chrom', 'start', 'end', 'pred', 'bound'))


#
# Assess quality of predictions
#
labels <- as.vector(predictions$bound == 'B')
scores <- predictions$pred
prg_curve <- prg::create_prg_curve(labels, scores)
auprg <- prg::calc_auprg(prg_curve)
message('AUPRG,', auprg, ',', predictions.path)
# convex_hull = prg::prg_convex_hull(prg_curve)
# fig = prg::plot_prg(prg_curve)
# print(fig)
# print(prg_curve)


# Use ROCR
pred <- ROCR::prediction(scores, labels)
perf <- ROCR::performance(pred, 'tpr', 'fpr')
