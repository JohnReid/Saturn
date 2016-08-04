#!/usr/bin/env Rscript
#
# Examine correlations between TF binding
#


#
# Set warnings as errors
#
options(warn = 2)


#
# Load libraries
devtools::load_all()
library(Saturn)
library(ggplot2)
library(ggthemes)


#' Convert a factor Rle to an integer Rle
Rle.factor.as.integer <- function(x) Rle(as.integer(runValue(x)), runLength(x))

x <- as.integer(chip[['GM12878']])
y <- as.integer(chip[['HeLa-S3']])
x <- chip[['GM12878']]
y <- chip[['HeLa-S3']]
x.rle <- Rle(as.integer(runValue(x)), runLength(x))
y.rle <- Rle(as.integer(runValue(y)), runLength(y))
pearson.cor <- function(x, y) {
  mean.x <- sum(x) / length(x)
  mean.y <- sum(y) / length(y)
  x. <- x - mean.x
  y. <- y - mean.y
  sum(x. * y.) / sqrt(sum(x.^2) * sum(y.^2))
}
system.time(pearson.cor(x.rle, y.rle))
cor(x, y)

#' Analyse each TF
#'
analyse.tf <- function(tf.meta) {
  tf <- tf.meta$TF[1]
  cell.all <- tf.meta$cell
  #
  # Check we have at least two cell types
  #
  if (2 > length(cell.all)) {
    return(data.frame(cell.1 = factor(levels = cell.levels), cell.2 = factor(levels = cell.levels), cor = integer()))
  }
  #
  # Load the binding data
  #
  chip <- load.chip.features(tf)
  #
  # Generate all pairs of cells
  #
  cell.pairs <- combn(cell.all, 2)
  #
  # Calculate correlation of each pair
  #
  cors <- sapply(1:ncol(cell.pairs), function(i) {
    cell.1 <- cell.pairs[1,i]
    cell.2 <- cell.pairs[2,i]
    pearson.cor(
      Rle.factor.as.integer(chip[[as.character(cell.1)]]),
      Rle.factor.as.integer(chip[[as.character(cell.2)]]))
  })
  data.frame(cell.1 = cell.pairs[1,], cell.2 = cell.pairs[2,], cor = cors)
}
cor.by.tf <- tfs %>% filter('train' == split) %>% group_by(TF) %>% do(analyse.tf(.)) %>% ungroup()
sapply(cor.by.tf, class)
cor.by.tf
devtools::use_data(cor.by.tf)


#
# Make a plot
#
ggplot(cor.by.tf, aes(x = cell.1, y = cell.2, fill = cor)) +
  # geom_boxplot() +
  geom_tile() +
  facet_wrap(~ TF, nrow = 5) +
  # scale_fill_few() +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path('..', 'Plots', 'cor-by-tf.pdf'))
