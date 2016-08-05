#!/usr/bin/env Rscript
#
# Examine DNase levels against binding
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


#' Analyse each TF
#'
analyse.tf <- function(tf.meta) {
  tf <- tf.meta$TF[1]
  cell.all <- tf.meta$cell
  #
  # Load the binding data
  #
  chip <- load.chip.features(tf)
  #
  # Load DNase features
  #
  dnase <- load.dnase.features(cell.all)
  #'
  #' Analyse the DNase levels in the cell aggregated by binding status
  #'
  cell.dnase.by.bound <- function(cell.row) {
    cell <- as.character(cell.row$cell)
    message(tf, ' : ', cell)
    cell.bound <- chip[[cell]]
    cell.dnase <- dnase[[cell]][training.region.test.idxs()]
    zero.dnase <- 0 == cell.dnase
    is.bound <- 'B' == cell.bound
    is.ambig <- 'A' == cell.bound
    is.unbnd <- 'U' == cell.bound
    by.bound <- function(FUN) c(
      FUN(cell.dnase[is.bound]),
      FUN(cell.dnase[is.ambig]),
      FUN(cell.dnase[is.unbnd]))
    data.frame(
      cell = factor(cell, cell.levels),
      bound = factor(c('B', 'A', 'U'), levels = binding.levels),
      prop.zero.dnase = c(
        sum(is.bound & zero.dnase) / sum(is.bound),
        sum(is.ambig & zero.dnase) / sum(is.ambig),
        sum(is.unbnd & zero.dnase) / sum(is.unbnd)),
      mean = by.bound(mean),
      sd = by.bound(sd))
  }
  tf.meta %>% rowwise() %>% do(cell.dnase.by.bound(.))
}
dnase.by.bound <- tfs %>% filter('train' == split) %>% group_by(TF) %>% do(analyse.tf(.)) %>% ungroup()
dnase.by.bound
# devtools::use_data(dnase.by.bound)

ggplot(dnase.by.bound, aes(x = bound, y = prop.zero.dnase, color = bound, fill = bound)) +
  # geom_boxplot() +
  geom_jitter(height = 0) +
  facet_wrap(~ TF, nrow = 4) +
  scale_fill_few() +
  theme_few()
ggsave(file.path('..', 'Plots', 'dnase-by-bound.pdf'))
