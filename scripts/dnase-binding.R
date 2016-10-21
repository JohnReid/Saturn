#!/usr/bin/env Rscript
#
# Examine DNase levels against binding
#


#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
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
      total = c(sum(is.bound), sum(is.ambig), sum(is.unbnd)),
      mean = by.bound(mean),
      sd = by.bound(sd))
  }
  tf.meta %>% rowwise() %>% do(cell.dnase.by.bound(.))
}
dnase.by.bound <- tfs %>% filter('train' == split) %>% group_by(TF) %>% do(analyse.tf(.)) %>% ungroup()
dnase.by.bound
# devtools::use_data(dnase.by.bound, overwrite = TRUE)

#
# Plot how the proportion of zero DNase regions varies as function of bound status by TF
#
ggplot(dnase.by.bound, aes(x = bound, y = prop.zero.dnase, color = bound, fill = bound)) +
  # geom_boxplot() +
  geom_jitter(height = 0, size = 2) +
  facet_wrap(~ TF, nrow = 4) +
  scale_fill_few() +
  theme_few()
ggsave(file.path('..', 'Plots', 'dnase-by-bound.pdf'))

#
# Plot the proportion bound by TF
#
ggplot(dnase.by.bound, aes(x = bound, y = total, color = bound, fill = bound)) +
  # geom_boxplot() +
  # geom_jitter(height = 0, size = 2) +
  geom_bar(stat = "identity") +
  facet_grid(TF ~ cell) +
  scale_fill_few() +
  scale_y_log10() +
  theme_few()
  # theme(
    # axis.text.x = element_text(angle = 90, hjust = 1),
    # axis.text.y = element_text(angle = 90, hjust = 1))
ggsave(file.path('..', 'Plots', 'bound-totals.pdf'))

#
# Calculate the total number bound, ambiguous and
#
totals <-
  dnase.by.bound %>%
  dcast(TF + cell ~ bound, value.var = 'total') %>%
  mutate(
    bound.to.unbnd = B / U,
    bound.to.ambig = B / A)
totals
ggplot(totals, aes(x = TF, y = cell, fill = log10(bound.to.unbnd))) +
    geom_tile(colour='white') +
    theme_tufte() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4))
ggsave(
  file.path('..', 'Plots', 'bound-to-unbnd.pdf'),
  width = 6, height = 4, units = 'in')
ggplot(totals, aes(x = TF, y = cell, fill = log10(bound.to.ambig))) +
    geom_tile(colour='white') +
    theme_tufte() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4))
ggsave(
  file.path('..', 'Plots', 'bound-to-ambig.pdf'),
  width = 6, height = 4, units = 'in')
