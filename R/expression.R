#' Load expression data.
#'
saturn.expr <- memoise::memoise(function(cell, biorep) {
  path <- file.path(
    saturn.data(), 'RNAseq',
    sprintf('gene_expression.%s.biorep%d.tsv', cell, biorep))
  message('Loading: ', path)
  readr::read_tsv(
    path, skip = 1, progress = FALSE,
    col_names = c(
      'gene_id',
      'transcript_ids',
      'length',
      'effective_length',
      'expected_count',
      'TPM',
      'FPKM'),
    col_types = readr::cols(
      gene_id = readr::col_character(),
      transcript_ids = readr::col_character(),
      .default = readr::col_double()))
})


#' Load average expression across both replicates
#'
load.avg.expr <- memoise::memoise(function(cell) {
    rep.1 <- saturn.expr(cell, 1) %>% mutate(gene_id = factor(gene_id))
    rep.2 <- saturn.expr(cell, 2) %>% mutate(gene_id = factor(gene_id))
    rbind(rep.1, rep.2) %>%
      group_by(gene_id) %>%
      summarise(log.TPM = mean(log10(TPM+1)), log.FPKM = mean(log10(FPKM+1))) %>%
      mutate(ensembl_gene_id = factor(stringr::str_split_fixed(gene_id, stringr::fixed('.'), 2)[,1]))
})


#' Get the expression features for the cell
#'
expr.features.for.cell <- function(.cell) {
  cluster.expr <- filter(expr.cell.cluster, .cell == cell)
  cell.expr <- cluster.expr$log.TPM
  names(cell.expr) <- stringr::str_c('expr.', cluster.expr$cluster_id)
  #
  # Add the genes of interest expression levels
  #
  c(cell.expr, expr.of.interest[,.cell])
}


#' The directory for DNase features
#'
expr.features.dir <- function() file.path(saturn.data(), 'Features', 'Expr')


#' The file for the nearest gene features
#'
nearest.gene.feature.file <- function()
  file.path(expr.features.dir(), 'nearest-gene.rds')


#' Load the expression features
#'
load.expr.feature <- memoise::memoise(function() {
  message('Loading nearest gene features')
  readRDS(nearest.gene.feature.file())
})


#' The file for the nearest gene expression feature
#'
nearest.expr.feature.file <- function(cell) file.path(
  expr.features.dir(), stringr::str_c('nearest-expr-', cell, '.rds'))


#' Load the DNase features for the cell
#'
load.nearest.expr.feature <- memoise::memoise(function(cell) {
  message('Loading nearest expression features for cell ', cell)
  readRDS(nearest.expr.feature.file(cell))
})


