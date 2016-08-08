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