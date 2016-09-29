
smoothing.params <- rbind(
  data.frame(
    TF      = c("ARID3A", "ATF2",   "ATF3",   "ATF7",   "CEBPB",  "CREB1",  "CTCF"),
    smooth  = c( TRUE,     TRUE,     TRUE,     TRUE,     FALSE,    TRUE,     FALSE),
    L       = c( 200,      200,      200,      200,      10,       50,       10),
    W       = 10,
    LO      = TRUE
  ),
  data.frame(
    TF      = c("E2F1",   "E2F6",   "EGR1",   "EP300",  "FOXA1",  "FOXA2",  "GABPA"),
    smooth  = c( TRUE,     TRUE,     TRUE,     FALSE,    TRUE,     TRUE,     TRUE),
    L       = c( 200,      200,      50,       10,       200,      200,      50),
    W       = 10,
    LO      = TRUE
  ),
  data.frame(
    TF      = c("GATA3",  "HNF4A",  "JUND",   "MAFK",   "MAX",    "MYC",    "NANOG"),
    smooth  = c( FALSE,    TRUE,     FALSE,    FALSE,    TRUE,     TRUE,     TRUE),
    L       = c( 200,      50,       10,       10,       50,       50,       50),
    W       = 10,
    LO      = TRUE
  ),
  data.frame(
    TF      = c("REST",   "RFX5",   "SPI1",   "SRF",    "STAT3",  "TAF1",   "TCF12"),
    smooth  = c( TRUE,     TRUE,     TRUE,     TRUE,     TRUE,     TRUE,     TRUE),
    L       = c( 10,       200,      200,      200,      200,      50,       200),
    W       = 10,
    LO      = TRUE
  ),
  data.frame(
    TF      = c("TCF7L2", "TEAD4",  "YY1",    "ZNF143"),
    smooth  = c( TRUE,     TRUE,     TRUE,     TRUE),
    L       = c( 200,      10,       50,       50),
    W       = 10,
    LO      = TRUE
  ))


#
# Logistic transform and its inverse logit
#
logistic <- function(x) 1/(1+exp(-x))
logit <- function(p) log(p/(1-p))

#' Create a function that calculates the similarity between two regions in the predictions indexed by row
#'
region.similarity <- function(preds, length.scale) {
  function(i, j) {
    ifelse(
      preds$chrom[i] == preds$chrom[j],
      exp(-((preds$start[i] - preds$start[j])/length.scale)**2/2),
      0)
  }
}


#' Smooth the predictions
#'
smooth.predictions <- function(preds, length.scale, max.width, log.transform = TRUE) {
  #
  # Distance between two rows
  #
  row.similarity <- region.similarity(preds, length.scale)

  #
  # Create sparse symmetric banded similarity matrix
  #
  # Calculate non-zero indices and values
  N <- nrow(preds)
  i <- as.vector(sapply(1:N, function(i) rep(i, max.width)))
  j <- rep(1:max.width, N) + i - 1
  x <- row.similarity(i, j)
  #
  # Discard indices outside matrix
  keep <- j <= N
  sum(! keep)
  i <- i[keep]
  j <- j[keep]
  x <- x[keep]
  #
  # Create matrix
  K <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(N, N), symmetric = TRUE)
  # object.size(K) / dim(K)[1] / dim(K)[2]
  # class(K)

  #
  # Normalise smoothing matrix by rows
  #
  row.sums <- Matrix::rowSums(K)
  K <- Matrix::Diagonal(x = 1 / row.sums) %*% K
  # stopifnot(all(abs(1 - Matrix::rowSums(K)) < 1e-12))
  # object.size(K) / dim(K)[1] / dim(K)[2]
  # class(K)

  #
  # Smooth predictions
  #
  predictions <- preds$prediction
  if (log.transform) {
    predictions <- logit(predictions)
  }
  predictions.smoothed <- as.vector(K %*% predictions)
  if (log.transform) {
    predictions.smoothed <- logistic(predictions.smoothed)
  }
  class(predictions.smoothed)

  #
  # Return result
  #
  stopifnot (! is.null(predictions.smoothed))
  return(predictions.smoothed)
}
