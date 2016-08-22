library(Saturn)
context("Rle construction")

test_that("both methods of Rle construction produce the same answers", {
  # For different lengths
  for (L in seq(1, 1000, by=5)) {
    # Sample a random number of random indexes
    .idxs <- sort(sample(L, sample(L, 1)))
    suppressWarnings(Rle.slow <- sparse.to.rle(L, .idxs, 1:length(.idxs)))
    # print(Rle.slow)
    Rle.fast <- Rle.from.sparse(L, .idxs, 1:length(.idxs))
    # print(Rle.fast)
    expect_equal(Rle.slow, Rle.fast)
  }
})
