library(Saturn)
context("Rle construction")

test_that("both methods of Rle construction produce the same answers", {
  for (L in seq(10, 10, step=10)) {
    .idxs <- sample(L, sample(L, 1))
    Rle.slow <- sparse.to.rle(L, .idxs, 1:length(.idxs))
    Rle.fast <- Rle.from.sparse(L, .idxs, 1:length(.idxs))
    message(Rle.slow)
    message(Rle.fast)
    expect_equal(Rle.slow, Rle.fast)
  }
})
