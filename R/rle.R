# Rle Tips and Tricks
# Patrick Aboyoun
# July 10, 2016

rollmeanRle <- function (x, k)
{
  n <- length(x)
  cumsum(c(Rle(sum(window(x, 1, k))), window(x, k + 1, n) - window(x, 1, n - k))) / k
}

rollvarRle <- function(x, k)
{
  n <- length(x)
  means <- rollmeanRle(x, k)
  nextMean <- window(means, 2, n - k + 1)
  cumsum(c(Rle(sum((window(x, 1, k) - means[1])^2)),
  k * diff(means)^2 -
  (window(x, 1, n - k) - nextMean)^2 +
  (window(x, k + 1, n) - nextMean)^2)) / (k - 1)
}

rollcovRle <- function(x, y, k)
{
  n <- length(x)
  meanX <- rollmeanRle(x, k)
  meanY <- rollmeanRle(y, k)
  nextMeanX <- window(meanX, 2, n - k + 1)
  nextMeanY <- window(meanY, 2, n - k + 1)
  cumsum(c(Rle(sum((window(x, 1, k) - meanX[1]) * (window(y, 1, k) - meanY[1]))),
  k * diff(meanX) * diff(meanY) -
  (window(x, 1, n - k) - nextMeanX) * (window(y, 1, n - k) - nextMeanY) +
  (window(x, k + 1, n) - nextMeanX) * (window(y, k + 1, n) - nextMeanY))) / (k - 1)
}

rollsdRle <- function(x, k)
{
  sqrt(rollvarRle(x, k))
}

rollcorRle <- function(x, y, k)
{
  rollcovRle(x, y, k) / (rollsdRle(x, k) * rollsdRle(y, k))
}
