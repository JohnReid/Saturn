#' The names of the memoised functions in the Saturn package
#'
memoised.fn.names <- function(environment.name = 'package:Saturn') {
  obj.names <- ls(environment.name)
  fn.names <- obj.names[sapply(lapply(obj.names, get), is.function)]
  fn.names[sapply(lapply(fn.names, get), memoise::is.memoised)]
}

#' Clear the caches of all the memoised functions in the package
#'
forget.all.memoised <- function(environment.name = 'package:Saturn') {
  fn.names <- memoised.fn.names(environment.name)
  result <- lapply(lapply(fn.names, get), memoise::forget)
  names(result) <- fn.names
  invisible(result)
}
