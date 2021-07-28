#' @export
max_lambda <- function(x) {
  Re((eigen(x)$values)[1])
}

#' @export
max_r <- function(x) {
  log(max_lambda(x))
}

#' @export
stable_stage_dist <- function(x) {
  ((x %*% (Re((eigen(x)$vectors)[, 1]))) / (sum((x %*% (Re((eigen(x)$vectors)[, 1]))))))[, 1]
}
