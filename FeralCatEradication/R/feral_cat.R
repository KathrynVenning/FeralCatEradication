#' @export
max.lambda <- function(x) {
  Re((eigen(x)$values)[1])
}

max.r <- function(x) {
  log(max.lambda(x))
}

stable_stage_dist <- function(x) {
  ((x %*% (Re((eigen(x)$vectors)[, 1]))) / (sum((x %*% (Re((eigen(x)$vectors)[, 1]))))))[, 1]
}
