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
  real_first_eigen_vector <- Re((eigen(x)$vectors)[, 1])
  parallel_matrix <- x %*% real_first_eigen_vector
  (parallel_matrix / (sum(parallel_matrix)))[, 1]
}
