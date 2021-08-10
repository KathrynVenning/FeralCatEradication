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

# Generation length function
# reproductive value (r_0) where leslie_matrix = Leslie matrix; age_max = maximum age of females
#' @export
total_female_offspring_per_female <- function(leslie_matrix, age_max) {
  # define the transition matrix
  transition_matrix <- leslie_matrix[1:age_max, 1:age_max]
  transition_matrix[1, 1:(age_max)] <- 0

  # define the fertility matrix
  fertility_matrix <- leslie_matrix[1:age_max, 1:age_max]
  diag(fertility_matrix[2:age_max, 1:(age_max - 1)]) <- 0

  # define the identity matrix
  identity_matrix <- matrix(data = 0, nrow = age_max, ncol = age_max)
  diag(identity_matrix) <- 1

  # define the fundamental matrix
  n_fund <- MASS::ginv(identity_matrix - transition_matrix)

  # define the reproductive matrix
  reproductive_matrix <- fertility_matrix %*% n_fund

  # define r_0 (number of female offspring produced per female during lifetime)
  r_0 <- Re((eigen(reproductive_matrix)$values)[1])

  return(r_0)
}

# Mean generation time function
# where leslie_matrix is a Leslie Matrix
#' @export
g_val <- function(leslie_matrix, age_max) {
  mean_generation_time <- (log(total_female_offspring_per_female(leslie_matrix, age_max))) / (log(Re((eigen(leslie_matrix)$values)[1])))
  return(mean_generation_time)
}

#' @export
coefficients_proportion_realized_survival <- function(k_vec, red_vec){
  k_red_dat <- data.frame(k_vec, red_vec)
  param_init <- c(1, 15000, 2.5)
  fit_lp <- nls(red_vec ~ a / (1 + (k_vec / b)^c),
    data = k_red_dat,
    algorithm = "port",
    start = c(a = param_init[1], b = param_init[2], c = param_init[3]),
    trace = TRUE,
    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1 / 1024)
  )
  a_lp <- coef(fit_lp)[1]
  b_lp <- coef(fit_lp)[2]
  c_lp <- coef(fit_lp)[3]
  coefficients <- list(a_lp, b_lp, c_lp)
}