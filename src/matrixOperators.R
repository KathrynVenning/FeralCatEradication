# Matrix operators for population models
library(MASS)

# Stable stage distribution
stable_stage_dist <- function(x) {
  ((x %*% (Re((eigen(x)$vectors)[, 1]))) / (sum((x %*% (Re((eigen(x)$vectors)[, 1]))))))[, 1]
}

# Generation length function
# reproductive value (r_0) where leslie_matrix = Leslie matrix; age_max = maximum age of females
r_val <- function(leslie_matrix, age_max) {
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

  # output
  print("number of female offspring produced per female during its lifetime")
  print("_________________________________________________________________")
  print(r_0)

}

# Mean generation time function
# where leslie_matrix is a Leslie Matrix
g_val <- function(leslie_matrix, age_max) {
  mean_generation_time <- (log(r_val(leslie_matrix, age_max))) / (log(Re((eigen(leslie_matrix)$values)[1])))
  print("mean generation time")
  print("____________________")
  print(mean_generation_time)
}
