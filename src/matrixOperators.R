# Matrix operators for population models
library(MASS)

# Mean generation time function
# where leslie_matrix is a Leslie Matrix
g_val <- function(leslie_matrix, age_max) {
  mean_generation_time <- (log(r_val(leslie_matrix, age_max))) / (log(Re((eigen(leslie_matrix)$values)[1])))
  print("mean generation time")
  print("____________________")
  print(mean_generation_time)
}
