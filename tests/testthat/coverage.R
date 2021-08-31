library(covr)
library(testthat)
cobertura <- file_coverage(
  c(
    "R/feral_cat.R",
    "R/monthly_matrix_leslie.R",
    "R/untreated_population.R"
  ),
  c(
    "tests/testthat/test_feral_cat.R",
    "tests/testthat/test_monthly_matrix_leslie.R",
    "tests/testthat/test_plots.R",
    "tests/testthat/test_untreated_population.R"
  )
)
covr::codecov(covertura = cobertura, token = "d40cba41-8ee3-414d-9e04-581d33a42b62")
