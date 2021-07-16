library(testthat)
library(FeralCatEradication)

test_that("all is ok", {
  expected_value = 1
  obtained_value = return_one()
  expect_equal(expected_value, obtained_value)
})

test_that("Eigenvalor", {
  expected_eigenvalue <- 3
  matriz <- matrix(c(2, 1, 1, 2), nrow = 2)
  obtained_eigenvalue <- max.lambda(matriz)
  expect_equal(expected_eigenvalue, obtained_eigenvalue)
})
