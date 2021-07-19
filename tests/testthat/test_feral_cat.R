library(testthat)
library(FeralCatEradication)

test_that("all is ok", {
  expected_value = 1
  obtained_value = return_one()
  expect_equal(expected_value, obtained_value)
})

describe("Get the first eigenvalue og the Leslie Matrix", {
  it("Matrix 2x2 real eigenvalues", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(2, 1, 1, 2), nrow = 2)
    obtained_eigenvalue <- max.lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
  it("Matrix 3x3 with complex eigenvalues", {
    expected_eigenvalue <- -0.5
    matriz <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3)
    obtained_eigenvalue <- max.lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
  it("Diagonal matrix 3x3", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_eigenvalue <- max.lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
})

describe("Get the Malthusian parameter from the Leslie Matrix", {
  it("Diagonal matrix 3x3", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(exp(1), 0, 0, 0, exp(2), 0, 0, 0, exp(3)), nrow = 3)
    obtained_eigenvalue <- max.r(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
})