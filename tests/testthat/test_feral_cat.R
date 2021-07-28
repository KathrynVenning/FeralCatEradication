library(testthat)
library(FeralCatEradication)

describe("Get the first eigenvalue og the Leslie Matrix", {
  it("Matrix 2x2 real eigenvalues", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(2, 1, 1, 2), nrow = 2)
    obtained_eigenvalue <- max_lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
  it("Matrix 3x3 with complex eigenvalues", {
    expected_eigenvalue <- -0.5
    matriz <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3)
    obtained_eigenvalue <- max_lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
  it("Diagonal matrix 3x3", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_eigenvalue <- max_lambda(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
})

describe("Get the Malthusian parameter from the Leslie Matrix", {
  it("Diagonal matrix 3x3", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(exp(1), 0, 0, 0, exp(2), 0, 0, 0, exp(3)), nrow = 3)
    obtained_eigenvalue <- max_r(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
})

describe("Get stable stage distribution", {
  it("Example 6.2.1 of Grossman", {
    expected_stable_stage_distribution <- c(0.22, 0.78)
    matriz <- matrix(c(0, 2, 0.3, 0.5), nrow = 2)
    obtained_stable_stage_distribution <- stable_stage_dist(matriz)
    expect_equal(expected_stable_stage_distribution, obtained_stable_stage_distribution, tolerance=1e-3)
  })
})
