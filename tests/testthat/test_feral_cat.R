library(testthat)
library(FeralCatEradication)

describe("Get version of the module", {
  it("The version is 0.1.3", {
    expected_version <- c("0.1.3")
    obtained_version <- packageVersion("FeralCatEradication")
    version_are_equal <- expected_version == obtained_version
    expect_true(version_are_equal)
  })
})

describe("Get the first eigenvalue of the Leslie Matrix", {
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

describe("total_female_offspring_per_female", {
  it("Maximum age: 3 years; matrix: 3x3", {
    expected_total_female_offspring_per_female <- c(1)
    maximum_age <- 3
    leslie_matrix <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_total_female_offspring_per_female <- total_female_offspring_per_female(leslie_matrix,maximum_age)
    expect_equal(expected_total_female_offspring_per_female, obtained_total_female_offspring_per_female, tolerance=1e-3)
  })
})

describe("Mean generation time function", {
  it("Maximum age: 3 years; Diagonal matrix: 3x3", {
    expected_mean_generation <- c(0)
    maximum_age <- 3
    leslie_matrix <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_mean_generation <- g_val(leslie_matrix,maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance=1e-3)
  })
  it("Kathryn example. Maximum age: 7 years; matrix: 7x7", {
    expected_mean_generation <- c(3.21)
    maximum_age <- 7
    m_vec <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
    s_vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
    popmat <- matrix(data = 0, nrow = maximum_age, ncol = maximum_age)
    diag(popmat[2:maximum_age, ]) <- s_vec
    popmat[maximum_age, maximum_age] <- 0
    popmat[1, ] <- m_vec
    leslie_matrix <- popmat
    obtained_mean_generation <- g_val(leslie_matrix,maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance=1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_mean_generation <- c(0)
    maximum_age <- 4
    leslie_matrix <- matrix(c(1.5, 1.5, 0.25, 0, 0.8, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.25, 0), nrow = 4, byrow=TRUE)
    obtained_mean_generation <- g_val(leslie_matrix,maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance=1e-3)
  })
})
