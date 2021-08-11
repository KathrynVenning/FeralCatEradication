library(testthat)
library(FeralCatEradication)

describe("Get version of the module", {
  it("The version is 0.1.6", {
    expected_version <- c("0.1.6")
    obtained_version <- packageVersion("FeralCatEradication")
    version_are_equal <- expected_version == obtained_version
    expect_true(version_are_equal)
  })
})

maximum_age <- 7
fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
popmat <- matrix(data = 0, nrow = maximum_age, ncol = maximum_age)
diag(popmat[2:maximum_age, ]) <- survival_probability
popmat[maximum_age, maximum_age] <- 0
popmat[1, ] <- fertility
leslie_matrix_kathryn <- popmat

leslie_matrix_gotelli <- matrix(c(1.5, 1.5, 0.25, 0, 0.8, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.25, 0), nrow = 4, byrow = TRUE)

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
  it("Kathryn example. Maximum age: 7 years; Matrix: 7x7", {
    expected_eigenvalue <- 1.25
    obtained_eigenvalue <- max_lambda(leslie_matrix_kathryn)
    expect_equal(expected_eigenvalue, obtained_eigenvalue, tolerance = 1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_eigenvalue <- 2.095
    obtained_eigenvalue <- max_lambda(leslie_matrix_gotelli)
    expect_equal(expected_eigenvalue, obtained_eigenvalue, tolerance = 1e-3)
  })
})

describe("Get the Malthusian parameter from the Leslie Matrix", {
  it("Diagonal matrix 3x3", {
    expected_eigenvalue <- 3
    matriz <- matrix(c(exp(1), 0, 0, 0, exp(2), 0, 0, 0, exp(3)), nrow = 3)
    obtained_eigenvalue <- max_r(matriz)
    expect_equal(expected_eigenvalue, obtained_eigenvalue)
  })
  it("Kathryn example. Maximum age: 7 years; Matrix: 7x7", {
    expected_eigenvalue <- 0.223
    obtained_eigenvalue <- max_r(leslie_matrix_kathryn)
    expect_equal(expected_eigenvalue, obtained_eigenvalue, tolerance = 1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_eigenvalue <- 0.7397
    obtained_eigenvalue <- max_r(leslie_matrix_gotelli)
    expect_equal(expected_eigenvalue, obtained_eigenvalue, tolerance = 1e-3)
  })
})

describe("Get stable stage distribution", {
  it("Example 6.2.1 of Grossman", {
    expected_stable_stage_distribution <- c(0.22, 0.78)
    matriz <- matrix(c(0, 2, 0.3, 0.5), nrow = 2)
    obtained_stable_stage_distribution <- stable_stage_dist(matriz)
    expect_equal(expected_stable_stage_distribution, obtained_stable_stage_distribution, tolerance = 1e-3)
  })
  it("Kathryn example. Maximum age: 7 years; Matrix: 7x7", {
    expected_stable_stage_distribution <- c(0.6026, 0.2219, 0.0817, 0.0458, 0.0256, 0.0144, 0.00733)
    obtained_stable_stage_distribution <- stable_stage_dist(leslie_matrix_kathryn)
    expect_equal(expected_stable_stage_distribution, obtained_stable_stage_distribution, tolerance = 1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_stable_stage_distribution <- c(0.67397, 0.25731, 0.06140, 0.00733)
    obtained_stable_stage_distribution <- stable_stage_dist(leslie_matrix_gotelli)
    expect_equal(expected_stable_stage_distribution, obtained_stable_stage_distribution, tolerance = 1e-3)
  })
})

describe("total_female_offspring_per_female", {
  it("Maximum age: 3 years; matrix: 3x3", {
    expected_total_female_offspring_per_female <- c(1)
    maximum_age <- 3
    leslie_matrix <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_total_female_offspring_per_female <- total_female_offspring_per_female(leslie_matrix, maximum_age)
    expect_equal(expected_total_female_offspring_per_female, obtained_total_female_offspring_per_female, tolerance = 1e-3)
  })
  it("Kathryn example. Maximum age: 7 years; Matrix: 7x7", {
    expected_total_female_offspring_per_female <- c(2.0427)
    maximum_age <- 7
    obtained_total_female_offspring_per_female <- total_female_offspring_per_female(leslie_matrix_kathryn, maximum_age)
    expect_equal(expected_total_female_offspring_per_female, obtained_total_female_offspring_per_female, tolerance = 1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_total_female_offspring_per_female <- c(2.8)
    maximum_age <- 4
    obtained_total_female_offspring_per_female <- total_female_offspring_per_female(leslie_matrix_gotelli, maximum_age)
    expect_equal(expected_total_female_offspring_per_female, obtained_total_female_offspring_per_female, tolerance = 1e-3)
  })
})

describe("Mean generation time function", {
  it("Maximum age: 3 years; Diagonal matrix: 3x3", {
    expected_mean_generation <- c(0)
    maximum_age <- 3
    leslie_matrix <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
    obtained_mean_generation <- g_val(leslie_matrix, maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance = 1e-3)
  })
  it("Kathryn example. Maximum age: 7 years; matrix: 7x7", {
    expected_mean_generation <- c(3.21)
    maximum_age <- 7
    obtained_mean_generation <- g_val(leslie_matrix_kathryn, maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance = 1e-3)
  })
  it("Gotelli example. Maximum age: 4 years; Diagonal matrix: 4x4", {
    expected_mean_generation <- c(1.392)
    maximum_age <- 4
    obtained_mean_generation <- g_val(leslie_matrix_gotelli, maximum_age)
    expect_equal(expected_mean_generation, obtained_mean_generation, tolerance = 1e-3)
  })
})

describe("Get parameters of survival modifier", {
  pop_found <- 1629
  k_max <- 2 * pop_found
  k_vec <- c(1, pop_found / 2, pop_found, 0.75 * k_max, k_max)
  red_vec <- c(1, 0.965, 0.89, 0.79, 0.71)
  obtained_coefficients <- coefficients_proportion_realized_survival(k_vec, red_vec)
  assert_survival_modifier <- function(i) {
    expect_equal(red_vec[i], survival_modifier(k_vec[i], obtained_coefficients), tolerance = 1e-2)
  }
  it("Kathryn example", {
    a_lp <- 1.001
    b_lp <- 5459.994
    c_lp <- 1.690
    expected_coefficients <- list(a_lp = a_lp, b_lp = b_lp, c_lp = c_lp)
    expect_equal(expected_coefficients, obtained_coefficients, tolerance = 1e-3)
    for (i in 1:length(k_vec)) {
      assert_survival_modifier(i)
    }
  })
})
