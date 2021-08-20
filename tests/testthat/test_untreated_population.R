source("../../R/untreated_population.R")

describe(" The function est_beta_params", {
  it("Return right answer", {
    expected <- list(alpha = NaN, beta = NaN)
    obtained <- est_beta_params(mu = 0, var = 1)
    expect_equal(expected, obtained)
  })
  it("Return right answer", {
    expected <- list(alpha = -1, beta = 0)
    obtained <- est_beta_params(mu = 1, var = 1)
    expect_equal(expected, obtained)
  })
  it("Return right answer", {
    expected <- list(alpha = -6.6, beta = 4.4)
    obtained <- est_beta_params(mu = 3, var = 5)
    expect_equal(expected, obtained)
  })
})

