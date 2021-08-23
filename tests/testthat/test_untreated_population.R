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

describe("get_stochastic_fertility", {
  it("limits are correts", {
    fertility <- c(3, 4)
    sd_fertility <- c(0.1, 0.2)
    obtained_stochastic_fertility <- get_stochastic_fertility(fertility, sd_fertility)
    all_less_than_five <- all(obtained_stochastic_fertility < 5)
    expect_true(all_less_than_five)
    all_greater_than_two <- all(obtained_stochastic_fertility > 2)
    expect_true(all_greater_than_two)
  })
})

describe("The class Survival_Fertility", {
  it("The method builder exist", {
    fertility <- seq(1,4)
    survival_probability <- rbeta(3,1,1)
    survival <- Survival_Fertility$new(fertility, survival_probability)
  })
})
