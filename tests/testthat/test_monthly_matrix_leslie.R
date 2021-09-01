#library(FeralCatEradication)
setwd("/workspaces/FeralCatEradication")
source("R/feral_cat.R")
source("R/monthly_matrix_leslie.R")
source("R/untreated_population.R")

describe("monthly_matrix_leslie", {
  fertility <- c((0.745 / 3), 2.52, 1.98)
  survival_probability <- c(0.46, 0.7)
  it(" Fertility for tree classes of age", {
    expected <- c(rep((0.745 / 3) / 12, 12), rep(2.52 / 12, 12), rep(1.98 / 12, 12))
    obtained <- monthly_matrix_leslie(fertility, survival_probability)
    expect_equal(expected, obtained[1, ])
  })
  it("Survival for tree classes of age", {
    expected <- c(rep((0.46)^(1 / 12), 12), rep(0.7^(1 / 12), 23))
    obtained <- monthly_matrix_leslie(fertility, survival_probability)
    expect_equal(expected, diag(obtained[2:36, ]))
  })
  fertility <- c((0.745 / 3), 2.52, 2.52, 1.98)
  survival_probability <- c(0.46, 0.7, 0.7)
  it(" Fertility for four classes of age", {
    expected <- c(rep((0.745 / 3) / 12, 12), rep(2.52 / 12, 24), rep(1.98 / 12, 12))
    obtained <- monthly_matrix_leslie(fertility, survival_probability)
    expect_equal(expected, obtained[1, ])
  })
})
