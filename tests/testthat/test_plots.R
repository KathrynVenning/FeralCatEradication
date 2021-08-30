library(testthat)
library(vdiffr)
library(ggplot2)
library(FeralCatEradication)

fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
initial_population <- 1629

# Project
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date
survival <- Survival_Fertility$new(fertility, survival_probability)
population_with_cc <- Population$new(survival)
population_with_cc$run_generations(yr_now, yr_end, initial_population = initial_population)


test_that("ggplot2 histogram works", {
  plotter <- Plotter_Population$new()
  p <- plotter$plot(population_with_cc)
  expect_doppelganger("default histogram", p)
})

test_that("base histogram works", {
  p <- function() hist(mtcars$disp)
  expect_doppelganger("base histogram", p)
})