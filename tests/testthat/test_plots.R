library(testthat)
library(vdiffr)
library(ggplot2)
library(FeralCatEradication)
source("../../R/feral_cat.R")

fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
initial_population <- 1629

# Project
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date
interval_time <- Interval_Time$new(initial_year = yr_now, final_year = yr_end)
survival <- Survival_Fertility$new(fertility, survival_probability)
population_with_cc <- Population$new(survival)
population_with_cc$run_generations(interval_time, initial_population = initial_population)


test_that("Population anualy from 2020 to 2030", {
  plotter <- Plotter_Population$new()
  p <- plotter$plot(population_with_cc)
  expect_doppelganger("default histogram", p)
})

test_that("base histogram works", {
  p <- function() hist(mtcars$disp)
  expect_doppelganger("base histogram", p)
})
