# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# libraries
library(plotly)
library(FeralCatEradication)
options(scipen = 1000)
library(tidyverse)

fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
initial_population <- 1629

# Project
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date

population_with_cc <- Population$new(fertility, survival_probability)
population_with_cc$run_generations(yr_now, yr_end, initial_population = initial_population)

plotter <- Plotter_Population$new()
plotter$plot(population_with_cc)
plotter$save("reports/figures/time_serie_individuals.jpg")

# Compensatory density feedback
# K = carry capacity
# Population rate of increase relative to carry capacity
# Larger distance between populationa and K = faster population growth
k_max <- 2 * pop_found
k_vec <- c(1, pop_found / 2, pop_found, 0.75 * k_max, k_max)
red_vec <- c(1, 0.965, 0.89, 0.79, 0.71)
coefficients <- coefficients_proportion_realized_survival(k_vec, red_vec)

# compensatory density-feedback deterministic model
# set population storage matrices
population_with_cc$run_generations(yr_now, yr_end, initial_population = initial_population, coefficients)

plotter$plot(population_with_cc)
plotter$save("reports/figures/time_serie_individuals_with_carry_capacity.jpg")
