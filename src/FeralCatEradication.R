# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# libraries
library(plotly)
library(FeralCatEradication)
source("R/feral_cat.R")
options(scipen = 1000)
library(tidyverse)

fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
initial_population <- 1629

# Project
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date
interval_time <- Interval_Time$new(initial_year = yr_now, final_year = yr_end)
survival <- Survival_Fertility$new(fertility, survival_probability)
population_with_cc <- Population$new(survival)
simulator <- Runner_Population$new(population_with_cc)
simulator$run_generations(interval_time, initial_population = initial_population)

plotter <- Plotter_Population$new()
plotter$plot(simulator)
plotter$save("reports/figures/time_serie_individuals.jpg")

# Compensatory density feedback
# K = carry capacity
# Population rate of increase relative to carry capacity
# Larger distance between populationa and K = faster population growth
capacity <- Carry_Capacity$new()
coefficients <- capacity$coefficients_model(half_capacity = initial_population)
# compensatory density-feedback deterministic model
# set population storage matrices
simulator2 <- Runner_Population_With_CC$new(population_with_cc, coefficients)
simulator2$run_generations(interval_time, initial_population = initial_population)

plotter$plot(simulator2)
plotter$plot_carry_capacity(capacity)
plotter$save("reports/figures/time_serie_individuals_with_carry_capacity.jpg")
