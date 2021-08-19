# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# libraries
library(plotly)
library(FeralCatEradication)
options(scipen = 1000)
library(tidyverse)

# create Leslie matrix

# Create vectors
# Fertility
# KI cat birth rates matrix, data for female offsping produced each year.
# Data from Budke, C & Slater, M (2009)
fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
age_max <- length(fertility)

# Fertility errors based on Budke & Slater
# Mean and standard deviations, juvenile fertility:
# Survival
# KI cat survival
# probability of surviving from one year to the next. e.g surviving fourth year of life
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
# create matrix
popmat <- FeralCatEradication::matrix_leslie(fertility, survival_probability)

# matrix properties
gen_l <- FeralCatEradication::g_val(popmat, age_max) # mean generation length

# initial population vector
pop_found <- 1629 # +/- 661 founding population size Hohnen et al 2020
ssd <- FeralCatEradication::stable_stage_dist(popmat)
init_vec <- ssd * pop_found # initial population vector

# Project
# set time limit for projection in 1-yr increments
yr_now <- 2020 # update if more data available post-2010
#************************
yr_end <- 2030 # set projection end date
#************************

# set population storage matrices
population_with_cc <- Population$new(fertility, survival_probability)
population_with_cc$run_generations(yr_now, yr_end, initial_population = init_vec)

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
population_with_cc$run_generations(yr_now, yr_end, initial_population = init_vec, coefficients)

plotter$plot(population_with_cc)
plotter$save("reports/figures/time_serie_individuals_with_carry_capacity.jpg")
