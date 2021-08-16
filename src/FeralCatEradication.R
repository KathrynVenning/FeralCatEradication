# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# remove everything
rm(list = ls())

# libraries
library(plotly)
library(FeralCatEradication)
options(scipen = 1000)
library(tidyverse)

# functions
# beta distribution shape parameter estimator function
est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# create Leslie matrix
age_max <- 7

# Create vectors
# Fertility
# KI cat birth rates matrix, data for female offsping produced each year.
# Data from Budke, C & Slater, M (2009)
fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)

# Fertility errors based on Budke & Slater
# Mean and standard deviations, juvenile fertility:
juv_m_sd <- mean(c(((0.745 / 3 - 0.352 / 3) / 2), ((1.58 / 3 - 0.745 / 3) / 2)))
fy_m_sd <- mean(c(((0.745 - 0.352) / 2), ((1.58 - 0.745) / 2))) # Mean and standard deviations, juvenile fertility
a_m_sd <- mean(c(((2.52 - 1.98) / 2), ((3.78 - 2.52) / 2))) # Mean and standard deviations, adult fertility
# Mean and standard deviations vector, juvenile and adult fertility:
m_sd_vec <- c(0.18 * fertility[1], 0.18 * fertility[2], a_m_sd, a_m_sd, a_m_sd, a_m_sd, a_m_sd)

# Survival
# KI cat survival
# probability of surviving from one year to the next. e.g surviving fourth year of life
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)

# survival errors based on Budke & Slater
y1_2_s_sd <- mean(c(((0.46 - 0.27) / 2), ((0.73 - 0.46) / 2))) # mean and standard deviations, juvenile survival
a_s_sd <- mean(c(((0.7 - 0.55) / 2), ((0.78 - 0.7) / 2))) # mean and standard deviations, adult survival
# Mean and standard deviations vector, juvenile and adult survival:
s_sd_vec <- c(y1_2_s_sd, y1_2_s_sd, a_s_sd, a_s_sd, a_s_sd, a_s_sd)

# create matrix
popmat <- matrix(data = 0, nrow = age_max, ncol = age_max)
diag(popmat[2:age_max, ]) <- survival_probability
popmat[age_max, age_max] <- 0
popmat[1, ] <- fertility
popmat_orig <- popmat # save original matrix

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
t <- (yr_end - yr_now) # timeframe

tot_f <- sum(popmat_orig[1, ])
popmat <- popmat_orig # resets matrix

# set population storage matrices
n_mat <- matrix(0, nrow = age_max, ncol = (t + 1)) # empty matrix
n_mat[, 1] <- init_vec # fill first matrix column with initial population vector

# set up projection loop
for (i in 1:t) {
  n_mat[, i + 1] <- popmat %*% n_mat[, i]
}

# Number of predators - cats - through time period, no density reduction treatment, no carry capacity
n_pred <- colSums(n_mat)
yrs <- seq(yr_now, yr_end, 1)
predators <- tibble(yrs, n_pred)
ggplot(data = predators, aes(x = yrs, y = n_pred)) +
  geom_point(shape=19) +
  geom_line(linetype="dashed") +
  labs(x = "year", y = "N")
ggsave("reports/figures/time_serie_predators.jpg")

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
n_mat <- matrix(0, nrow = age_max, ncol = (t + 1))
n_mat[, 1] <- init_vec
popmat <- popmat_orig

# set up projection loop
for (i in 1:t) {
  tot_n_i <- sum(n_mat[, i])
  pred_red <- FeralCatEradication::survival_modifier(tot_n_i, coefficients)
  diag(popmat[2:age_max, ]) <- survival_probability * pred_red
  popmat[age_max, age_max] <- 0
  n_mat[, i + 1] <- popmat %*% n_mat[, i]
}

n_pred <- colSums(n_mat)
capacity <- tibble(yrs, n_pred)
# Untreated population increases, rate of increase relative to K, no stochastic sampling:
ggplot(data=capacity, aes(yrs, n_pred)) +
  geom_point() +
  labs(x = "year", y = "N") +
  lims(y = c(0, 1.05 * k_max)) +
  geom_hline(yintercept=k_max, linetype="dashed", color = "red")
ggsave("reports/figures/something_with_Carry_capacity.jpg")
