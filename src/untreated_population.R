library(FeralCatEradication)
source("src/parameters_of_fertility_and_survival.R")
####################################################
## iterations and quasi ext for each following model
####################################################
iter <- 10000 # final model run at 10 000
itdiv <- iter / 1000 # final model rate at iter/1000
################################################################################################################
## untreated population
###############################################################################################################
## stochatic projection with density feedback
## set storage matrices & vectors

initial_population <- 1629
capacity <- Carry_Capacity$new()
coefficients <- capacity$coefficients_model(half_capacity = initial_population)
survival <- Stochastic_Survival_Fertility$new(fertility, survival_probability)
survival$set_standard_desviations(std_fertility, std_survival_probability)
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date
number_year <- yr_end - yr_now + 1
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = number_year)
population <- Population$new(survival)
for (simulation in seq(1, iter)) {
  population$run_generations(yr_now, yr_end, initial_population = initial_population, coefficients)
  n.sums.mat[simulation, ] <- colSums(population$n_mat) / initial_population
}

yrs <- seq(yr_now, yr_end)
n.md <- apply(n.sums.mat, MARGIN = 2, median, na.rm = T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN = 2, quantile, probs = 0.975, na.rm = T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN = 2, quantile, probs = 0.025, na.rm = T) # lower over all iterations
plot(yrs, n.md, type = "l", main = "Min N with SD for untr pop", xlab = "year", ylab = "Minimum population", lwd = 2, ylim = c(0.95 * min(n.lo), 1.05 * max(n.up)))
lines(yrs, n.lo, lty = 2, col = "red", lwd = 1.5)
lines(yrs, n.up, lty = 2, col = "red", lwd = 1.5)
untreated <- data.frame(yrs, n.md, n.lo, n.up)
