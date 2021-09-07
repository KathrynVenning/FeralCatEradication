# harvest rate 200-280
harv.prop.consist <- seq(0.2, 0.99, 0.05) # sequence harvest/culling quotas, e.g remove 0.5-.99 porportion of founding pop

# define our quasi-extinction probability storage vector
min.med.n <- min.lo.n <- min.up.n <- rep(0, length(harv.prop.consist))

library("ggplot2")
library(FeralCatEradication)
source("R/feral_cat.R")
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
interval_time <- Interval_Time$new(initial_year = yr_now, final_year = yr_end)
number_year <- yr_end - yr_now + 1
n_sums_mat <- matrix(data = 0, nrow = iter, ncol = number_year)
population <- Population$new(survival)

for (s in 1:length(harv.prop.consist)) {

  # set storage matrices & vectors
  simulator2 <- Runner_Population_With_CC_harvest$new(population, coefficients, harv.prop.consist[s])
  for (simulation in seq(1, iter)) {
    simulator2$run_generations(interval_time, initial_population = initial_population)
    n_sums_mat[simulation, ] <- colSums(simulator2$n_mat) / initial_population
  }
  min.pop.vec <- apply(n_sums_mat, MARGIN = 1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm = T)
  min.lo.n[s] <- quantile(min.pop.vec, probs = 0.025, na.rm = T)
  min.up.n[s] <- quantile(min.pop.vec, probs = 0.975, na.rm = T)

  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep = ""))
  print("##############")
}

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)
ggplot(minn.prop.pop, aes(x = harv.prop.consist, y = min.med.n)) +
  geom_line(colour = "blue") +
  geom_ribbon(aes(ymin = min.lo.n, ymax = min.up.n), alpha = 0.2) +
  labs(x = "constant proportional cull", y = "proportion of N1")
ggsave("reports/figures/constant_proportional_annual_cull.jpg")
