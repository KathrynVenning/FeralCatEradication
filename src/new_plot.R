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

yrs <- interval_time$get_time_sequence()
for (s in 1:length(harv.prop.consist)) {
  
  # set storage matrices & vectors
  simulator2 <- Runner_Population_With_CC_harvest$new(population, coefficients, harv.prop.consist[s])
  for (simulation in seq(1, iter)) {
    simulator2$run_generations(interval_time, initial_population = initial_population)
    n_sums_mat[simulation, ] <- colSums(simulator2$n_mat) / initial_population
  }
  min.pop.vec <- apply(n_sums_mat, MARGIN=1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  n.md <- apply((n_sums_mat), MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply((n_sums_mat), MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply((n_sums_mat), MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  print("##############")
  
} # ends S loop

plot(harv.prop.consist, min.med.n, type="l", pch=19, xlab="harvest proportion", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col="red", lty=2)
lines(harv.prop.consist, min.up.n, col="red", lty=2)

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)
