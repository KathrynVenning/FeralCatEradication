library(FeralCatEradication)
source("R/monthly_matrix_leslie.R")
source("R/untreated_population.R")

fertility <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)
survival_probability <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)
initial_population <- 1629
yr_now <- 2020 # update if more data available post-2010
yr_end <- 2030 # set projection end date
interval_time <- Monthly_Interval_Time$new(initial_year = yr_now, final_year = yr_end)

ml <- matrix_leslie(fertility, survival_probability)
mml <- monthly_matrix_leslie(fertility, survival_probability)
max_lambda(ml)
max_lambda(mml)^12
survival <- Monthly_Survival_Fertility$new(fertility, survival_probability)
population <- Population$new(survival)
population$run_generations(interval_time, initial_population = initial_population)

plotter <- Plotter_Population$new()
plotter$plot(population)
plotter$save("reports/figures/monthly_time_serie_individuals.jpg")
