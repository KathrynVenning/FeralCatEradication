library(comprehenr)
library(tidyverse)

monthly_matrix_leslie <- function(fertility, survival) {
  fertility <- comprehenr::to_vec(for (f in fertility) rep(f / 12, 12))
  survival_probability <- comprehenr::to_vec(for (s in survival) rep(s^(1 / 12), 12))
  survival_probability <- append(survival_probability, rep(survival_probability[length(survival_probability)], 11))
  ml <- matrix_leslie(fertility, survival_probability)
  return(ml)
}

Monthly_Population <- R6::R6Class("Monthly_Population",
  public = list(
    fertility = NULL,
    survival_probability = NULL,
    n_mat = NULL,
    sequence_years = NULL,
    initialize = function(fertility, survival_probability) {
      self$fertility <- fertility
      self$survival_probability <- survival_probability
    },
    run_generations = function(initial_year, final_year, initial_population, coefficients = list(a_lp = 2, b_lp = 4, c_lp = 0)) {
      n_mat <- private$setup_variables(initial_year, final_year, initial_population)
      for (year in 1:private$years) {
        tot_n_i <- sum(n_mat[, year])
        modified_survival_probability <- modifier_survival_probability(tot_n_i, coefficients, survival_probability)
        popmat <- monthly_matrix_leslie(self$fertility, modified_survival_probability)
        n_mat[, year + 1] <- popmat %*% n_mat[, year]
      }
      self$n_mat <- n_mat
    }
  ),
  private = list(
    years = NULL,
    setup_variables = function(initial_year, final_year, initial_population) {
      private$setup_temporal_variables(initial_year, final_year)
      n_mat <- private$setup_matrix_population(initial_population)
      return(n_mat)
    },
    setup_temporal_variables = function(initial_year, final_year) {
      private$years <- (final_year - initial_year) * 12
      self$sequence_years <- seq(initial_year, final_year, 1 / 12)
    },
    setup_matrix_population = function(initial_population) {
      month_max <- length(self$fertility) * 12
      n_mat <- matrix(0, nrow = month_max, ncol = (private$years + 1))
      popmat <- monthly_matrix_leslie(self$fertility, self$survival_probability)
      ssd <- FeralCatEradication::stable_stage_dist(popmat)
      classes_age_population <- ssd * initial_population
      n_mat[, 1] <- classes_age_population
      return(n_mat)
    }
  )
)
