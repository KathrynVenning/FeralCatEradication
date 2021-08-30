library(comprehenr)
library(tidyverse)

#' @export
monthly_matrix_leslie <- function(fertility, survival) {
  fertility <- comprehenr::to_vec(for (f in fertility) rep(f / 12, 12))
  survival_probability <- comprehenr::to_vec(for (s in survival) rep(s^(1 / 12), 12))
  survival_probability <- append(survival_probability, rep(survival_probability[length(survival_probability)], 11))
  ml <- matrix_leslie(fertility, survival_probability)
  return(ml)
}

Monthly_Population <- R6::R6Class("Monthly_Population",
  inherit = Population,
  public = list(),
  private = list(
    setup_temporal_variables = function(initial_year, final_year) {
      private$years <- (final_year - initial_year) * 12
      self$sequence_years <- seq(initial_year, final_year, 1 / 12)
    }
  )
)
