#' @export
est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

get_stochastic_fertility <- function(fertility, sd_fertility) {
  fert.stch <- rnorm(length(fertility), fertility, sd_fertility)
  fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
  return(fert.stoch)
}

get_stochastic_survival <- function(survival, sd_survival) {
  parameter <- estBetaParams(survival, sd_survival^2)
  s.stoch <- rbeta(length(survival), parameter[["alpha"]], parameter[["beta"]])
  return(s.stoch)
}

Survival_Fertility <- R6::R6Class("Survival_Fertility",
  public = list(
    initialize = function(fertility, survival_probability) {
      private$survival_probability <- survival_probability
      private$fertility <- fertility
    },
    get_fertility = function() {
      return(private$fertility)
    }
  ),
  private = list(
    fertility = NULL,
    survival_probability = NULL
  )
)
