#' @export
max_lambda <- function(x) {
  Re((eigen(x)$values)[1])
}

#' @export
max_r <- function(x) {
  log(max_lambda(x))
}

#' @export
stable_stage_dist <- function(x) {
  real_first_eigen_vector <- Re((eigen(x)$vectors)[, 1])
  parallel_matrix <- x %*% real_first_eigen_vector
  (parallel_matrix / (sum(parallel_matrix)))[, 1]
}

# Generation length function
# reproductive value (r_0) where leslie_matrix = Leslie matrix; age_max = maximum age of females
#' @export
total_female_offspring_per_female <- function(leslie_matrix, age_max) {
  # define the transition matrix
  transition_matrix <- leslie_matrix[1:age_max, 1:age_max]
  transition_matrix[1, 1:(age_max)] <- 0

  # define the fertility matrix
  fertility_matrix <- leslie_matrix[1:age_max, 1:age_max]
  diag(fertility_matrix[2:age_max, 1:(age_max - 1)]) <- 0

  # define the identity matrix
  identity_matrix <- matrix(data = 0, nrow = age_max, ncol = age_max)
  diag(identity_matrix) <- 1

  # define the fundamental matrix
  n_fund <- MASS::ginv(identity_matrix - transition_matrix)

  # define the reproductive matrix
  reproductive_matrix <- fertility_matrix %*% n_fund

  # define r_0 (number of female offspring produced per female during lifetime)
  r_0 <- Re((eigen(reproductive_matrix)$values)[1])

  return(r_0)
}

# Mean generation time function
# where leslie_matrix is a Leslie Matrix
#' @export
g_val <- function(leslie_matrix, age_max) {
  mean_generation_time <- (log(total_female_offspring_per_female(leslie_matrix, age_max))) / (log(Re((eigen(leslie_matrix)$values)[1])))
  return(mean_generation_time)
}

#' @export
coefficients_proportion_realized_survival <- function(k_vec, red_vec) {
  k_red_dat <- data.frame(k_vec, red_vec)
  param_init <- c(1, 15000, 2.5)
  fit_lp <- nls(red_vec ~ a / (1 + (k_vec / b)^c),
    data = k_red_dat,
    algorithm = "port",
    start = c(a = param_init[1], b = param_init[2], c = param_init[3]),
    trace = TRUE,
    nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1 / 1024)
  )
  coefficients <- clean_coefficients(fit_lp)
}

clean_coefficients <- function(fit_lp) {
  a_lp <- coef(fit_lp)[1]
  b_lp <- coef(fit_lp)[2]
  c_lp <- coef(fit_lp)[3]
  names(a_lp) <- NULL
  names(b_lp) <- NULL
  names(c_lp) <- NULL
  coefficients <- list(a_lp = a_lp, b_lp = b_lp, c_lp = c_lp)
}

#' @export
survival_modifier <- function(tot_n_i, coefficients) {
  pred_red <- coefficients$a_lp / (1 + (tot_n_i / coefficients$b_lp)^coefficients$c_lp)
}

#' @export
modifie_survival_probability <- function(tot_n_i, coefficients, survival_probability) {
  pred_red <- FeralCatEradication::survival_modifier(tot_n_i, coefficients)
  modified_survival_probability <- survival_probability * pred_red
  return(modified_survival_probability)
}

#' @export
dont_modifie_survival_probability <- function(tot_n_i, coefficients, survival_probability) {
  return(survival_probability)
}

#' @export
matrix_leslie <- function(fertility, survival_probability) {
  age_max <- length(fertility)
  popmat <- matrix(data = 0, nrow = age_max, ncol = age_max)
  diag(popmat[2:age_max, ]) <- survival_probability
  popmat[age_max, age_max] <- 0
  popmat[1, ] <- fertility
  return(popmat)
}

#' @export
Population <- R6::R6Class("Population",
  public = list(
    fertility = NULL,
    survival_probability = NULL,
    popmat = NULL,
    n_mat = NULL,
    initialize = function(fertility, survival_probability) {
      self$fertility <- fertility
      self$survival_probability <- survival_probability
      self$popmat <- matrix_leslie(fertility, survival_probability)
    },
    run_generations = function(years, initial_population) {
      age_max <- length(self$fertility)
      n_mat <- matrix(0, nrow = age_max, ncol = (years + 1))
      n_mat[, 1] <- initial_population
      for (year in 1:years) {
        tot_n_i <- sum(n_mat[, year])
        modified_survival_probability <- survival_probability
        popmat <- matrix_leslie(fertility, modified_survival_probability)
        n_mat[, year + 1] <- popmat %*% n_mat[, year]
      }
      self$n_mat <- n_mat
    }
  ),
  private = list()
)

#' @export
Plotter_Population <- R6::R6Class("Plotter_Population",
  public = list(
    initialize = function(population) {
    },
    plot = function(years, population) {
      individuals <- private$setup_variables(years, population)
      y_ticks <- private$setup_y_ticks(individuals)
      private$make_plot(individuals, y_ticks)
    },
    save = function(path) {
      ggsave(path)
    }
  ),
  private = list(
    setup_variables = function(years, population) {
      n_pred <- colSums(population$n_mat)
      individuals <- tibble(yrs = as.character(years), n_pred)
      return(individuals)
    },
    setup_y_ticks = function(individuals) {
      marcasEjeY <- pretty(c(0, max(individuals$n_pred)))
      return(marcasEjeY)
    },
    make_plot = function(individuals, y_ticks) {
      ggplot(data = individuals, aes(x = yrs, y = n_pred)) +
        geom_point(shape = 19) +
        geom_line(linetype = "dashed") +
        theme_classic() +
        scale_y_continuous(
          expand = c(0, 0),
          limits = range(y_ticks),
          breaks = y_ticks
        ) +
        labs(x = "", y = "Number of individuals (cats)")
    }
  )
)
