library(comprehenr)
return_one <- function() {
  return(1)
}

monthly_matrix_leslie <- function(fertility, survival) {
  fertility <- comprehenr::to_vec(for (f in fertility) rep(f / 12, 12))
  survival_probability <- comprehenr::to_vec(for (s in survival) rep(s^(1 / 12), 12))
  survival_probability <- append(survival_probability, rep(survival_probability[length(survival_probability)], 11))
  ml <- matrix_leslie(fertility, survival_probability)
  return(ml)
}
