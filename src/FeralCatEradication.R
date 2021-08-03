# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly

# remove everything
rm(list = ls())

# libraries
library(plotly)
library(FeralCatEradication)
source("src/matrixOperators.R")
options(scipen = 1000)

# functions
# beta distribution shape parameter estimator function
est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# create Leslie matrix
age_max <- 7

# Create vectors
# Fertility
# KI cat birth rates matrix, data for female offsping produced each year.
# Data from Budke, C & Slater, M (2009)
m_vec <- c((0.745 / 3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98)

# Fertility errors based on Budke & Slater
# Mean and standard deviations, juvenile fertility:
juv_m_sd <- mean(c(((0.745 / 3 - 0.352 / 3) / 2), ((1.58 / 3 - 0.745 / 3) / 2)))
fy_m_sd <- mean(c(((0.745 - 0.352) / 2), ((1.58 - 0.745) / 2))) # Mean and standard deviations, juvenile fertility
a_m_sd <- mean(c(((2.52 - 1.98) / 2), ((3.78 - 2.52) / 2))) # Mean and standard deviations, adult fertility
# Mean and standard deviations vector, juvenile and adult fertility:
m_sd_vec <- c(0.18 * m_vec[1], 0.18 * m_vec[2], a_m_sd, a_m_sd, a_m_sd, a_m_sd, a_m_sd)

# Survival
# KI cat survival
# probability of surviving from one year to the next. e.g surviving fourth year of life
s_vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7)

# survival errors based on Budke & Slater
y1_2_s_sd <- mean(c(((0.46 - 0.27) / 2), ((0.73 - 0.46) / 2))) #mean and standard deviations, juvenile survival
a_s_sd <- mean(c(((0.7 - 0.55) / 2), ((0.78 - 0.7) / 2))) #mean and standard deviations, adult survival
# Mean and standard deviations vector, juvenile and adult survival:
s_sd_vec <- c(y1_2_s_sd, y1_2_s_sd, a_s_sd, a_s_sd, a_s_sd, a_s_sd)

# create matrix
popmat <- matrix(data = 0, nrow = age_max, ncol = age_max)
diag(popmat[2:age_max, ]) <- s_vec
popmat[age_max, age_max] <- 0
popmat[1, ] <- m_vec
popmat_orig <- popmat # save original matrix

# matrix properties
FeralCatEradication::max_lambda(popmat) # 1-yr lambda
FeralCatEradication::max_r(popmat) # rate of population change, 1-yr
FeralCatEradication::stable_stage_dist(popmat) # stable stage distribution
FeralCatEradication::total_female_offspring_per_female(popmat, age_max) # reproductive value
gen_l <- g_val(popmat, age_max) # mean generation length

# initial population vector
pop_found <- 1629 # +/- 661 founding population size Hohnen et al 2020
ssd <- FeralCatEradication::stable_stage_dist(popmat)
init_vec <- ssd * pop_found #initial population vector

# Project
# set time limit for projection in 1-yr increments
yr_now <- 2020 # update if more data available post-2010
#************************
yr_end <- 2030 # set projection end date
#************************
t <- (yr_end - yr_now) #timeframe

tot_f <- sum(popmat_orig[1, ])
popmat <- popmat_orig #resets matrix
yr_vec <- seq(yr_now, yr_end) #year vector, 2020, 2021, 2022...

# set population storage matrices
n_mat <- matrix(0, nrow = age_max, ncol = (t + 1)) #empty matrix
n_mat[, 1] <- init_vec #fill first matrix column with initial population vector

# set up projection loop
for (i in 1:t) {
  n_mat[, i + 1] <- popmat %*% n_mat[, i]
}

# Number of predators - cats - through time period, no density reduction treatment, no carry capacity
n_pred <- colSums(n_mat)
yrs <- seq(yr_now, yr_end, 1)
plot(yrs, n_pred, type = "b", lty = 2, pch = 19, xlab = "year", ylab = "N")

# Compensatory density feedback
# K = carry capacity
# Population rate of increase relative to carry capacity
# Larger distance between populationa and K = faster population growth
k_max <- 2 * pop_found
k_min <- 1 #not used
k_vec <- c(1, pop_found / 2, pop_found, 0.75 * k_max, k_max) #1= k_min, .75 = red_thresh??
red_thresh <- 0.75 #not used
red_vec <- c(1, 0.965, 0.89, 0.79, 0.71)
jpeg("reports/figures/k_vec.jpg")
plot(k_vec, red_vec, pch = 19, type = "b")
dev.off()
k_red_dat <- data.frame(k_vec, red_vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
param_init <- c(1, 15000, 2.5)
fit_lp <- nls(red_vec ~ a / (1 + (k_vec / b) ^ c),
              data = k_red_dat,
              algorithm = "port",
              start = c(a = param_init[1], b = param_init[2], c = param_init[3]),
              trace = TRUE,
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1 / 1024))
fit_lp_summ <- summary(fit_lp)
jpeg("reports/figures/reduction_factor.jpg")
plot(k_vec, red_vec, pch = 19, xlab = "N", ylab = "reduction factor")
dev.off()
k_vec_cont <- seq(1, 2 * pop_found, 1)
pred_lp_fx <- coef(fit_lp)[1] / (1 + (k_vec_cont / coef(fit_lp)[2]) ^ coef(fit_lp)[3])
lines(k_vec_cont, pred_lp_fx, lty = 2, lwd = 3, col = "red")

a_lp <- coef(fit_lp)[1]
b_lp <- coef(fit_lp)[2]
c_lp <- coef(fit_lp)[3]

print(a_lp)
print(b_lp)
print(c_lp)

# compensatory density-feedback deterministic model
# set population storage matrices
n_mat <- matrix(0, nrow = age_max, ncol = (t + 1))
n_mat[, 1] <- init_vec
popmat <- popmat_orig

# set up projection loop
for (i in 1:t) {
  tot_n_i <- sum(n_mat[, i])
  pred_red <- a_lp / (1 + (tot_n_i / b_lp) ^ c_lp)
  diag(popmat[2:age_max, ]) <- s_vec * pred_red
  popmat[age_max, age_max] <- 0
  n_mat[, i + 1] <- popmat %*% n_mat[, i]
}

n_pred <- colSums(n_mat)
jpeg("reports/figures/something_with_Carry_capacity.jpg")
# Untreated population increases, rate of increase relative to K, no stochastic sampling:
plot(yrs, n_pred, type = "b", lty = 2, pch = 19, xlab = "year", ylab = "N", ylim = c(0, 1.05 * k_max))
abline(h = k_max, lty = 2, col = "red") #carry capacity
legend(yrs[2], n_pred[6], legend = c("N", "Carry capacity"),
       col = c("black", "red"), lty = 1:2, cex = 0.8)
dev.off()
