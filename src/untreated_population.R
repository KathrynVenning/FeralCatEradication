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

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t + 1)) # storage matrix

for (e in 1:iter) {
  popmat <- popmat.orig

  n.mat <- matrix(0, nrow = age.max, ncol = (t + 1))
  n.mat[, 1] <- init.vec

  for (i in 1:t) {
    # stochastic survival values
    s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
    s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)

    # stochastic fertilty sampler (gaussian)
    fert.stch <- rnorm(length(popmat[, 1]), popmat[1, ], m.sd.vec)
    fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)

    totN.i <- sum(n.mat[, i])
    pred.red <- a.lp / (1 + (totN.i / b.lp)^c.lp)

    popmat[1, ] <- fert.stoch
    diag(popmat[2:age.max, ]) <- s.stoch * pred.red
    # popmat[age.max,age.max] <- 0

    n.mat[, i + 1] <- popmat %*% n.mat[, i]
  } # end i loop

  n.sums.mat[e, ] <- ((as.vector(colSums(n.mat)) / pop.found))

  if (e %% itdiv == 0) print(e)
} # end e loop

n.md <- apply(n.sums.mat, MARGIN = 2, median, na.rm = T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN = 2, quantile, probs = 0.975, na.rm = T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN = 2, quantile, probs = 0.025, na.rm = T) # lower over all iterations

plot(yrs, n.md, type = "l", main = "Min N with SD for untr pop", xlab = "year", ylab = "Minimum population", lwd = 2, ylim = c(0.95 * min(n.lo), 1.05 * max(n.up)))
lines(yrs, n.lo, lty = 2, col = "red", lwd = 1.5)
lines(yrs, n.up, lty = 2, col = "red", lwd = 1.5)

untreated <- data.frame(yrs, n.md, n.lo, n.up)
