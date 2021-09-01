# harvest rate 200-280
harv.prop.consist <- seq(0.2, 0.99, 0.05) # sequence harvest/culling quotas, e.g remove 0.5-.99 porportion of founding pop

# define our quasi-extinction probability storage vector
min.med.n <- min.lo.n <- min.up.n <- rep(0, length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {

  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t + 1))

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

      ## harvest things here
      n.mat[, i + 1] <- n.mat[, i + 1] - round(stable.stage.dist(popmat) * round(sum(n.mat[, i + 1]) * harv.prop.consist[s], 0), 0)


      if (length(which(n.mat[, i + 1] < 0)) > 0) {
        n.mat[which(n.mat[, i + 1] < 0), i + 1] <- 0
      }
    } # end i loop

    n.sums.mat[e, ] <- as.vector((colSums(n.mat)) / pop.found) # / pop.mat for min proportion remaining population

    if (e %% itdiv == 0) print(e)
  } # end e loop

  # calculate minimum population size

  min.pop.vec <- apply(n.sums.mat, MARGIN = 1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm = T)
  min.lo.n[s] <- quantile(min.pop.vec, probs = 0.025, na.rm = T)
  min.up.n[s] <- quantile(min.pop.vec, probs = 0.975, na.rm = T)

  n.md <- apply((n.sums.mat), MARGIN = 2, mean, na.rm = T) # minimum over all iterations
  n.up <- apply((n.sums.mat), MARGIN = 2, quantile, probs = 0.975, na.rm = T) # upper over all iterations
  n.lo <- apply((n.sums.mat), MARGIN = 2, quantile, probs = 0.025, na.rm = T) # lower over all iterations

  plot(yrs, n.md, type = "l", xlab = "year", ylab = "minimum N", lwd = 2, ylim = c(0.95 * min(n.lo), 1.05 * max(n.up)))
  lines(yrs, n.lo, lty = 2, col = "red", lwd = 1.5)
  lines(yrs, n.up, lty = 2, col = "red", lwd = 1.5)

  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep = ""))
  print("##############")
} # ends S loop

plot(harv.prop.consist, min.med.n, type = "l", pch = 19, xlab = "harvest proportion", ylab = "min N", ylim = c(min(min.lo.n), max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col = "red", lty = 2)
lines(harv.prop.consist, min.up.n, col = "red", lty = 2)

minn.prop.pop <- data.frame(harv.prop.consist, min.med.n, min.lo.n, min.up.n)
