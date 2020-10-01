###### Supplementary code for 
##### Predicting feral cat-reduction targets and costs on large islands using stochastic population models
###### can only be used after running lines 9-136 of main code

rm(list = ls())

#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 10000 #final model run at 10 000
itdiv <- iter/100 #final model rate at iter/1000

################################################################################################################
## untreated population with leakage
###############################################################################################################
## stochatic projection with density feedback
## set storage matrices & vectors
stray.cat.vec <- seq(0,100,10) #stray cats added 0 - 100 cats increasing by 10
final.md.out <- final.n.up.out <- final.n.lo.out <- rep(NA,length(stray.cat.vec)) #storage matrix

for (s in 1:length(stray.cat.vec)) {
  
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix
  
  for (e in 1:iter) {
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertilty sampler (gaussian)
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
      diag(popmat[2:age.max,]) <- s.stoch[-age.max]*pred.red
      popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
      
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      n.mat[,i+1] <- n.mat[,i+1] + round(rnorm(1,mean=stray.cat.vec[s],sd=(0.05*stray.cat.vec[s]))*ssd, 0) #adding stray cats into the population vector
      
    } # end i loop
    
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat))/pop.found))
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  n.md <-  apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  
  final.md.out[s] <- n.md[t+1]
  final.n.up.out[s] <- n.up[t+1]
  final.n.lo.out[s] <- n.lo[t+1]

  print(s)
  
} # s vec

  untreated.leakage <- matrix(data=0, nrow = 4, ncol = length(yrs)) #results, proportional pop increase relative to year, empty matrix
  untreated.leakage[1,] <- yrs #fill matrix, column 1 = years
  untreated.leakage[2,] <- final.md.out #column 2 = median final N
  untreated.leakage[3,] <- final.n.up.out #column 3 = upper 95% confidence 
  untreated.leakage[4,] <- final.n.lo.out# column 4 = lower 95% confidence 

  plot(stray.cat.vec, n.md, pch=19, type="l", xlab="mean number of stray cats added/year", ylab="final N at end of projection interval",  lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up))) #plot
  lines(stray.cat.vec,final.n.lo.out,lty=2,col="red",lwd=1.5)
  lines(stray.cat.vec,final.n.up.out,lty=2,col="red",lwd=1.5)
  
  

##############################################################################################################
### main two-phase cull with leakage
###############################################################################################################
## stochatic projection with density feedback
## set storage matrices & vectors
stray.cat.vec <- seq(0,100,10) #stray cats added 0 - 100 cats increasing by 10
final.md.out <- rep(NA,length(stray.cat.vec)) #storage matrix

for (s in 1:length(stray.cat.vec)) {
  
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix
  
  for (e in 1:iter) {
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertilty sampler (gaussian)
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
      diag(popmat[2:age.max,]) <- s.stoch[-age.max]*pred.red
      popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
      
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      n.mat[,i+1] <- n.mat[,i+1] + round(rnorm(1,mean=stray.cat.vec[s],sd=(0.05*stray.cat.vec[s]))*ssd, 0) #add stray cats into population vector
      
      #harvest
      if (i <= 3) {
        n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*0.6, 0), 0) 
      } else {
        n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*0.45, 0), 0)
      }
      
    } # end i loop
    
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat))/pop.found))
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  final.md.out[s] <- n.md[t+1]
  final.n.up.out[s] <- n.up[t+1]
  final.n.lo.out[s] <- n.lo[t+1]
  print(s)
  
} # s vec

  plot(stray.cat.vec, final.md.out, pch=19, type="l", xlab="mean number of stray cats added/year", ylab="final N at end of projection interval")
  lines(stray.cat.vec, final.n.up.out, lty=2, col="red")
  lines(stray.cat.vec, final.n.lo.out, lty=2 , col="red")
  
totalN <- final.md.out*pop.found 

final.md.table <- matrix(0, nrow = 5, ncol = length(stray.cat.vec)) #results matrix
final.md.table[1,] <- stray.cat.vec #how many stray cats added 
final.md.table[2,] <- final.md.out #final median N proportionate
final.md.table[3,] <- totalN #final median N
final.md.table[4,] <- final.n.lo.out #final lower 95% confidence
final.md.table[5,] <- final.n.up.out #final upper 95% confidence
final.md.table #display table

###############################################################################################################################
##stopping early
###############################################################################################################################

##############################################################################################################
## main two-phase cull 
###############################################################################################################
## stochatic projection with density feedback
## set storage matrices & vectors

n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1)) #storage matrix

for (e in 1:iter) {
  popmat <- popmat.orig
  
  n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
  n.mat[,1] <- init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
    s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
    s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
    
    # stochastic fertilty sampler (gaussian)
    fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
    fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
    
    totN.i <- sum(n.mat[,i])
    pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
    diag(popmat[2:age.max,]) <- s.stoch[-age.max]*pred.red
    popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
    
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    
    #harvest 
    if (i <3 ) {
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*0.6, 0), 0)
    } else {
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*0.45, 0), 0)
    }
    
    if (length(which(n.mat[,i+1] < 0)) > 0) {
      n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
    }
  } # end i loop
  
  n.sums.mat[e,] <- ((as.vector(colSums(n.mat))/pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(yrs,n.md,type="l", main = "Min N with SD for untr pop", xlab="year", ylab="Minimum population", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,n.up,lty=2,col="red",lwd=1.5)

two.phase <- matrix(data=0, nrow = 2, ncol = length(yrs)) #results, proportional pop relative to year
two.phase[1,] <- yrs
two.phase[2,] <- n.md
two.phase 

two.phase.N <- matrix (data = 0, nrow = 2, ncol = length(yrs)) #results, total pop relative to year
two.phase.N[1,] <- yrs
two.phase.N[2,] <- n.md*pop.found
two.phase.N


#################################################### 
## Stopping culling early 
####################################################
iter <- 10000
itdiv <- iter/10
stopped.yrs.vec <- seq(3,11,1) #stopping cull after years 3-11, intervals of 1
t <- seq(2020,2100,1) #time

recov.med <- recov.lo <- recov.up <- rep(NA,length(stopped.yrs.vec)) #storage

for (s in 1:length(stopped.yrs.vec)) {
  
  recovery.time <- rep(NA,iter) #storage for loop
  
  for (e in 1:iter) {
    
    n.mat <- matrix(0, nrow=age.max, ncol=(length(t)+1))
    n.mat[,1] <- init.vec
    popmat <- popmat.orig
    
    #for (i in 1:4) {
    for (i in 1:length(t)) {
      # stochastic survival values
      s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
      s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
      s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
      
      # stochastic fertility sampler (Gaussian)
      fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
      fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
      
      totN.i <- sum(n.mat[,i])
      pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
      diag(popmat[2:age.max,]) <- s.stoch[-age.max]*pred.red
      popmat[age.max,age.max] <- s.stoch[age.max]*pred.red
      
      # projection
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      #harvest
      if (i <= 11) {
        harv.initial <- round(stable.stage.dist(popmat.orig) * round(sum(n.mat[,i+1])*0.6, 0), 0)
        if (i <= 2) {
          n.mat[,i+1] <- ifelse(n.mat[,i+1] - harv.initial < 0, 0, n.mat[,i+1] - harv.initial) #apply initial harvest
        }
        if (i > 2) {
          harv.maint <- round(stable.stage.dist(popmat.orig) * round(sum(n.mat[,i+1])*0.45, 0), 0)
          if (i <= stopped.yrs.vec[s]) {
            n.mat[,i+1] <- ifelse(n.mat[,i+1] - harv.maint < 0, 0, n.mat[,i+1] - harv.maint)}
        } #stop maintenance harvest at years 3-11 (stooped.yrs.vec)
      } # end if
      
    } # end i loop
    
    recovery.time[e] <- which(colSums(n.mat[,-1]) >= sum(init.vec))[1]
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  recov.med[s] <- median(recovery.time)
  recov.lo[s] <- quantile(recovery.time, probs=0.025, na.rm=T)
  recov.up[s] <- quantile(recovery.time, probs=0.975, na.rm=T)
  
  print("#################")
  print(paste("stop year = ", s, sep=""))
  print("#################")
  
}

plot(stopped.yrs.vec, recov.med, type="l", xlab="stopped year", ylab="years to initial recovery")
lines(stopped.yrs.vec, recov.lo, lty=2, col="red")
lines(stopped.yrs.vec, recov.up, lty=2, col="red")
