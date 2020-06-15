# Kathryn Venning 
# Flinders University
# feral cat reduction on Kangaroo Island
# requires library - Plotly

## Remove everything
rm(list = ls())

# libraries
library(plotly)

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
#available on GitHub with this file, labelled "matrix operators"

# create Leslie matrix
age.max = 7

## create vectors 
#fertility 
m.vec <- c(0, 0.745, 2.52, 2.52, 2.52, 2.52, 1.98) ## KI cat birth rates matrix, data for female offsping produced each year. Data from Budke, C & Slater, M (2009)

# fertility errors based on Budke & Slater
juv.m.sd <- mean(c(((0.745-0.352)/2),((1.58-0.745)/2))) #mean and standard deviations, juvenile fertility
A.m.sd <- mean(c(((2.52-1.98)/2),((3.78-2.52)/2))) #mean and standard deviations, adult fertility
m.sd.vec <- c(0,juv.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd) #mean and standard deviations vector, juvenile and adult fertility 

#survival 
s.vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7, 0.55) ##KI cat survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors based on Budke & Slater
y1.2.S.sd <- mean(c(((0.46-0.27)/2),((0.73-0.46)/2))) #mean and standard deviations, juvenile survival
A.S.sd <- mean(c(((0.7-0.55)/2),((0.78-0.7)/2))) #mean and standard deviations, adult survival
s.sd.vec <- c(y1.2.S.sd,y1.2.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd) #mean and standard deviations vector, juvenile and adult survival

# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max)
diag(popmat[2:age.max,]) <- s.vec[-age.max]
popmat[age.max,age.max] <- s.vec[age.max]
popmat[1,] <- m.vec
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 3083 # founding population size
init.vec <- stable.stage.dist(popmat) * pop.found #initial population vector

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2020 # update if more data available post-2010
#************************
yr.end <- 2030 # set projection end date
#************************
t <- (yr.end - yr.now) #timeframe

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig #resets matrix 
yr.vec <- seq(yr.now,yr.end) #year vector, 2020, 2021, 2022...

## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1)) #empty matrix
n.mat[,1] <- init.vec #fill first matrix column with initial population vector

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat) #number of predators - cats - through time period, no density reduction treatment, no carry capacity
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback # K = carry capacity
#population rate of increase relative to carry capacity. Larger distance between populationa and K = faster population growth
K.max <- 2*pop.found
K.min <- 1 #not used
K.vec <- c(1,pop.found/2,pop.found,0.75*K.max,K.max) #1= K.min, .75 = red.thresh??
red.thresh <- 0.75 #not used
red.vec <- c(1,0.99,0.95,0.875,0.75)
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

# logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
param.init <- c(1, 15000, 2.5)
fit.lp <- nls(red.vec ~ a/(1+(K.vec/b)^c), 
              data = Kred.dat,
              algorithm = "port",
              start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
              trace = TRUE,      
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.lp.summ <- summary(fit.lp)
plot(K.vec,red.vec,pch=19,xlab="N",ylab="reduction factor")
K.vec.cont <- seq(1,2*pop.found,1)
pred.lp.fx <- coef(fit.lp)[1]/(1+(K.vec.cont/coef(fit.lp)[2])^coef(fit.lp)[3])
lines(K.vec.cont,pred.lp.fx,lty=2,lwd=3,col="red")

a.lp <- coef(fit.lp)[1]
b.lp <- coef(fit.lp)[2]
c.lp <- coef(fit.lp)[3]


## compensatory density-feedback deterministic model
## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## set up projection loop
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
  diag(popmat[2:age.max,]) <- s.vec[-age.max]*pred.red
  popmat[age.max,age.max] <- s.vec[age.max]*pred.red
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N") #untreated population increases, rate of increase relative to K, no stochastic sampling

#################################################### 
## iterations and quasi ext for each following model
####################################################
iter <- 10000 #final model run at 10 000
itdiv <- iter/10 #final model rate at iter/1000

# quasi-extinction threshold
q.ext <- 20 #lower threshold for viable population
prop.q.ext <- pop.found*0.004 #lower proportional threshold for viable population

################################################################################################################
## untreated population
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

untreated <- matrix(data=0, nrow = 2, ncol = length(yrs)) #results, proportional pop increase relative to year
untreated[1,] <- yrs
untreated[2,] <- n.md

###############################################################################################################################
##consistent yearly harvest
###############################################################################################################################

# harvest rate
harv.prop.consist <- seq(0.2,0.99,0.05) #sequence harvest/culling quotas, e.g remove 0.5-.99 porportion of founding pop

# define our quasi-extinction probability storage vector
q.ext.vec <-  min.med.n <- min.lo.n <- min.up.n <- rep(0,length(harv.prop.consist))

for (s in 1:length(harv.prop.consist)) {
  
  # set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
  
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
      
      ## harvest things here
      n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.consist[s], 0), 0)
      
      
      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
      }
      
    } # end i loop
    
    n.sums.mat[e,] <- as.vector((colSums(n.mat))/pop.found) # / pop.mat for min proportion remaining population 
    
    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # calculate minimum population size
  
  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  q.ext.vec[s] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
  q.ext.vec
  
  n.md <- apply((n.sums.mat), MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply((n.sums.mat), MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("harvest proportion = ", harv.prop.consist[s], sep=""))
  print("##############")
  
} # ends S loop

plot(harv.prop.consist, q.ext.vec, type="b", pch=19, xlab="harvest proportion", ylab="extinction probability")

plot(harv.prop.consist, min.med.n, type="l", pch=19, xlab="harvest proportion", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(harv.prop.consist, min.lo.n, col="red", lty=2)
lines(harv.prop.consist, min.up.n, col="red", lty=2)

minn.prop.pop <- data.frame(harv.prop.consist, (min.med.n/pop.found))


##################################################################################################################################################
## high harvest for initial 2 years, consistent harvest remaining years
#####################################################################################################################################################

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)

# storage
qext.mat <- minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init)) #storage matrices

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {

    # storage
    n.sums.mat <- p.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
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
        
        ## harvest 
        if (i < 3) {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
        } else {
          n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)

    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)  
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    # quasi-extinction
    qext.mat[n, m] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
   
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
    
} # end m loop (maintenance harvest rate)
  
## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# minimum proportional population size (median)
par(mar=c(5,5,2,8))
pminmed3d <- plot_ly(z = ~pmin.med.mat, autocontour=F, type="contour", line = list(smoothing = 0.90), contours = list(start=0.01, end=0.32, size=0.025, showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "med min pN1", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
pminmed3d

pmin3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~pmin.med.mat) %>%
  add_surface(z = ~pmin.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~pmin.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min pN1", tickfont=f3, titlefont=f1)))
pmin3d

# quasi ext (median)
par(mar=c(5,5,2,8))
minmed3d <- plot_ly(z = ~qext.mat, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "qE", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
minmed3d

min3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~minn.med.mat) %>%
  add_surface(z = ~minn.lo.mat, opacity = 0.55) %>%
  add_surface(z = ~minn.up.mat, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="min N1", tickfont=f3, titlefont=f1)))
min3d

# quasi-extinction
par(mar=c(5,5,2,8))
minmed3d <- plot_ly(z = ~pmin.med.mat, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "med min N1", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
minmed3d



##################################################################################################################################################
## COSTS
## high harvest for initial 2 years, consistent harvest remaining years
#####################################################################################################################################################
## contributed by CJA Bradshaw
###########################################################################################
## cost parameters
felixer.unit <- 13000 # AU$ #cost from felix vs felixer report 
trap.unit <- c(157,297) # AU$ # cost per trap from the internet 
shoot.ph <- 518.54/20 # ammo & labour (total AU$ over 20 hours) # Holmes et al 2015
bait.ha <- 24.87 # NZ$26/ha # Nugent et al 2012

# Felixers
# 6 felixers for Dudley (n = 262), killed 13 cats over 38 days
num.felixer <- round((6/262) * pop.found, 0)
pfelixer.killr <- (13/262 * (1/(227/365)))

Dudley.area <- 375*100 # ha
felixer.dens <- 6/Dudley.area
KI.area <- 4405*100 # ha
KI.felixer.num <- round(KI.area * felixer.dens, 0)
KI.felixer.num

# traps
# 40 traps killed 21 over 148 days
ptrap.killr <- (21/262 * (1/(148/365)))
trap.dens <- 40/Dudley.area
KI.trap.num <- round(KI.area * trap.dens, 0)
KI.trap.num

# shooting
# 14725 person-hours killed 872 (+ 172 from wounds) cats (Marion)
# assume cats not killed by Felixers & traps shot by hunters
cats.pph <- (872+172)/14725
###########################################################################################

###########################################################################################
## Type III functional response (reduction in capture efficiency with decreasing density)
max.eff <- 1 # max efficiency 
min.eff <- 0 #min efficiency
max.pN <- 1 #max population proportion
min.pN <- 0 #min population proportion
infl.eff <- 0.5

pN.vec <- c(min.pN, 0.2, 0.4, 0.5, 0.7, 0.8, max.pN) 
eff.vec <- c(min.eff, 0.05, 0.3, infl.eff, 0.85, 0.95, max.eff)
plot(pN.vec, eff.vec, type="b", pch=19)
eff.dat <- data.frame(pN.vec, eff.vec)
colnames(eff.dat) <- c("pN", "eff")

# a/(1 + b*e^(-cx)) (logistic)
param.init <- c(1, 85, 8.9)
fit.eff <- nls(eff ~ a/(1+(b*exp(-c*pN))), 
              data = eff.dat,
              algorithm = "port",
              start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
              trace = TRUE,      
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.eff.summ <- summary(fit.eff)
plot(pN.vec,eff.vec,pch=19,xlab="pN",ylab="efficiency")
pN.vec.cont <- seq(0,1,0.01)
pred.eff.fx <- coef(fit.eff)[1]/(1+(coef(fit.eff)[2]*exp(-coef(fit.eff)[3]*pN.vec.cont)))
lines(pN.vec.cont,pred.eff.fx,lty=2,lwd=3,col="red")

a.eff <- coef(fit.eff)[1]
b.eff <- coef(fit.eff)[2]
c.eff <- coef(fit.eff)[3]

###########################################################################################


## run choices
## make up shortfall in kill by ...
#shortfall.method <- "F" # adding Felixer units
#shortfall.method <- "T" # adding traps
shortfall.method <- "H" # increasing hunting pressure

# harvest rate
harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)

# storage
qext.mat <- minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- totcost.med <- totcost.lo <- totcost.up <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init))

for (m in 1:length(harv.prop.maint)) {
  
  for (n in 1:length(harv.prop.init)) {
    
    # storage
    k.sums.mat <- n.sums.mat <- p.sums.mat <- totalcost.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
    
    for (e in 1:iter) {
      
      popmat <- popmat.orig
      
      n.mat <- k.mat <- matrix(0, nrow=age.max,ncol=(t+1))
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
        
        ## harvest things here
        if (i < 3) {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
        } else {
          k.mat[,i+1] <- round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
          n.mat[,i+1] <- n.mat[,i+1] - k.mat[,i+1]
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
          n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        if (length(which(k.mat[,i+1] < 0)) > 0) {
          k.mat[which(k.mat[,i+1] < 0), i+1] <- 0
        }
        
      } # end i loop
      
      k.sums.mat[e,] <- as.vector(colSums(k.mat))
      n.sums.mat[e,] <- as.vector(colSums(n.mat))
      p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
      
      # cost of cats killed here
      eff.vec.iter <- a.eff/(1+(b.eff*exp(-c.eff*p.sums.mat[e,]))) # efficiency this iteration
      
      # calculate numbers killed per year using Felixers + traps at initial roll-out
      felixer.kill.base <- round(k.sums.mat[e,] * (eff.vec.iter*pfelixer.killr), 0)
      trap.kill.base <- round(k.sums.mat[e,] * (eff.vec.iter*ptrap.killr), 0)
      ft.kill.base <- felixer.kill.base + trap.kill.base
      shortfall <- k.sums.mat[e,] - ft.kill.base # how many cats not being killed by these methods?
      
      base.cost <- (KI.felixer.num*felixer.unit) + (KI.trap.num*runif(1,min=trap.unit[1],max=trap.unit[2])) # at initial roll-out numbers
      
      # make up shortfall
      if (shortfall.method == "H") {
        makeup.iter <- shoot.ph*(shortfall / (cats.pph*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "F") {
        makeup.iter <- felixer.unit*(shortfall / (pfelixer.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      if (shortfall.method == "T") {
        makeup.iter <- (runif(1,min=trap.unit[1],max=trap.unit[2]))*(shortfall / (ptrap.killr*eff.vec.iter)) # how many person-hours required to make up shortfall?
      }
      
      totalcost.mat[e,] <- base.cost + makeup.iter 
      
      if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)
    
    min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
    min.ppop.vec <- apply(p.sums.mat, MARGIN=1, min, na.rm=T)
    
    # median, lower & upper minimum population sizes
    minn.med.mat[n, m] <- median(min.pop.vec, na.rm=T) 
    minn.lo.mat[n, m] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
    minn.up.mat[n, m] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
    
    # median, lower & upper minimum proportional population sizes
    pmin.med.mat[n, m] <- median(min.ppop.vec, na.rm=T)  
    pmin.lo.mat[n, m] <- quantile(min.ppop.vec, probs=0.025, na.rm=T) 
    pmin.up.mat[n, m] <- quantile(min.ppop.vec, probs=0.975, na.rm=T)
    
    # quasi-extinction
    qext.mat[n, m] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
    
    ## costs
    totcost.vec <- apply(totalcost.mat, MARGIN=1, sum, na.rm=T)
    totcost.med[n, m] <- median(totcost.vec, na.rm=T)
    totcost.lo[n, m] <- quantile(totcost.vec, probs=0.025, na.rm=T)
    totcost.up[n, m] <- quantile(totcost.vec, probs=0.975, na.rm=T)
    
    print("##############################")
    print(paste("init harvest proportion = ", harv.prop.init[n], sep=""))
    print("##############################")
    
  } # end n loop (initial harvest rate)
  
  print("##############################")
  print(paste("maint harvest proportion = ", harv.prop.maint[m], sep=""))
  print("##############################")
  
} # end m loop (maintenance harvest rate)

## plot 3D surfaces
f1 <- list(
  family = "Avenir Light",
  size = 26,
  color = "black"
)
f2 <- list(
  family = "Avenir Light",
  size = 18,
  color = "black"
)
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)

# total cost (median)
par(mar=c(5,5,2,8))
costcontmed3d <- plot_ly(z = ~totcost.med, autocontour=T, type="contour", line = list(smoothing = 0.90), contours = list(showlabels = TRUE, labelfont=list(
  size=18, family="Avenir Light", face="bold", color="white"))) %>%
  colorbar(title = "tot $", titlefont=f2, tickfont=f2) %>%
  layout(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)))
costcontmed3d

cost3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~totcost.med) %>%
  add_surface(z = ~totcost.lo, opacity = 0.55) %>%
  add_surface(z = ~totcost.up, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="maintenance cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,0.5,0.1)), tickvals=seq(0,8,2)),
    yaxis = list(title="initial cull", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.5,0.9,0.1)), tickvals=seq(0,8,2)),
    zaxis = list(title="tot $", tickfont=f3, titlefont=f1)))
cost3d



#############################################################################################################################
###### TNR 
#############################################################################################################################
### SAME METHODS AS ABOVE, ALTERED FERTILITY INSTEAD OF SURVIVAL
# Trap Neuter Release #
# Kathryn Venning 
# Flinders University

#clear#
rm(list=ls())

## functions##

# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
#avaiable with this file on GitHub


#Matrix#
age.max = 7

## create vectors 
#fertility
TNR <- seq(.1,.9,.1)

f.vec <- c(0, 0.745, 2.52, 2.52, 2.52, 2.52, 1.98) ## KI cat birth rates matrix, data for female offsping produced each year. Data from Budke, C & Slater, M (2009)

# fertility errors based on Budke & Slater
juv.m.sd <- mean(c(((0.745-0.352)/2),((1.58-0.745)/2))*.5)
A.m.sd <- mean(c(((2.52-1.98)/2),((3.78-2.52)/2))*.5)
m.sd.vec <- c(0,juv.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd)

#survival 
s.vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7, 0.55) ##KI cat survival

# survival errors based on Budke & Slater
y1.2.S.sd <- mean(c(((0.46-0.27)/2),((0.73-0.46)/2)))
A.S.sd <- mean(c(((0.7-0.55)/2),((0.78-0.7)/2)))
s.sd.vec <- c(y1.2.S.sd,y1.2.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd)

# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max)
diag(popmat[2:age.max,]) <- s.vec[-age.max]
popmat[age.max,age.max] <- s.vec[age.max]
popmat[1,] <- f.vec
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 5000 # founding population size
init.vec <- stable.stage.dist(popmat) * pop.found


#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2018 # update if more data available post-2010
#************************
yr.end <- 2030 # set projection end date
#************************
t <- (yr.end - yr.now)

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig
yr.vec <- seq(yr.now,yr.end)

## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec

## set up projection loop
for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
yrs <- seq(yr.now, yr.end, 1)
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

# compensatory density feedback
K.max <- 2*pop.found
K.min <- 1
K.vec <- c(1,pop.found/2,pop.found,0.75*K.max,K.max)
red.thresh <- 0.75
red.vec <- c(1,0.99,0.95,0.875,0.75)
plot(K.vec,red.vec,pch=19,type="b")
Kred.dat <- data.frame(K.vec,red.vec)

# logistic power function a/(1+(x/b)^c)
param.init <- c(1, 15000, 2.5)
fit.lp <- nls(red.vec ~ a/(1+(K.vec/b)^c), 
              data = Kred.dat,
              algorithm = "port",
              start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
              trace = TRUE,      
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
fit.lp.summ <- summary(fit.lp)
plot(K.vec,red.vec,pch=19,xlab="N",ylab="reduction factor")
K.vec.cont <- seq(1,2*pop.found,1)
pred.lp.fx <- coef(fit.lp)[1]/(1+(K.vec.cont/coef(fit.lp)[2])^coef(fit.lp)[3])
lines(K.vec.cont,pred.lp.fx,lty=2,lwd=3,col="red")

a.lp <- coef(fit.lp)[1]
b.lp <- coef(fit.lp)[2]
c.lp <- coef(fit.lp)[3]


###################
# TNR
##################
# stochatic projection with density feedback
iter <- 1000
itdiv <- iter/10

TNR <- seq(.01,.9,.01)
q.ext <- 20
q.ext.vec <-  min.med.n <- min.lo.n <- min.up.n <- rep(0,length(TNR))

for (s in 1:length(TNR)) {
  
  #storage matrix
  n.sums.mat <- matrix(0, nrow = iter, ncol = (t+1))
  
  for (e in 1:iter){
    popmat <- popmat.orig
    
    n.mat <- matrix(0, nrow = age.max, ncol = (t+1))
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
      popmat[1,] <- fert.stoch  # add new stochastically resampled fertilities
      diag(popmat[2:age.max,]) <- s.stoch[-age.max]*pred.red # add new stochastically resampled survivals
      popmat[age.max,age.max] <- s.stoch[age.max]*pred.red # add new stochastically resampled survivals
      
      #fertility reduction 
      popmat[1,] <- popmat[1,]*TNR[s]
      
      # project
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      if (length(which(n.mat[,i+1] < 0)) > 0) {
        n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
      }
      
    } #end i loop
    
    n.sums.mat[e,] <- as.vector(colSums(n.mat))
    
    if (e %% itdiv==0) print(e) 
    
  } #end e loop
  
  min.pop.vec <- apply(n.sums.mat, MARGIN=1, min, na.rm=T)
  min.med.n[s] <- median(min.pop.vec, na.rm=T)
  min.lo.n[s] <- quantile(min.pop.vec, probs=0.025, na.rm=T)
  min.up.n[s] <- quantile(min.pop.vec, probs=0.975, na.rm=T)
  
  q.ext.vec[s] <- (sum(ifelse(round(min.pop.vec, 0) < q.ext, 1, 0)) / iter)
  q.ext.vec
  
  n.md <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) # minimum over all iterations
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yrs,n.md,type="l",xlab="year", ylab="minimum N", lwd=2, ylim=c(0.95*min(n.lo),1.05*max(n.up)))
  lines(yrs,n.lo,lty=2,col="red",lwd=1.5)
  lines(yrs,n.up,lty=2,col="red",lwd=1.5)
  
  print("##############")
  print(paste("TNR = ", TNR[s], sep=""))
  print("##############")
  
} #end s loop

plot(1-TNR, q.ext.vec, type="l", pch=19, xlab="proportion spayed each year", ylab="extinction probability")

plot(1-TNR, min.med.n, type="l", pch=19, xlab="proportion spayed each year", ylab="min N", ylim=c(min(min.lo.n),max(min.up.n)))
lines(1-TNR, min.lo.n, col="red", lty=2)
lines(1-TNR, min.up.n, col="red", lty=2)







