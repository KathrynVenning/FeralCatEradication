# Kathryn Venning, Corey Bradshaw, Frédérik Saltré
# Global Ecology, Flinders University — globalecologyflinders.com
# feral cat reduction on Kangaroo Island
# requires library - Plotly
### update 07/02/2021
## update includes: first year fertility, final year survival, predator reduction feedback, removed quasi extinction, previous version 'OFFICIAL cat eradication models GitHub'

## remove everything
rm(list = ls())

# libraries
library(plotly)
source("matrixOperators.r")
options(scipen = 1000)

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## source/matrix operators
#source("matrixOperators.r")

# create Leslie matrix
age.max = 7

## create vectors 
#fertility 
m.vec <- c((0.745/3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98) ## KI cat birth rates matrix, data for female offsping produced each year. Data from Budke, C & Slater, M (2009)

# fertility errors based on Budke & Slater
juv.m.sd <- mean(c(((0.745/3-0.352/3)/2),((1.58/3-0.745/3)/2))) #mean and standard deviations, juvenile fertility
fy.m.sd <- mean(c(((0.745-0.352)/2),((1.58-0.745)/2))) #mean and standard deviations, juvenile fertility
A.m.sd <- mean(c(((2.52-1.98)/2),((3.78-2.52)/2))) #mean and standard deviations, adult fertility
m.sd.vec <- c(0.18*m.vec[1],0.18*m.vec[2],A.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd) #mean and standard deviations vector, juvenile and adult fertility 

#survival
s.vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7) ##KI cat survival # probability of surviving from one year to the next. e.g surviving fourth year of life

# survival errors based on Budke & Slater
y1.2.S.sd <- mean(c(((0.46-0.27)/2),((0.73-0.46)/2))) #mean and standard deviations, juvenile survival
A.S.sd <- mean(c(((0.7-0.55)/2),((0.78-0.7)/2))) #mean and standard deviations, adult survival
s.sd.vec <- c(y1.2.S.sd,y1.2.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd) #mean and standard deviations vector, juvenile and adult survival

# create matrix
popmat <- matrix(data = 0, nrow=age.max, ncol=age.max)
diag(popmat[2:age.max,]) <- s.vec
popmat[age.max,age.max] <- 0
popmat[1,] <- m.vec
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat, age.max) # reproductive value
gen.l <- G.val(popmat, age.max) # mean generation length

## initial population vector
pop.found <- 1629 # +/- 661 founding population size Hohnen et al 2020 
ssd <- stable.stage.dist(popmat)
init.vec <- ssd * pop.found #initial population vector

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
red.vec <- c(1,0.965,0.89,0.79,0.71)
jpeg("k_vec.jpg")
plot(K.vec,red.vec,pch=19,type="b")
dev.off()
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
jpeg("reduction_factor.jpg")
plot(K.vec,red.vec,pch=19,xlab="N",ylab="reduction factor")
dev.off()
K.vec.cont <- seq(1,2*pop.found,1)
pred.lp.fx <- coef(fit.lp)[1]/(1+(K.vec.cont/coef(fit.lp)[2])^coef(fit.lp)[3])
lines(K.vec.cont,pred.lp.fx,lty=2,lwd=3,col="red")

a.lp <- coef(fit.lp)[1]
b.lp <- coef(fit.lp)[2]
c.lp <- coef(fit.lp)[3]

print(a.lp)
print(b.lp)
print(c.lp)

## compensatory density-feedback deterministic model
## set population storage matrices
n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
n.mat[,1] <- init.vec
popmat <- popmat.orig

## set up projection loop
for (i in 1:t) {
  totN.i <- sum(n.mat[,i])
  pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
  diag(popmat[2:age.max,]) <- s.vec*pred.red
  popmat[age.max,age.max] <- 0
  n.mat[,i+1] <- popmat %*% n.mat[,i]
}

n.pred <- colSums(n.mat)
jpeg("something_with_Carry_capacity.jpg")
plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N",ylim=c(0,1.05*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
abline(h=K.max, lty=2, col="red") #carry capacity
legend(yrs[2], n.pred[6], legend=c("N", "Carry capacity"),
       col=c("black", "red"), lty=1:2, cex=0.8)
dev.off()