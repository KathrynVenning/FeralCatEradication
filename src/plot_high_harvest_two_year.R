harv.prop.init <- seq(0.5,0.9,0.05)
harv.prop.maint <- seq(0.1,0.5,0.05)

# storage
minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init)) #storage matrices

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
        
        popmat[1,] <- fert.stoch
        diag(popmat[2:age.max,]) <- s.stoch*pred.red
        #popmat[age.max,age.max] <- 0
        
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

twophase.med <- data.frame(pmin.med.mat)
colnames(twophase.med) <- harv.prop.maint
rownames(twophase.med) <- harv.prop.init

twophase.lo <- data.frame(pmin.lo.mat)
colnames(twophase.lo) <- harv.prop.maint
rownames(twophase.lo) <- harv.prop.init

twophase.up <- data.frame(pmin.up.mat)
colnames(twophase.up) <- harv.prop.maint
rownames(twophase.up) <- harv.prop.init

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