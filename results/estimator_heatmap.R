

# chose the true curvature. 
set.seed(5)
kappa.true = 1
library(ggpubr)
source("00_functions.R")

max.dist.radius = 2.5
dyz = 2
location.scale = dyz/2
ell = 8

n.sims = 1
n.grid = 60 # granularity of the grid spacing
sim.parallel = F
bias.correction.plot = F

#reverse.bias.correction.plot = T
censoring.range = Inf # knowledge of the true curvature being somewhere in that range. 
B.boot <- 2

y = 1
z = 2
m = 3
x = 4


bias.scale = 1/2
sd.scale = 4

breaks = sd.scale*c(0.0,0.25,0.75,1.5,2.5)
bias.breaks = bias.scale*c(-9,-3,-1,0,1,3,9)
bias.limits = bias.scale*c(-9,9)
limits = sd.scale*c(0,2.5)


# fixed parameters
p = 2
kappa.set <- c(-2,-1,0,0.5,1)
bias.plots.list <- list()
sd.plots.list <- list()
rmse.plots.list <- list()
med.plots.list <- list()

cor.bias.plots.list <- list()
cor.sd.plots.list <- list()
cor.rmse.plots.list <- list()
cor.med.plots.list <- list()


reference.widths <- c()
datasets <- list()
for(kappa.idx in seq_along(kappa.set)){
  kappa.true <- kappa.set[kappa.idx]
  epsilon <- 10**(-6)
  # ideally we would like a large derivative in order to have a small
  # variance. 
  # we use this to make a heatmap of the magnitude of the derivative
  # showing the regions of where we would ideally like to take points 
  
  #y, z, m , x
  
  if(kappa.true == 0){
    Z.init <- c(c(1,rep(0,p-1)),
                c(-1,rep(0,p-1)),
                c(0,rep(0,p-1)))
    Z.init <- location.scale*Z.init
    Z.init <- matrix(Z.init, ncol = p, byrow = T)
  } else if(kappa.true > 0){
    #
    # location scale might need to be small in the spherical case
    Z.init <- c(c(cos(location.scale*sqrt(kappa.true)),sin(location.scale*sqrt(kappa.true)),rep(0,p-1)),
                c(cos(location.scale*sqrt(kappa.true)),-sin(location.scale*sqrt(kappa.true)),rep(0,p-1)),
                c(1,rep(0,p)))
    Z.init <- matrix(Z.init, ncol = p + 1, byrow = T)
    
  } else if(kappa.true < 0){
    #
    Z.init <- c(c(cosh(location.scale*sqrt(-kappa.true)),sinh(location.scale*sqrt(-kappa.true)),rep(0,p-1)),
                c(cosh(location.scale*sqrt(-kappa.true)),-sinh(location.scale*sqrt(-kappa.true)),rep(0,p-1)),
                c(1,rep(0,p)))
    Z.init <- matrix(Z.init, ncol = p + 1, byrow = T)
  }
  
  D.init <- pos_to_dist(Z.init, kappa.true)
  #print(D.init)
  
  
  if(kappa.true > 0 & pi/(sqrt(abs(kappa.true))) < max.dist.radius ){
    # turns into uniform over the whole sphere if it exceeds the range of the latent space. 
    latent.grid.radius = pi/(sqrt(abs(kappa.true))) - 10**(-6) # simulation radius of data 
  } else {
    latent.grid.radius = max.dist.radius
  }
  
  x.grid <- seq(-latent.grid.radius,latent.grid.radius, length.out = n.grid)
  
  dat = expand.grid(x.grid,x.grid)
  
  filtered.plot.idx = dat[,1]**2 + dat[,2]**2 <= latent.grid.radius**2
  
  
  dat = dat[filtered.plot.idx,]
  J <- nrow(dat)
  dat[,3] <- rep(NA,J) # bias 
  dat[,4] <- rep(NA,J) # sd
  dat[,5] <- rep(NA,J) # rmse
  dat[,6] <- rep(NA,J) # median
  if(bias.correction.plot){
    dat[,7] <- rep(NA,J) # bias 
    dat[,8] <- rep(NA,J) # sd
    dat[,9] <- rep(NA,J) # rmse
    dat[,10] <- rep(NA,J) # median
  }
  
  for(j in 1:J){
    cat(paste("Locations:", j,"/",J), end = '\r')
    if(kappa.true > 0 ){
      d0 = sqrt(dat[j,1]**2 + dat[j,2]**2)
      z0 = cos(sqrt(kappa.true)*d0)
      scl <- sqrt((1 - z0^2)/(dat[j,1]^2 + dat[j,2]^2))
      Z.new <- c(z0, scl*dat[j,1],scl*dat[j,2])
      
    } else if(kappa.true < 0) { 
      d0 = sqrt(dat[j,1]**2 + dat[j,2]**2)
      z0 = cosh(sqrt(-kappa.true)*d0)
      scl <- sqrt((z0^2 - 1)/(dat[j,1]^2 + dat[j,2]^2))
      Z.new <- c(z0, scl*dat[j,1],scl*dat[j,2])
      
    } else if(kappa.true == 0){
      Z.new <- c(dat[j,1],dat[j,2])
      
    }
    #print(Z.new)
    Z.block <- rbind(Z.init, Z.new)
    D.small <- pos_to_dist(Z.block,kappa.true)
    #print(D.small)
    Z <- Z.block[rep(1:nrow(Z.block), each = ell), ]
    D.block <- pos_to_dist(Z.block,kappa.true)
    D <- pos_to_dist(Z,kappa.true)
    
    # dxy = D.small[4,1]
    # dxz = D.small[4,2]
    # dyz = D.small[2,1]
    # dxm = D.small[4,3]
    # 
    # d.vec <- c(dxy,dxz,dyz,dxm)
    # print(d.vec)
    # x.seq = seq(-300,2.3,0.02)
    # y = c()
    # for(s in x.seq){
    #   y <- c(y,g_ee(s,dxy, dxz,dyz,dxm))
    # }
    # max.kappa = (pi/(2*(max(d.vec))))**2
    # plot(x.seq,y)
    # abline(h = 0)
    # abline(v = max.kappa)
    
    rand.eff <- rep(0,4*ell)
    kappa.hat.set <- rep(0,n.sims)
    bias.corrected.kappa.hat.set <- kappa.hat.set
    #reverse.bias.corrected.kappa.hat.set <- kappa.hat.set
    clique.set = matrix(1:(4*ell), ncol = ell, byrow = T)
    K <- nrow(clique.set)
    nus.hat.mat <- matrix(nrow = K, ncol = ell)
    
    for(cl in 1:K){
      clique.idx <- clique.set[cl,]
      nus.hat <- rep(0,ell)
      nus.hat.mat[cl,] <- nus.hat
    }
    
    if(sim.parallel){
      results.list <- mclapply(1:n.sims, function(sim){
        A <- sim_ls_network(rand.eff,D)
        D.est = D_estimate(A,clique.set,fixed.effects = nus.hat.mat)
        D.hat = D.est$estimates
        
        dxy.hat = D.hat[4,1]
        dxz.hat = D.hat[4,2]
        dyz.hat = D.hat[2,1]
        dxm.hat = D.hat[4,3]
        if(!any(is.infinite(D.hat))){
          kappa.hat <- tryCatch({
            kappa.tmp <- estimate_kappa(dxy.hat,dxz.hat,dyz.hat,dxm.hat)
            if(abs(kappa.tmp) > censoring.range){
              kappa.tmp = NA 
            }
            kappa.tmp
          },
          error = function(e) {
            return(NA)
          })
          if(bias.correction.plot){
            bias.estimate <- parametric_bootstrap_dev_estimate(A,nus.hat.mat,clique.set,B.boot,y,z,m,x,parallel = F)
            bias.corrected.kappa.hat <- kappa.hat - bias.estimate
          } 
          # d.hat.vec <- c(dxy.hat, dxz.hat,dyz.hat,dxm.hat)
          # print(d.hat.vec)
          # max.kappa = (pi/(2*(max(d.hat.vec))))**2 
          # x.seq = seq(-900,-850,0.02)
          #  
          # y = c()
          # for(s in x.seq){
          #   y <- c(y,g_ee(s,dxy.hat, dxz.hat,dyz.hat,dxm.hat))
          # }
          # 
          # plot(x.seq,y)
          # abline(h = 0)
          # abline(v = max.kappa)
          # print(kappa.hat)
        } else {
          kappa.hat <- NA
          bias.corrected.kappa.hat <- NA
        }
        
        out.set <- list(kappa.hat,bias.corrected.kappa.hat)
        return(out.set)
      })
      
      
      
      for(sim in 1:n.sims){
        kappa.hat.set[sim] <- results.list[[sim]][[1]]
        
        if(bias.correction.plot){
          bias.corrected.kappa.hat.set[sim] <- results.list[[sim]][[2]]
        }
      }
    } else {
      for(sim in 1:n.sims){
        A <- sim_ls_network(rand.eff,D)
        
        
        D.est = D_estimate(A,clique.set,fixed.effects = nus.hat.mat)
        D.hat = D.est$estimates
        
        dxy.hat = D.hat[4,1]
        dxz.hat = D.hat[4,2]
        dyz.hat = D.hat[2,1]
        dxm.hat = D.hat[4,3]
        if(!any(is.infinite(D.hat))){
          kappa.hat <- tryCatch({
            kappa.tmp <- estimate_kappa(dxy.hat, dxz.hat,dyz.hat,dxm.hat)
            if(abs(kappa.tmp) > censoring.range){
              kappa.tmp = NA 
            }
            kappa.tmp
          },
          error = function(e) {
            return(NA)
          })
          if(bias.correction.plot){

            bias.estimate <- parametric_bootstrap_dev_estimate(A,nus.hat.mat,clique.set,B.boot,y,z,m,x,parallel = F)
            bias.corrected.kappa.hat <- kappa.hat - bias.estimate
          } 
        } else {
          kappa.hat <- NA
          if(bias.correction.plot){
            bias.corrected.kappa.hat <- NA
          }
        }
        
        
        kappa.hat.set[sim] <- kappa.hat
        if(bias.correction.plot){
          bias.corrected.kappa.hat.set[sim] <- bias.corrected.kappa.hat
        }
      }
    }
    #plot(density(kappa.hat.set, na.rm = T), xlim = c(-20,5))
    #abline(v = kappa.true, lty = 3, col = "purple")
    dat[j,3] <- mean(kappa.hat.set,na.rm = T) - kappa.true # bias
    dat[j,4] <- sd(kappa.hat.set,na.rm = T)# sd 
    dat[j,5] <- sqrt(mean((kappa.hat.set - kappa.true)**2,na.rm = T))# rmse 
    dat[j,6] <- median(kappa.hat.set,na.rm = T) - kappa.true # bias
    if(bias.correction.plot){
      dat[j,7] <- mean(bias.corrected.kappa.hat.set,na.rm = T) - kappa.true # bias
      dat[j,8] <- sd(bias.corrected.kappa.hat.set,na.rm = T)# sd 
      dat[j,9] <- sqrt(mean((bias.corrected.kappa.hat.set - kappa.true)**2,na.rm = T))# rmse 
      dat[j,10] <- median(bias.corrected.kappa.hat.set,na.rm = T) - kappa.true # bias
    } 

  }
 

  reference.widths <- dyz/2
  dat <- data.frame(dat)
  
  
  if(bias.correction.plot){
    colnames(dat) <- c("x","y","bias","sd","rmse","median",
                       "cor.bias","cor.sd","cor.rmse","cor.median")  
  } else{
    colnames(dat) <- c("x","y","bias","sd","rmse","median")
  }
  
  
  datasets[[kappa.idx]] <- dat

  bias.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=bias)) + 
    geom_raster() + theme_classic() +
    scale_fill_gradientn(limits = bias.limits,
                         colours=c("navyblue","navyblue","#0000FFFF","#FFFFFFFF","#FF0000FF","darkred","darkred"),
                         breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    # scale_fill_gradientn(limits = c(0,0.8),
    #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
    #                      breaks=b, labels=format(b)) + 
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths,y=0),colour="black") + 
    geom_point(aes(x=-reference.widths,y=0),colour="black") + 
    ggtitle(paste("Bias Heatmap, Curvature: ", kappa.true))
  
  sd.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=sd)) + 
    geom_raster() + theme_classic() +
    # scale_fill_gradientn(limits = c(-1,1),
    #                      colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
    #                      breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    scale_fill_gradientn(limits = limits,
                         colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
                         breaks=breaks, labels=format(breaks)) +
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths,y=0),colour="black") + 
    geom_point(aes(x=-reference.widths,y=0),colour="black") + 
    ggtitle(paste("SD Heatmap, Curvature: ", kappa.true))
  
  rmse.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=rmse)) + 
    geom_raster() + theme_classic() +
    # scale_fill_gradientn(limits = c(-1,1),
    #                      colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
    #                      breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    scale_fill_gradientn(limits = limits,
                         colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
                         breaks=breaks, labels=format(breaks)) +
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths,y=0),colour="black") + 
    geom_point(aes(x=-reference.widths,y=0),colour="black") + 
    ggtitle(paste("RMSE Heatmap, Curvature: ", kappa.true))
  
  med.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=median)) + 
    geom_raster() + theme_classic() +
    scale_fill_gradientn(limits = bias.limits,
                         colours=c("navyblue","navyblue","#0000FFFF","#FFFFFFFF","#FF0000FF","darkred","darkred"),
                         breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    # scale_fill_gradientn(limits = c(0,0.8),
    #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
    #                      breaks=b, labels=format(b)) + 
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths,y=0),colour="black") + 
    geom_point(aes(x=-reference.widths,y=0),colour="black") + 
    ggtitle(paste("Median Heatmap, Curvature: ", kappa.true))
  
  
  
  if(bias.correction.plot){
    cor.bias.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=cor.bias)) + 
      geom_raster() + theme_classic() +
      scale_fill_gradientn(limits = bias.limits,
                           colours=c("navyblue","navyblue","#0000FFFF","#FFFFFFFF","#FF0000FF","darkred","darkred"),
                           breaks=bias.breaks, labels=format(bias.breaks)) +
      #scale_fill_viridis_c(option = "magma") + 
      # scale_fill_gradientn(limits = c(0,0.8),
      #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
      #                      breaks=b, labels=format(b)) + 
      guides(fill=guide_legend(reverse=TRUE))  + 
      geom_point(aes(x=0,y=0),colour="black") + 
      geom_point(aes(x=reference.widths,y=0),colour="black") + 
      geom_point(aes(x=-reference.widths,y=0),colour="black") + 
      ggtitle(paste("Bias Heatmap, Curvature: ", kappa.true))
    
    cor.sd.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=cor.sd)) + 
      geom_raster() + theme_classic() +
      # scale_fill_gradientn(limits = c(-1,1),
      #                      colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
      #                      breaks=bias.breaks, labels=format(bias.breaks)) +
      #scale_fill_viridis_c(option = "magma") + 
      scale_fill_gradientn(limits = limits,
                           colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
                           breaks=breaks, labels=format(breaks)) +
      guides(fill=guide_legend(reverse=TRUE))  + 
      geom_point(aes(x=0,y=0),colour="black") + 
      geom_point(aes(x=reference.widths,y=0),colour="black") + 
      geom_point(aes(x=-reference.widths,y=0),colour="black") + 
      ggtitle(paste("SD Heatmap, Curvature: ", kappa.true))
    
    cor.rmse.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=cor.rmse)) + 
      geom_raster() + theme_classic() +
      # scale_fill_gradientn(limits = c(-1,1),
      #                      colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
      #                      breaks=bias.breaks, labels=format(bias.breaks)) +
      #scale_fill_viridis_c(option = "magma") + 
      scale_fill_gradientn(limits = limits,
                           colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
                           breaks=breaks, labels=format(breaks)) +
      guides(fill=guide_legend(reverse=TRUE))  + 
      geom_point(aes(x=0,y=0),colour="black") + 
      geom_point(aes(x=reference.widths,y=0),colour="black") + 
      geom_point(aes(x=-reference.widths,y=0),colour="black") + 
      ggtitle(paste("RMSE Heatmap, Curvature: ", kappa.true))
    
    cor.med.plot.tmp <- ggplot(data = datasets[[kappa.idx]], aes(x = x,y = y, fill=cor.median)) + 
      geom_raster() + theme_classic() +
      scale_fill_gradientn(limits = bias.limits,
                           colours=c("navyblue","navyblue","#0000FFFF","#FFFFFFFF","#FF0000FF","darkred","darkred"),
                           breaks=bias.breaks, labels=format(bias.breaks)) +
      #scale_fill_viridis_c(option = "magma") + 
      # scale_fill_gradientn(limits = c(0,0.8),
      #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
      #                      breaks=b, labels=format(b)) + 
      guides(fill=guide_legend(reverse=TRUE))  + 
      geom_point(aes(x=0,y=0),colour="black") + 
      geom_point(aes(x=reference.widths,y=0),colour="black") + 
      geom_point(aes(x=-reference.widths,y=0),colour="black") + 
      ggtitle(paste("Median Heatmap, Curvature: ", kappa.true))
  }
  bias.plots.list[[kappa.idx]] <- bias.plot.tmp
  sd.plots.list[[kappa.idx]] <- sd.plot.tmp
  rmse.plots.list[[kappa.idx]] <- rmse.plot.tmp
  med.plots.list[[kappa.idx]] <- med.plot.tmp
  
  if(bias.correction.plot){
    cor.bias.plots.list[[kappa.idx]] <- cor.bias.plot.tmp
    cor.sd.plots.list[[kappa.idx]] <- cor.sd.plot.tmp
    cor.rmse.plots.list[[kappa.idx]] <- cor.rmse.plot.tmp
    cor.med.plots.list[[kappa.idx]] <- cor.med.plot.tmp
  }
  
  print(paste("Curvature", kappa.idx,"of", length(kappa.set)))
}



text <- paste("Bias Map, Clique Size:",ell)

# Create a text grob
bias.tgrob <- text_grob(text,size = 20)
bias.plot_text <- as_ggplot(bias.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

bias.plot <- ggarrange(bias.plot_text,bias.plots.list[[1]], bias.plots.list[[2]], bias.plots.list[[3]], 
          bias.plots.list[[4]], bias.plots.list[[5]],  common.legend = TRUE)


text <- paste("SD Map, Clique Size:",ell)

# Create a text grob
sd.tgrob <- text_grob(text,size = 20)
sd.plot_text <- as_ggplot(sd.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

sd.plot <- ggarrange(sd.plot_text,sd.plots.list[[1]], sd.plots.list[[2]], sd.plots.list[[3]], 
                       sd.plots.list[[4]], sd.plots.list[[5]],  common.legend = TRUE)




text <- paste("RMSE map, Clique Size:",ell)

# Create a text grob
rmse.tgrob <- text_grob(text,size = 20)
rmse.plot_text <- as_ggplot(rmse.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

rmse.plot <- ggarrange(rmse.plot_text,rmse.plots.list[[1]], rmse.plots.list[[2]], rmse.plots.list[[3]], 
                     rmse.plots.list[[4]], rmse.plots.list[[5]],  common.legend = TRUE)

text <- paste("Median map, Clique Size:",ell)

# Create a text grob
med.tgrob <- text_grob(text,size = 20)
med.plot_text <- as_ggplot(med.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

med.plot <- ggarrange(med.plot_text,med.plots.list[[1]], med.plots.list[[2]], med.plots.list[[3]], 
                      med.plots.list[[4]], med.plots.list[[5]],  common.legend = TRUE)







#############################
# Create a text grob
cor.bias.tgrob <- text_grob(text,size = 20)
cor.bias.plot_text <- as_ggplot(cor.bias.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

cor.bias.plot <- ggarrange(cor.bias.plot_text,cor.bias.plots.list[[1]], cor.bias.plots.list[[2]], cor.bias.plots.list[[3]], 
                           cor.bias.plots.list[[4]], cor.bias.plots.list[[5]],  common.legend = TRUE)


text <- paste("SD Map, Clique Size:",ell)

# Create a text grob
cor.sd.tgrob <- text_grob(text,size = 20)
cor.sd.plot_text <- as_ggplot(cor.sd.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

cor.sd.plot <- ggarrange(cor.sd.plot_text,cor.sd.plots.list[[1]], cor.sd.plots.list[[2]], cor.sd.plots.list[[3]], 
                         cor.sd.plots.list[[4]], cor.sd.plots.list[[5]],  common.legend = TRUE)




text <- paste("RMSE map, Clique Size:",ell)

# Create a text grob
cor.rmse.tgrob <- text_grob(text,size = 20)
cor.rmse.plot_text <- as_ggplot(cor.rmse.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

cor.rmse.plot <- ggarrange(cor.rmse.plot_text,rmse.plots.list[[1]], cor.rmse.plots.list[[2]], cor.rmse.plots.list[[3]], 
                           cor.rmse.plots.list[[4]], cor.rmse.plots.list[[5]],  common.legend = TRUE)


text <- paste("Median map, Clique Size:",ell)

# Create a text grob
cor.med.tgrob <- text_grob(text,size = 20)
cor.med.plot_text <- as_ggplot(cor.med.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

cor.med.plot <- ggarrange(cor.med.plot_text,med.plots.list[[1]], cor.med.plots.list[[2]], cor.med.plots.list[[3]], 
                          cor.med.plots.list[[4]], cor.med.plots.list[[5]],  common.legend = TRUE)







#############################
# Create a text grob


# bias.plot
# sd.plot
# rmse.plot

# cor.bias.plot
# cor.sd.plot
# cor.rmse.plot


# saveRDS(bias.plot, "results/bias_map.rds")
# saveRDS(sd.plot, "results/sd_map.rds")
# saveRDS(rmse.plot, "results/rmse_map.rds")
# saveRDS(med.plot, "results/median_map.rds")
# 
# 
# saveRDS(cor.bias.plot, "results/cor_bias_map.rds")
# saveRDS(cor.sd.plot, "results/cor_sd_map.rds")
# saveRDS(cor.rmse.plot, "results/cor_rmse_map.rds")
# saveRDS(cor.med.plot, "results/cor_median_map.rds")
# 
# 
# ggsave(paste0("results/bias_full_ell_",ell,".png"),bias.plot)
# ggsave(paste0("results/sd_full_ell_",ell,".png"),sd.plot)
# ggsave(paste0("results/rmse_full_ell_",ell,".png"),rmse.plot)
# ggsave(paste0("results/med_full_ell_",ell,".png"),med.plot)
# ggsave(paste0("results/cor_bias_full_ell_",ell,".png"),cor.bias.plot)
# ggsave(paste0("results/cor_sd_full_ell_",ell,".png"),cor.sd.plot)
# ggsave(paste0("results/cor_rmse_full_ell_",ell,".png"),cor.rmse.plot)
# ggsave(paste0("results/cor_med_full_ell_",ell,".png"),cor.med.plot)

# ggarrange(bias.plots.list[[1]], bias.plots.list[[2]], bias.plots.list[[3]], 
#           bias.plots.list[[4]], bias.plots.list[[5]], bias.plots.list[[6]],  common.legend = TRUE)

# 
# bias.plot
# cor.bias.plot
# 
# med.plot
# cor.med.plot
# 
# rmse.plot
# cor.rmse.plot