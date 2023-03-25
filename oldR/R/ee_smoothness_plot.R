

# chose the true curvature. 

kappa.true = -1


max.dist.radius = 1.5
dyz = 2
location.scale = dyz/2
n.grid = 60
kappa.set <- c(-1,-0.5,0,0.5,1)
plots.list <- list()
datasets <- list()
reference.widths <- c()
for(kappa.idx in seq_along(kappa.set)){

  kappa.true <- kappa.set[kappa.idx]
  epsilon <- 10**(-6)
  # ideally we would like a large derivative in order to have a small
  # variance. 
  # we use this to make a heatmap of the magnitude of the derivative
  # showing the regions of where we would ideally like to take points 
  
  #y, z, m , x
  p = 2
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
  
  
  if(kappa.true > 0 & pi/(2*sqrt(abs(kappa.true))) < max.dist.radius){
    # turns into uniform over the whole sphere if it exceeds the range of the latent space. 
    latent.grid.radius = sin(pi/(2*sqrt(kappa.true))) - 10**(-6) # simulation radius of data 
    
  } else if(kappa.true > 0){
    latent.grid.radius = sin(sqrt(kappa.true)*max.dist.radius)
  } else if(kappa.true < 0){
    latent.grid.radius = sinh(sqrt(-kappa.true)*max.dist.radius)
  } else {
    latent.grid.radius = max.dist.radius
  }
  
  x.grid <- seq(-latent.grid.radius,latent.grid.radius, length.out = n.grid)
  
  dat = expand.grid(x.grid,x.grid)
  filtered.plot.idx = dat[,1]**2 + dat[,2]**2 <= latent.grid.radius**2
  dat = dat[filtered.plot.idx,]
  J <- nrow(dat)
  dat[,3] <- rep(NA,J)
  dat[,4] <- rep(NA,J)
  
  for(j in 1:J){
    cat(paste("Locations:", j,"/",J), end = '\r')
    if(kappa.true > 0){
      Z.new <- c(1 - sqrt(dat[j,1]^2 + dat[j,2]^2), dat[j,1],dat[j,2])
      
    } else if(kappa.true < 0) { 
      Z.new <- c(1 + sqrt(dat[j,1]^2 + dat[j,2]^2), dat[j,1],dat[j,2])
      
    } else if(kappa.true == 0){
      Z.new <- c(dat[j,1],dat[j,2])
      
    }
    Z.block <- rbind(Z.init, Z.new)
    D.block <- pos_to_dist(Z.block,kappa.true)
    
    dxy.block = D.block[4,1]
    dxz.block = D.block[4,2]
    dyz.block = D.block[2,1]
    dxm.block = D.block[4,3]
    grad.out <- abs(g_grad_kappa(kappa.true,
                                 dxy.block,dxz.block,
                                 dyz.block,dxm.block))
    dat[j,3] <- grad.out
    manual.grad <- (1/epsilon)*(g_ee(kappa.true + epsilon,
                                     dxy.block,dxz.block,
                                     dyz.block,dxm.block) - g_ee(kappa.true,
                                                                 dxy.block,dxz.block,
                                                                 dyz.block,dxm.block))
    dat[j,4] <- manual.grad
  }
  
  if(kappa.true == 0){
    ref.x = dyz/2
  } else if(kappa.true > 0){
    ref.x <- Re(sin(sqrt(kappa.true*dxy/2 + 0i)))
  } else if(kappa.true < 0){
    ref.x <- Im(sin(sqrt(kappa.true*dxy/2 + 0i)))
  }
  reference.widths[kappa.idx] <- ref.x
  dat <- data.frame(dat)
  colnames(dat) <- c("x","y","grad","manual_grad")
  
  datasets[[kappa.idx]] <- dat
  b = c(0,0.2,0.4,0.6,0.8,1.2)
  plot.tmp <- ggplot(datasets[[kappa.idx]], aes(x, y, fill=grad)) + 
    geom_raster() + theme_classic() +
    #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) + 
    #scale_fill_viridis_c(option = "magma") + 
    scale_fill_gradientn(limits = c(0,1.2),
                         colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow","chartreuse1"),
                         breaks=b, labels=format(b)) + 
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="grey") + 
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="grey") + 
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="grey") + 
    ggtitle(paste("Gradient Heatmap, Curvature: ", kappa.true))
  #plot.tmp
  plots.list[[kappa.idx]] <- plot.tmp
  # ggplot(dat, aes(x, y, fill=manual_grad)) + 
  #   geom_raster() + theme_classic() +
  #   scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) + 
  #   guides(fill=guide_legend(reverse=TRUE))  
}


library(ggpubr)
ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]], 
          plots.list[[4]], plots.list[[5]],  common.legend = TRUE)
