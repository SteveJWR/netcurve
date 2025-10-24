

# chose the true curvature.
source("R/00_functions.R")

library(ggpubr)
library(lolaR)
save.plot = T
fig.height = 2600 # 1800
fig.width = 3200 # 3000
fig.res = 350


scale.plot = 1
max.dist.radius = 2
dyz = 2
location.scale = dyz/2
n.grid = 200 # TODO: return to 200
kappa.set <- c(-2,-1,-0.5,0,0.5,1)
plots.list <- list()
datasets <- list()
reference.widths <- c()
distances <- list()

# Magnitude of the sd capped for plotting
sd.plot.cap = 50 # TODO: check whether 15 is reasonable

p = 2
for(kappa.idx in seq_along(kappa.set)){

  kappa.true <- kappa.set[kappa.idx]

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
  Z.init.euclidean <- c(c(1,rep(0,p-1)),
              c(-1,rep(0,p-1)),
              c(0,rep(0,p-1)))
  Z.init.euclidean <- location.scale*Z.init.euclidean
  Z.init.euclidean <- matrix(Z.init.euclidean, ncol = p, byrow = T)

  #plotting based on a projection on for the same radius and angle as the Euclidean Space.
  latent.grid.radius = max.dist.radius
  x.grid <- seq(-latent.grid.radius,latent.grid.radius, length.out = n.grid)

  dat = expand.grid(x.grid,x.grid)
  filtered.plot.idx = dat[,1]**2 + dat[,2]**2 <= latent.grid.radius**2
  dat = dat[filtered.plot.idx,]
  J <- nrow(dat)
  dat[,3] <- rep(NA,J)
  dat[,4] <- rep(NA,J)

  for(j in 1:J){
    cat(paste("Locations:", j,"/",J), end = '\r')

    Z.euclidean.row <- c(dat[j,1],dat[j,2])
    Z.euclidean <- rbind(Z.init.euclidean, Z.euclidean.row)
    D.euclidean <- pos_to_dist(Z.euclidean,0)


    dxy = D.euclidean[4,1]
    dxz = D.euclidean[4,2]
    dyz = D.euclidean[2,1]
    dxm = lolaR::MidDist(kappa.true, dxy, dxz, dyz)

    var.theory.out <- sum(g_grad_d(kappa.true, dxy, dxz, dyz, dxm)^2)/g_grad_kappa(kappa.true,
                                                                                   dxy,dxz,
                                                                                   dyz,dxm,
                                                                                   manual.grad = F)^2
    var.theory.out <- Re(var.theory.out)

    # capping the sd for plotting purposes
    dat[j,3] <- min(sqrt(var.theory.out), sd.plot.cap)
  }

  if(kappa.true == 0){
    ref.x = dyz/2
  } else if(kappa.true > 0){
    ref.x <- dyz/2
  } else if(kappa.true < 0){
    ref.x <- dyz/2
  }
  reference.widths[kappa.idx] <- ref.x
  dat <- data.frame(dat)
  colnames(dat) <- c("x","y","SD","manual_grad")

  datasets[[kappa.idx]] <- dat
  b = round(seq(0,sd.plot.cap, length.out = 6))
  plot.tmp <- ggplot(datasets[[kappa.idx]], aes(x, y, fill=SD)) +
    geom_raster() + theme_classic() +
    #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
    #scale_fill_viridis_c(option = "magma") +
    scale_fill_gradientn(limits = c(0,sd.plot.cap),
                         colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow","chartreuse1"),
                         breaks=b, labels=format(b)) +
    guides(fill=guide_legend(reverse=TRUE))  +
    geom_point(aes(x=0,y=0),colour="black") +
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") +
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") +
    ggtitle(paste("SD Heatmap, Curvature: ", round(kappa.true,1)))
  #plot.tmp
  plots.list[[kappa.idx]] <- plot.tmp

}



save.plot.data <- F
if(save.plot.data){
  save(datasets, file = "plots/theoretical_variance_plot_datasets.rda")
}


load.plot.data <- T
if(load.plot.data){
  datasets <- load("plots/theoretical_variance_plot_datasets.rda")
}


dev.off()
ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
          plots.list[[4]], plots.list[[5]], plots.list[[6]],  common.legend = TRUE)


if(save.plot){
  ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
            plots.list[[4]], plots.list[[5]], plots.list[[6]],  common.legend = TRUE) %>%
    ggexport(filename = "plots/theoretical_variance.png", width = fig.width, height = fig.height, res = fig.res)
}










