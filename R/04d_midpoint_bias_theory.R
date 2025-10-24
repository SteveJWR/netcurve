## bias as a function of the midpoint position.

source("R/00_functions.R")


library(ggpubr)
library(lolaR)
save.plot = T
fig.height = 2600 # 1800
fig.width = 3200 # 3000
fig.res = 350

# TODO: add a common plot theme.
scale.plot = 1
max.dist.radius = 1.2
dyz = 2
location.scale = dyz/2
n.grid = 200 # TODO: return to 200
kappa.set <- c(-2,-1,-0.5,0,0.5,1)
# kappa.set <- c(-1,-0.8, -0.5, -0.4,-0.2,-0.1)
# kappa.set <- c(1,0.8, 0.5, 0.4,0.2,0.1)

# TODO: double check -0.5
plots.list <- list()
datasets <- list()
reference.widths <- c()
distances <- list()

# Magnitude of the bias capped.
bias.plot.cap = 5

# dimension of the space
p = 2

for(kappa.idx in seq_along(kappa.set)){

  kappa.true <- kappa.set[kappa.idx]

  if(kappa.true == 0){
    Z.init <- c(c(1,rep(0,p-1)),
                c(-1,rep(0,p-1)),
                c(0,sqrt(3),rep(0,p-2)))

    Z.init <- location.scale*Z.init
    Z.init <- matrix(Z.init, ncol = p, byrow = T)
  } else if(kappa.true > 0){

    pos1 = cos(2*location.scale*sqrt(kappa.true))/cos(location.scale*sqrt(kappa.true))
    pos2 = sqrt(1 - pos1**2)
    Z.init <- c(c(cos(location.scale*sqrt(kappa.true)),sin(location.scale*sqrt(kappa.true)),rep(0,p-1)),
                c(cos(location.scale*sqrt(kappa.true)),-sin(location.scale*sqrt(kappa.true)),rep(0,p-1)),
                c(pos1,0,pos2,rep(0,p-2)))
    Z.init <- matrix(Z.init, ncol = p + 1, byrow = T)

  } else if(kappa.true < 0){
    #
    pos1 = cosh(2*location.scale*sqrt(-kappa.true))/cosh(location.scale*sqrt(-kappa.true))
    pos2 = sqrt(1 + pos1**2)
    Z.init <- c(c(cosh(location.scale*sqrt(-kappa.true)),sinh(location.scale*sqrt(-kappa.true)),rep(0,p-1)),
                c(cosh(location.scale*sqrt(-kappa.true)),-sinh(location.scale*sqrt(-kappa.true)),rep(0,p-1)),
                c(pos1,0,pos2,rep(0,p-2)))
    Z.init <- matrix(Z.init, ncol = p + 1, byrow = T)
  }
  Z.init.euclidean <- c(c(1,rep(0,p-1)),
                        c(-1,rep(0,p-1)),
                        c(0,sqrt(3),rep(0,p-2)))
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

  for(j in 1:J){
    cat(paste("Locations:", j,"/",J), end = '\r')
    # Euclidean projection
    Z.euclidean.row <- c(dat[j,1],dat[j,2])
    # Z.zero.row <- c(0,0)
    Z.euclidean <- rbind(Z.init.euclidean, Z.euclidean.row)#, Z.zero.row)
    D.euclidean <- pos_to_dist(Z.euclidean,0)
    #dx0 = D.euclidean[4,5]

    if(kappa.true == 0){
      dxy = D.euclidean[3,1]
      dxz = D.euclidean[3,2]
      dyz = D.euclidean[2,1]
      dxm = D.euclidean[3,4]
    } else if(kappa.true > 0){

      if(sum(Z.euclidean.row**2) == 0){
        plotting.scale = 0
      } else {
        #plotting.scale = sqrt((1 - cos(dx0*sqrt(kappa.true)))/(sum(Z.euclidean.row**2)))
        z0 = cos(d.sur.mid*sqrt(kappa.true))
        plotting.scale = sqrt(1 - z0**2)/d.sur.mid
      }

      #z0 = sqrt(1 - plotting.scale^2*(sum(Z.euclidean.row**2)))
      Z.embedding.surrogate.mid = c(z0, plotting.scale*Z.euclidean.row)

      Z.embedding = rbind(Z.init, Z.embedding.surrogate.mid)
      D.embedding = pos_to_dist(Z.embedding, kappa.true)
      dxy = D.embedding[3,1]
      dxz = D.embedding[3,2]
      dyz = D.embedding[2,1]
      dxm = D.embedding[3,4]


    } else if(kappa.true < 0){
      if(sum(Z.euclidean.row**2) == 0){
        plotting.scale = 0

      } else {
        d.sur.mid = sqrt(sum(Z.euclidean.row**2))
        z0 = cosh(d.sur.mid*sqrt(-kappa.true))
        plotting.scale = sqrt(z0**2 - 1)/d.sur.mid

      }

      #z0 = sqrt(1 + plotting.scale^2*(sum(Z.euclidean.row**2)))
      Z.embedding.surrogate.mid = c(z0, plotting.scale*Z.euclidean.row)

      # Z.origin = c(1,0,0)
      # Z.block = rbind(Z.embedding.surrogate.mid,
      #                 Z.origin)
      # D.block = pos_to_dist(Z.block, kappa.true)
      # print(D.block)

      Z.embedding = rbind(Z.init, Z.embedding.surrogate.mid)
      D.embedding = pos_to_dist(Z.embedding, kappa.true)
      dxy = D.embedding[3,1]
      dxz = D.embedding[3,2]
      dyz = D.embedding[2,1]
      dxm = D.embedding[3,4]
      # TODO: Check for each of the embedding plots.
    }

    min.dm = MidDist(-bias.plot.cap + kappa.true, dxy, dxz, dyz)
    if(min.dm > dxm){
      kappa.est = -bias.plot.cap + kappa.true
    } else {
      # start.time <- Sys.time()
      kappa.est = EstimateKappa(dxy, dxz, dyz, dxm,
                                max.iter = 4,
                                ee.thresh = 0.001,
                                abs.ee.thresh = 10^(-6),
                                kappa.thresh = 1e-02,
                                min.curvature = -5000,
                                init.gridsize = 1000,
                                gridsearch = T)
      # end.time <- Sys.time()
      # time.taken <- end.time - start.time
      # print(time.taken)
      # print(kappa.est)
    }

    kappa.last = kappa.est
    bias = kappa.est - kappa.true
    bias[is.na(bias)] <- -bias.plot.cap
    # capping the sd for plotting purposes
    is_large = abs(bias) > bias.plot.cap
    dat[j,3] <- bias*(!is_large) + sign(bias)*bias.plot.cap*(is_large)
  }
  reference.widths[kappa.idx] <- dyz/2
  dat <- data.frame(dat)
  colnames(dat) <- c("x","y","bias")

  datasets[[kappa.idx]] <- dat
  b = round(seq(-bias.plot.cap,bias.plot.cap, length.out = 6))
  plot.tmp <- ggplot(datasets[[kappa.idx]], aes(x, y, fill=bias)) +
    geom_raster() + theme_classic() +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0)+
    guides(fill=guide_legend(reverse=TRUE))  +
    geom_point(aes(x=0,y=sqrt(3)),colour="black") +
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") +
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") +
    ggtitle(paste("Bias Heatmap, Curvature: ", round(kappa.true,1)))

  plot.tmp
  plots.list[[kappa.idx]] <- plot.tmp

}


save.plot.data <- T
if(save.plot.data){
  save(datasets, file = "plots/theoretical_bias_plot_datasets.rda")
}


load.plot.data <- F
if(load.plot.data){
  datasets <- load("plots/theoretical_bias_plot_datasets.rda")
}

dev.off()
ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
          plots.list[[4]], plots.list[[5]], plots.list[[6]],  common.legend = TRUE)



if(save.plot){
  ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
            plots.list[[4]], plots.list[[5]], plots.list[[6]],  common.legend = TRUE) %>%
    ggexport(filename = "plots/theoretical_bias.png", width = fig.width, height = fig.height, res = fig.res)
}






