

# chose the true curvature.
#TODO: Fix the theoretical variance here.
source("R/00_functions.R")
kappa.true = -1

scale.plot = 1
max.dist.radius = 2
dyz = 2
location.scale = dyz/2
n.grid = 300 # 200
kappa.set <- c(-2,-1,0,0.5,1)
plots.list <- list()
datasets <- list()
reference.widths <- c()
distances <- list()

sd.plot.cap = 15

for(kappa.idx in seq_along(kappa.set)){

  kappa.true <- kappa.set[kappa.idx]
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
  Z.init.euclidean <- c(c(1,rep(0,p-1)),
              c(-1,rep(0,p-1)),
              c(0,rep(0,p-1)))
  Z.init.euclidean <- location.scale*Z.init.euclidean
  Z.init.euclidean <- matrix(Z.init.euclidean, ncol = p, byrow = T)

  # if(kappa.true > 0 & pi/(2*sqrt(abs(kappa.true))) < max.dist.radius){
  #   # turns into uniform over the whole sphere if it exceeds the range of the latent space.
  #   latent.grid.radius = sin(pi/(2*sqrt(kappa.true))) - 10**(-6) # simulation radius of data
  #
  # } else if(kappa.true > 0){
  #   latent.grid.radius = sin(sqrt(kappa.true)*max.dist.radius)
  # } else if(kappa.true < 0){
  #   latent.grid.radius = sinh(sqrt(-kappa.true)*max.dist.radius)
  # } else {
  #   latent.grid.radius = max.dist.radius
  # }

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

    # if(kappa.true > 0){
    #
    #   M = matrix(c(Z.init[1,1], Z.init[1,2], Z.init[2,1], Z.init[2,2]), nrow = 2, ncol =2, byrow = 2)
    #   v = c(cos(sqrt(kappa.true)*D.euclidean[1,4]), cos(sqrt(kappa.true)*D.euclidean[2,4]))
    #   z.short = solve(M) %*% v
    #   Z.new <-c(z.short[1], z.short[2], sqrt(1 - z.short[1]^2 +  z.short[2]^2))
    #
    # } else if(kappa.true < 0) {
    #
    #   M = matrix(c(-Z.init[1,1], Z.init[1,2], -Z.init[2,1], Z.init[2,2]), nrow = 2, ncol =2, byrow = 2)
    #   v = c(cosh(sqrt(-kappa.true)*D.euclidean[1,4]), cosh(sqrt(-kappa.true)*D.euclidean[2,4]))
    #   z.short = solve(M) %*% v
    #   Z.new <-c(z.short[1], z.short[2], sqrt(1 + z.short[1]^2 +  z.short[2]^2))
    #
    # } else if(kappa.true == 0){
    #   Z.new <- c(dat[j,1],dat[j,2])
    #
    # }
    #
    # Z.block <- rbind(Z.init, Z.new)
    # D.block <- pos_to_dist(Z.block,kappa.true) #TODO: Fix this so that the corrseponding euclideean distance to the midpoint is matched


    # dxy.block = D.block[4,1]
    # dxz.block = D.block[4,2]
    # dyz.block = D.block[2,1]
    # dxm.block = D.block[4,3]

    dxy = D.euclidean[4,1] #1
    dxz = D.euclidean[4,2] #3
    dyz = D.euclidean[2,1] #2
    dxm = lolaR::MidDist(kappa.true, dxy, dxz, dyz)

    # y = c()
    # x = seq(-2,1,length.out = 100)
    # for(k in x){
    #   y <- c(y,lolaR::MidDist(k, dxy, dxz, dyz))
    # }
    # plot(x,y)

#
#     if(length(distances[[kappa.idx]]) > 0 ){
#       distances[[kappa.idx]] <- c(dxy.block - dxy.euclidean, dxz.block - dxz.euclidean,
#                                   dyz.block - dyz.euclidean)#, dxm.block - dxm.euclidean)
#     } else {
#       distances[[kappa.idx]] <- c(distances[[kappa.idx]],
#                                   dxy.block - dxy.euclidean, dxz.block - dxz.euclidean,
#                                   dyz.block - dyz.euclidean)#, dxm.block - dxm.euclidean)
#     }


    var.theory.out <- sum(g_grad_d(kappa.true, dxy, dxz, dyz, dxm)^2)/g_grad_kappa(kappa.true,
                                                                                   dxy,dxz,
                                                                                   dyz,dxm,
                                                                                   manual.grad = F)^2
    var.theory.out <- Re(var.theory.out)

    # capping the sd for plotting purposes
    dat[j,3] <- min(sqrt(sqrt(var.theory.out)), sd.plot.cap)
  }

  if(kappa.true == 0){
    ref.x = dyz/2
  } else if(kappa.true > 0){
    ref.x <- dyz/2 #Re(sin(sqrt(kappa.true*dyz/2 + 0i)))
  } else if(kappa.true < 0){
    ref.x <- dyz/2 #Im(sin(sqrt(kappa.true*dyz/2 + 0i)))
  }
  reference.widths[kappa.idx] <- ref.x
  dat <- data.frame(dat)
  colnames(dat) <- c("x","y","sqrt_sd","manual_grad")

  datasets[[kappa.idx]] <- dat
  b = round(seq(0,sd.plot.cap, length.out = 6))
  plot.tmp <- ggplot(datasets[[kappa.idx]], aes(x, y, fill=sqrt_sd)) +
    geom_raster() + theme_classic() +
    #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
    #scale_fill_viridis_c(option = "magma") +
    scale_fill_gradientn(limits = c(0,sd.plot.cap),
                         colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow","chartreuse1"),
                         breaks=b, labels=format(b)) +
    guides(fill=guide_legend(reverse=TRUE))  +
    geom_point(aes(x=0,y=0),colour="grey") +
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="grey") +
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="grey") +
    ggtitle(paste("SD Heatmap, Curvature: ", round(kappa.true,1)))
  #plot.tmp
  plots.list[[kappa.idx]] <- plot.tmp
  # ggplot(dat, aes(x, y, fill=manual_grad)) +
  #   geom_raster() + theme_classic() +
  #   scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  #   guides(fill=guide_legend(reverse=TRUE))
}


library(ggpubr)
dev.off()
ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
          plots.list[[4]], plots.list[[5]],  common.legend = TRUE)

save.plot = T
fig.height = 1000
fig.width = 2000
fig.res = 500
if(save.plot){
  ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]],
            plots.list[[4]], plots.list[[5]],  common.legend = TRUE) %>%
    ggexport(filename = "plots/theoretical_variance.eps", width = fig.width, height = fig.height, res = fig.res)
}




# max(datasets[[1]]$sqrt_sd)
# max(datasets[[2]]$sqrt_sd)
# max(datasets[[3]]$sqrt_sd)
# max(datasets[[4]]$sqrt_sd)
# max(datasets[[5]]$sqrt_sd)
#
#
# min(datasets[[1]]$sqrt_sd)
# min(datasets[[2]]$sqrt_sd)
# min(datasets[[3]]$sqrt_sd)
# min(datasets[[4]]$sqrt_sd)
# min(datasets[[5]]$sqrt_sd)
#
#
# hist(datasets[[1]]$sqrt_sd)
# hist(datasets[[2]]$sqrt_sd)
# hist(datasets[[3]]$sqrt_sd)
# hist(datasets[[4]]$sqrt_sd)
# hist(datasets[[5]]$sqrt_sd)
#
#
# plot(density(distances[[1]]))
# plot(density(distances[[2]]))
# plot(density(distances[[3]]))
# plot(density(distances[[4]]))
# plot(density(distances[[5]]))




