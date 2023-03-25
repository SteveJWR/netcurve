



# chose the true curvature. 
set.seed(6)
# kappa.true = 1
library(ggpubr)
source("00_functions.R")

max.dist.radius = 2.5
dyz = 2
location.scale = dyz/2

# TODO: Split up the simulation 
# over the different positions. 
# write a file which merges them all
# maybe have a more elegent way of doing this. 
n.sims = 40
n.grid = 70 # granularity of the grid spacing
sim.parallel = F
verbose = F 
soft.censoring.range = 10 # soft thresholds the result using Tanh
# this monotone transformation doesn't really affect the median. But allows for the limiting distribution to be 
# a bit more tractible. 
#B.boot <- 5

# some possibly reasonable values for the selection of the points going into the median estimator. 
c1 = .5
c2 = 2.5
c3 = .3


ell = 16

y = 1
z = 2
m = 3
x = 4


bias.scale = 2
sd.scale = 5
curve.scale = 10

breaks = sd.scale*c(0.0,0.25,0.75,1.5,2.5)
bias.breaks = bias.scale*c(-3,-1,-.5,0,.5,1,3)
bias.limits = bias.scale*c(-3,3)
limits = sd.scale*c(0,2.5)
censoring.range = Inf

# fixed parameters
p = 2
kappa.set <- c(-1,-0.5,0,0.5,1)
reference.widths <- rep(NA,length(kappa.set))

missing.js <- c()

dat.list <- list()

for(kappa.idx in seq(length(kappa.set))){
  data.block <- matrix(data = NA, nrow = 0, ncol = n.sims) 
  for(j in 1:3720){
    cat(paste("file:", j, "/", 3720), end = "\r")
    kappa.true = kappa.set[kappa.idx]
    file = paste0("results/sims/kappa_set_",j,"_kappa_",kappa.true,"_ell_",ell,".csv")
    data.tmp <- tryCatch({
      read.csv(file, header = T)
    }, error = function(e) {
      missing.js <<- c(missing.js,j)
      return(rep(NA,n.sims))})
    data.block <- rbind(data.block,as.vector(data.tmp))
  }
  dat.list[[kappa.idx]] = data.block 
}


bias.plots.list <- list()
sd.plots.list <- list()
rmse.plots.list <- list()
med.plots.list <- list()

#[30,31,32,33,60,61,70,71,85,86,91,92,119,120,144,145,150,151,163,164,172,173,179,180,205,206,226,227,259,260,275,276,296,297,304,305,321,322,335,336,344,345,357,358,370,371,398,399,408,409,418,419,424,425,432,433,437,438,440,441,443,444,445,454,455,460,461,462,474,475,476,491,492,496,497,501,502,507,508,520,521,527,528,536,537,538,539,544,545,546,547,555,556,569,570,576,577,578,585,586,602,603,604,605,612,613,622,623,624,625,641,642,656,657,661,662,666,667,674,675,683,684,695,696,698,699,700,721,722,723,724,736,737,742,743,751,752,755,756,762,763,774,775,782,783,800,801,812,813,814,831,832,834,835,844,845,846,848,849,857,858,859,860,868,869,875,876,877,885,886,895,896,908,909,915,916,921,922]

for(kappa.idx in seq(length(kappa.set))){
  kappa.true = kappa.set[kappa.idx]
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
  if(kappa.true == 0){
    ref.x = dyz/2
  } else if(kappa.true > 0){
    ref.x <- Z.init[1,2]
  } else if(kappa.true < 0){
    ref.x <- Z.init[2,2]
  }
  
  reference.widths[kappa.idx] <- ref.x
  
  
  
  if(kappa.true > 0 & pi/(sqrt(abs(kappa.true))) < max.dist.radius ){
    # turns into uniform over the whole sphere if it exceeds the range of the latent space. 
    latent.grid.radius = pi/(sqrt(abs(kappa.true))) - 10**(-6) # simulation radius of data 
  } else {
    latent.grid.radius = max.dist.radius
  }
  
  x.grid <- seq(-latent.grid.radius,latent.grid.radius, length.out = n.grid)
  
  # data points
  plot.dat = expand.grid(x.grid,x.grid)
  
  filtered.plot.idx = plot.dat[,1]**2 + plot.dat[,2]**2 <= latent.grid.radius**2
  
  plot.dat = plot.dat[filtered.plot.idx,]
  J <- nrow(plot.dat)
  data.block <- as.matrix(dat.list[[kappa.idx]])
  

  plot.dat[,3] <- rowMeans(scale_curvature(data.block,curve.scale) - scale_curvature(kappa.true,curve.scale), na.rm = T) # bias
  plot.dat[,4] <- rowSDs(scale_curvature(data.block,curve.scale)) # sd
  plot.dat[,5] <- sqrt(rowMeans((scale_curvature(data.block,curve.scale) - scale_curvature(kappa.true,curve.scale))^2,na.rm = T)) # rmse
  plot.dat[,6] <- rowMedians(scale_curvature(data.block,curve.scale) - scale_curvature(kappa.true,curve.scale)) # median-bias

  plot.dat <- data.frame(plot.dat)
  colnames(plot.dat) <- c("x","y","bias","sd","rmse","median_devn")
  
  #datasets[[kappa.idx]] <- dat
  
  bias.plot.tmp <- ggplot(data = plot.dat, aes(x = x,y = y, fill=bias)) + 
    geom_raster() + theme_classic() +
    scale_fill_gradientn(limits = bias.limits,
                         colours=c("navyblue","blue","#0000FFFF","#FFFFFFFF","#FF0000FF","red","darkred"),
                         breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    # scale_fill_gradientn(limits = c(0,0.8),
    #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
    #                      breaks=b, labels=format(b)) + 
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
    ggtitle(paste("Bias Heatmap, Curvature: ", kappa.true))
  
  sd.plot.tmp <- ggplot(data = plot.dat, aes(x = x,y = y, fill=sd)) + 
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
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
    ggtitle(paste("SD Heatmap, Curvature: ", kappa.true))
  
  rmse.plot.tmp <- ggplot(data = plot.dat, aes(x = x,y = y, fill=rmse)) + 
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
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
    ggtitle(paste("RMSE Heatmap, Curvature: ", kappa.true))
  
  med.plot.tmp <- ggplot(data = plot.dat, aes(x = x,y = y, fill=median_devn)) + 
    geom_raster() + theme_classic() +
    scale_fill_gradientn(limits = bias.limits,
                         colours=c("navyblue","blue","#0000FFFF","#FFFFFFFF","#FF0000FF","red","darkred"),
                         breaks=bias.breaks, labels=format(bias.breaks)) +
    #scale_fill_viridis_c(option = "magma") + 
    # scale_fill_gradientn(limits = c(0,0.8),
    #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
    #                      breaks=b, labels=format(b)) + 
    guides(fill=guide_legend(reverse=TRUE))  + 
    geom_point(aes(x=0,y=0),colour="black") + 
    geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
    geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
    ggtitle(paste("Median Heatmap, Curvature: ", kappa.true))
  
# 
#   selected.plot.tmp <- ggplot(data = plot.dat, aes(x = x,y = y, fill=selection)) + 
#     geom_raster() + theme_classic() +
#     scale_fill_gradientn(limits = c(0,1),
#                          colours=c("darkred","navyblue"),
#                          breaks=c(0,1), labels=c(0,1)) +
#     #scale_fill_viridis_c(option = "magma") + 
#     # scale_fill_gradientn(limits = c(0,0.8),
#     #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
#     #                      breaks=b, labels=format(b)) + 
#     guides(fill=guide_legend(reverse=TRUE))  + 
#     geom_point(aes(x=0,y=0),colour="black") + 
#     geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
#     geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
#     ggtitle(paste("Selection Plot, Curvature: ", kappa.true))
  
  
  
  bias.plots.list[[kappa.idx]] <- bias.plot.tmp
  sd.plots.list[[kappa.idx]] <- sd.plot.tmp
  rmse.plots.list[[kappa.idx]] <- rmse.plot.tmp
  med.plots.list[[kappa.idx]] <- med.plot.tmp
  
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

text <- paste("Median-Deviation map", "\n", "Clique Size:",ell)

# Create a text grob
med.tgrob <- text_grob(text,size = 20)
med.plot_text <- as_ggplot(med.tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))

med.plot <- ggarrange(med.plot_text,med.plots.list[[1]], med.plots.list[[2]], med.plots.list[[3]], 
                      med.plots.list[[4]], med.plots.list[[5]],  common.legend = TRUE)




# ggsave(paste0("results/plots/bias_kappa_",kappa.true,"_ell_",ell,".png"),bias.plot.tmp)
# ggsave(paste0("results/plots/sd_kappa_",kappa.true,"_ell_",ell,".png"),sd.plot.tmp)
# ggsave(paste0("results/plots/rmse_kappa_",kappa.true,"_ell_",ell,".png"),rmse.plot.tmp)
# ggsave(paste0("results/plots/med_kappa_",kappa.true,"_ell_",ell,".png"),med.plot.tmp)
# ggsave(paste0("results/plots/selection_kappa_",kappa.true,"_ell_",ell,".png"),selected.plot.tmp)
# 



