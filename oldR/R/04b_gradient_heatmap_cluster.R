

# chose the true curvature. 
set.seed(6)
# kappa.true = 1
library(ggpubr)
source("00_functions.R")

max.dist.radius = 6
dyz = 2
location.scale = dyz/2

# TODO: Split up the simulation 
# over the different positions. 
# write a file which merges them all
# maybe have a more elegent way of doing this. 
n.sims = 25
n.grid = 50 # granularity of the grid spacing
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


ell = 8

y = 1
z = 2
m = 3
x = 4


bias.scale = 2
sd.scale = 5

breaks = sd.scale*c(0.0,0.25,0.75,1.5,2.5)
bias.breaks = bias.scale*c(-3,-1,-.5,0,.5,1,3)
bias.limits = bias.scale*c(-3,3)
limits = sd.scale*c(0,2.5)
censoring.range = Inf

# fixed parameters
p = 2
kappa.set <- c(-1,-0.5,0,0.5,1)
reference.widths <- rep(NA,length(kappa.set))

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  kappa.idx = 1
  
} else {
  # coerce the value to an integer
  kappa.idx <- as.numeric(slurm_arrayid)
}

if(kappa.idx > 5){
  ell = 16
  kappa.idx = kappa.idx - 5
}

kappa.true <- kappa.set[kappa.idx]

epsilon <- 10**(-6) # rounding piece for projection of the sphere 

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
dat[,7] <- rep(NA,J) # Accept region 



clique.set = matrix(1:(4*ell), ncol = ell, byrow = T)
clique.vec <- rep(1:4, each = ell)
K <- nrow(clique.set)
nus.hat.mat <- matrix(nrow = K, ncol = ell)

for(cl in 1:K){
  clique.idx <- clique.set[cl,]
  nus.hat <- rep(0,ell)
  nus.hat.mat[cl,] <- nus.hat
}
nus.hat.vec <- as.vector(t(nus.hat.mat))
for(j in 1:J){
  #cat(end = "\n")
  cat(paste("Locations:", j,"/",J), end = '\r')
  #cat(end = "\n")
  if(kappa.true > 0 ){
    d0 = sqrt(dat[j,1]**2 + dat[j,2]**2)
    z0 = cos(sqrt(kappa.true)*d0)
    if(d0 != 0 ){
      scl <- sqrt((1 - z0^2)/(dat[j,1]^2 + dat[j,2]^2))
    } else {
      scl <- 0
    }
    
    Z.new <- c(z0, scl*dat[j,1],scl*dat[j,2])
    
  } else if(kappa.true < 0) { 
    d0 = sqrt(dat[j,1]**2 + dat[j,2]**2)
    z0 = cosh(sqrt(-kappa.true)*d0)
    if(d0 != 0 ){
      scl <- sqrt((z0^2 - 1)/(dat[j,1]^2 + dat[j,2]^2))
    } else {
      scl <- 0
    }
    
    Z.new <- c(z0, scl*dat[j,1],scl*dat[j,2])
    
  } else if(kappa.true == 0){
    Z.new <- c(dat[j,1],dat[j,2])
    
  }
  Z.block <- rbind(Z.init, Z.new)
  Z <- Z.block[rep(1:nrow(Z.block), each = ell), ]
  D.block <- pos_to_dist(Z.block,kappa.true)
  # min.off.diag = min(D.block + diag(100,4) )
  # if(min.off.diag < 0.1){
  #   next 
  # }
  D.true <- pos_to_dist(Z,kappa.true)
  
  rand.eff <- rep(0,4*ell)
  kappa.hat.set <- rep(0,n.sims)
  #bias.corrected.kappa.hat.set <- kappa.hat.set
  #reverse.bias.corrected.kappa.hat.set <- kappa.hat.set
  
  i1 = (c1)*D.block[y,z] <= D.block[x,y]
  i2 = (c1)*D.block[y,z] <= D.block[x,z]
  i3 = D.block[x,z] <= (c2)*D.block[y,z]
  i4 = D.block[x,y] <= (c2)*D.block[y,z]
  i5 = abs(D.block[x,y] - D.block[x,z]) <= 2*c3*(D.block[y,z])
  
  # selection region
  is.selected <- ifelse(i1*i2*i3*i4*i5 == 1, T, F)

  for(sim in 1:n.sims){
    A <- sim_ls_network(rand.eff,D.true)

    
    # D.est = D_estimate(A,clique.set,fixed.effects = nus.hat.mat)
    # D.hat = D.est$estimates
    
    # sometimes the initial point is good, sometimes it is bad 
    # is there a good way to do this 
    #TODO: discuss with Tyler
    # how to initialize this. 
    # A heuristic based on average number of connections 
    # across cliques? 
    # this is a heuristic approach based on connections.
    # a better way is to trim all large distances 
    
    #D0 = matrix(2,4,4)
    
    D0 = init_D0(A,clique.vec,nus.hat.vec)
    
    D.hat = estimate_D_restricted_fast(A,clique.vec,nus.hat.vec,thresh = 10**(-2),D0 = D0,verbose = verbose, solver = "OSQP")
    
    dxy.hat = D.hat[4,1]
    dxz.hat = D.hat[4,2]
    dyz.hat = D.hat[2,1]
    dxm.hat = D.hat[4,3]
    cat(paste("Locations:", j,"/",J), end = '\r')
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
    } else {
      kappa.hat <- NA
    }
    kappa.hat.set[sim] <- kappa.hat
    
  }
  
  
  # dat[j,3] <- mean(scale_curvature(kappa.hat.set,soft.censoring.range) ,na.rm = T) - kappa.true # bias
  # dat[j,4] <- sd(scale_curvature(kappa.hat.set,soft.censoring.range) ,na.rm = T)# sd 
  # dat[j,5] <- sqrt(mean((scale_curvature(kappa.hat.set,soft.censoring.range)  - scale_curvature(kappa.true,soft.censoring.range))**2,na.rm = T))# rmse 
  # dat[j,6] <- median(scale_curvature(kappa.hat.set,soft.censoring.range) ,na.rm = T) - kappa.true # bias
  # dat[j,7] <- is.selected
  dat = kappa.hat.set
}

if(kappa.true == 0){
  ref.x = dyz/2
} else if(kappa.true > 0){
  ref.x <- Z.init[1,2]
} else if(kappa.true < 0){
  ref.x <- Z.init[2,2]
}
reference.widths[kappa.idx] <- ref.x
dat <- data.frame(dat)
colnames(dat) <- c("x","y","bias","sd","rmse","median","selection")

#datasets[[kappa.idx]] <- dat

bias.plot.tmp <- ggplot(data = dat, aes(x = x,y = y, fill=bias)) + 
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

sd.plot.tmp <- ggplot(data = dat, aes(x = x,y = y, fill=sd)) + 
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

rmse.plot.tmp <- ggplot(data = dat, aes(x = x,y = y, fill=rmse)) + 
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

med.plot.tmp <- ggplot(data = dat, aes(x = x,y = y, fill=median)) + 
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

dat$selection <- 1*dat$selection
selected.plot.tmp <- ggplot(data = dat, aes(x = x,y = y, fill=selection)) + 
  geom_raster() + theme_classic() +
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("darkred","navyblue"),
                       breaks=c(0,1), labels=c(0,1)) +
  #scale_fill_viridis_c(option = "magma") + 
  # scale_fill_gradientn(limits = c(0,0.8),
  #                      colours=c("navyblue", "darkmagenta", "red", "darkorange1", "yellow"),
  #                      breaks=b, labels=format(b)) + 
  guides(fill=guide_legend(reverse=TRUE))  + 
  geom_point(aes(x=0,y=0),colour="black") + 
  geom_point(aes(x=reference.widths[!!kappa.idx],y=0),colour="black") + 
  geom_point(aes(x=-reference.widths[!!kappa.idx],y=0),colour="black") + 
  ggtitle(paste("Selection Plot, Curvature: ", kappa.true))

plot.list <- list(bias.plot.tmp,sd.plot.tmp,
                  rmse.plot.tmp,med.plot.tmp,
                  selected.plot.tmp)

ggsave(paste0("results/plots/bias_kappa_",kappa.true,"_ell_",ell,".png"),bias.plot.tmp)
ggsave(paste0("results/plots/sd_kappa_",kappa.true,"_ell_",ell,".png"),sd.plot.tmp)
ggsave(paste0("results/plots/rmse_kappa_",kappa.true,"_ell_",ell,".png"),rmse.plot.tmp)
ggsave(paste0("results/plots/med_kappa_",kappa.true,"_ell_",ell,".png"),med.plot.tmp)
ggsave(paste0("results/plots/selection_kappa_",kappa.true,"_ell_",ell,".png"),selected.plot.tmp)



