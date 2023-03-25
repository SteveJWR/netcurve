

rm(list = ls())

source("00_functions.R")
source("clique_finder.R")
set.seed(1)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}
RNGkind("L'Ecuyer-CMRG")


n.sims = 10 

plot.width = 30
plot.height = 14
plot.dpi = 600
#units = "cm"




# probability of ending up on which side of the sphere 

mu = -3
sd = 3

kappa1 = 0.8
kappa2 = 0.5

res = 1.1

curve.scale <- 10 

Time.steps = 50
t.change1 = 15
t.change2 = 38

approximate.variance <- 0.25**2
d.yz.min <- 1.5
centers.radius = 2.5
rho = 0.5


p1 = 2
p2 = 2
q = .5

scale <- 1
scale.set <- c(1/sqrt(2),1,2,4,8)

ell.set <- round(8 + 4*log2(scale.set))
c1 = .5
c2 = 2
c3 = .25

# rule of thumb for number of connections




# the problem is that this seems to be a hard geometry to search for cliques
#mu <- -4  # mu for spherical ish geometry 
p.val.results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
colnames(p.val.results) = paste0("CliqueSize_", ell.set)

time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){
  
  scale <- scale.set[scale.idx]
  
  n <- round(5000*scale)
  # n <- round(1000*sqrt(scale))
  
  n.centers1 <- round(50*sqrt(scale))
  n.centers2 <- round(50*sqrt(scale))
  
  ell = round(8 + 4*log2(scale)) # min clique-size 
  
  d.yz.min = 1.5
  if(ell < 8){
    d.yz.min = 1
  }
  #d.yz.max = -log(10/ell^2)
  d.yz.max = max(log(ell),log(ell^2/10)) 
  for(sim in seq(n.sims)){
    
    PI1 <- as.numeric(rdirichlet(1, rep(2,n.centers1)))
    PI2 <- as.numeric(rdirichlet(1, rep(2,n.centers2)))
    
    D <- connected_spheres_lpcm(n, n.centers1, n.centers2, 
                                p1,p2,PI1, PI2, 
                                q, kappa1,kappa2, 
                                approximate.variance, max.rad = centers.radius)
    
    
    nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
    nu.vec <- nu.vec*(nu.vec < 0 )
    nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 

    A.sim <- sim_ls_network(nu.vec, D)
    
    if(scale > 1){
      clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels, 
                                       min_clique_size = ell)
    } else {
      clique.set <- cluster_clique_search_4(A.sim,res = 1.2, min_clique_size = ell)
    }
    
    print(paste("Number of Cliques of size,",ell,":", length(clique.set)))
    if(length(clique.set) > 60 ){
      clique.set <- clique.set[1:60]
    }
    estimates = estimate_curvature(A.sim, clique.set)

    
    D.hat <- estimates$D
    mid.search <- estimates$midpoints
    
    
    y.opt = mid.search[1,1]
    z.opt = mid.search[1,2]
    m.opt = mid.search[1,3]
    
    x.set <- filter_indices(D.hat, y.opt,
                            z.opt,m.opt, 
                            c1 = c1,c2 = c2,c3 = c3)

    kappa.set.1 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)

    y.opt = mid.search[2,1]
    z.opt = mid.search[2,2]
    m.opt = mid.search[2,3]

    x.set <- filter_indices(D.hat, y.opt,
                            z.opt,m.opt, 
                            c1 = c1,c2 = c2,c3 = c3)

    kappa.set.2 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)

    y.opt = mid.search[3,1]
    z.opt = mid.search[3,2]
    m.opt = mid.search[3,3]
    
    x.set <- filter_indices(D.hat, y.opt,
                            z.opt,m.opt, 
                            c1 = c1,c2 = c2,c3 = c3)

    kappa.set.3 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)
    
    
    kappa.set.1 <- kappa.set.1[!is.na(kappa.set.1)]
    kappa.set.2 <- kappa.set.2[!is.na(kappa.set.2)]
    kappa.set.3 <- kappa.set.3[!is.na(kappa.set.3)]
    
    index <- c(rep(1,length(kappa.set.1)),
               rep(2,length(kappa.set.2)),
               rep(3,length(kappa.set.3)))
    kappa.vec <- c(kappa.set.1,
                   kappa.set.2,
                   kappa.set.3)
    
    y = scale_curvature(kappa.vec,curve.scale)
    x = index
    
    # good, this test basically rejects every time
    
    test.dat <- data.frame("loc" = x, "est" = kappa.vec)
    test <- kruskal.test(est ~ loc, data = test.dat) 
    p.val.results[sim,scale.idx] = test$p.value
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
    
  }
  
  
  
  # Example plot of a changepoint problem: 
  # save one for each sample size 
  
  ### 
 
  ## simple ploting of changes and smoothed profile
  # plot(y.clean, pch=20, col="black")
  # lines(res.l1$smt, col="red", lwd=2)
  # lines(y.true.clean, lty=4, col="green", lwd=2)
  # abline(v=cpt, lty=2, col="red")
  # abline(v=cpt.true, lty=3, col="green")
  
  plot.file <- paste0("plots/AdjSpheres_CliqueSize_",ell,"_block_",id,".png")
  plot.title <- paste0("Adjacent Spheres Curvature Clique Size: ",ell)
  
  
  labels <- c(rep(0,length(kappa.set.1)),rep(1,length(kappa.set.2)),rep(2,length(kappa.set.3)))
  bp.dat <- data.frame(matrix(c(labels, 
                                scale_curvature(kappa.set.1,curve.scale),
                                scale_curvature(kappa.set.2,curve.scale),
                                scale_curvature(kappa.set.3,curve.scale)), ncol = 2))
  colnames(bp.dat) <- c("group", "curvature")
  

  
  med.1 <- median(scale_curvature(kappa.set.1,curve.scale),na.rm = T)
  med.2 <- median(scale_curvature(kappa.set.2,curve.scale),na.rm = T)
  med.3 <- median(scale_curvature(kappa.set.3,curve.scale),na.rm = T)
  
  kappa.vec <- c(kappa.set.1,
                 kappa.set.2,
                 kappa.set.3)
  med.group <- median(scale_curvature(kappa.vec,curve.scale),na.rm = T)
  
  plt <- ggplot(bp.dat, aes(y = curvature, x = group)) + 
    geom_jitter() + 
    geom_vline(xintercept = 0.5, col = "red") + 
    geom_vline(xintercept = 1.5, col = "red") + 
    geom_segment(aes(x=-0.5,xend=0.5,y=med.1,yend=med.1), col = "blue", linetype = "dashed") + 
    geom_segment(aes(x=0.5,xend=1.5,y=med.2,yend=med.2), col = "blue", linetype = "dashed") + 
    geom_segment(aes(x=1.5,xend=2.5,y=med.3,yend=med.3), col = "blue", linetype = "dashed") + 
    geom_segment(aes(x=-0.5,xend=2.5,y=med.group,yend=med.group), col = "green", linetype = "dashed") + 
    xlab("Location") + 
    ylab("Soft Threshold Curvature") + 
    ggtitle(plot.title)
  
 
  
  ggsave(
    plot.file,
    plot = plt,
    device = NULL,
    path = NULL,
    scale = 1,
    width = plot.width,
    height = plot.height,
    units = "cm",
    dpi = plot.dpi
  )
}

csv.file <- paste0("results/adjacent_spheres_results","_block_",id,".csv")
write.csv(p.val.results, file = csv.file)

time.2 <- Sys.time()

print(paste("Time Difference:", time.2 - time.1))









