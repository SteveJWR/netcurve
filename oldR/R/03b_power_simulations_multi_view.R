
#load("3b_coverage.RData")
#save.image(file = "3b_coverage.RData")
rm(list = ls())

source("00_functions.R")
source("clique_finder.R")



slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(id)

n.sims = 10 
n.sims <- 2
sim.idx <- 1:n.sims



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

sim.avg.variance <- 0.25**2
d.yz.min <- 1.5
centers.radius = 2.5
rho = 0.5


p = 3
q = .5

scale <- 1
scale.set <- c(1/sqrt(2),1,2,4)



ell.set <- round(8 + 4*log2(scale.set))

num.midpoints = 3
tri.const.seq = (seq(0, 1, length.out = 21)) + 1
curve.scale = 10

tri.const = 1.5



# rule of thumb for number of connections




# the problem is that this seems to be a hard geometry to search for cliques
#mu <- -4  # mu for spherical ish geometry 
p.val.results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
colnames(p.val.results) = paste0("CliqueSize_", ell.set)

time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){
  
  scale <- scale.set[scale.idx]
  
  n <- round(5000*scale)
  # n <- round(4000*sqrt(scale))
  
  n.centers <- round(100*sqrt(scale))
  #n.centers <- round(100*scale)
  
  ell = round(8 + 4*log2(scale)) # min clique-size 
  
  #approximate.variance <- sim.avg.variance/max(scale^2,1)
  approximate.variance <- sim.avg.variance        
  
  d.yz.min = 1.5
  if(ell < 8){
    d.yz.min = 1
  }
  #d.yz.max = -log(10/ell^2)
  d.yz.max = max(log(ell),log(ell^2/10)) 
  for(sim in seq(n.sims)){
    
    PI <- as.numeric(rdirichlet(1, rep(2,n.centers)))
    
    
    
    cluster.model.variance = rgamma(n.centers, shape = approximate.variance)
    
    lpcm <- latent_position_cluster_model(n,n.centers, p, 
                                          centers.radius, 
                                          kappa1, 
                                          cluster.model.variance, 
                                          PI = PI)
    
    Z <- lpcm$Z
    
    nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
    nu.vec <- nu.vec*(nu.vec < 0)
    nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 
    
    # D1 <- pos_to_dist(Z, kappa1)
    # D2 <- pos_to_dist(Z, kappa2)
    # 
    # A.sim1 <- sim_ls_network(nu.vec, D1)
    # A.sim2 <- sim_ls_network(nu.vec, D2)
    # 
    # A.sim1 <- as(A.sim1, "sparseMatrix")
    # A.sim2 <- as(A.sim2, "sparseMatrix")
    
    
    A.sim1 <- sim_ls_network_fast_2(nu.vec, Z, kappa1)
    A.sim2 <- sim_ls_network_fast_2(nu.vec, Z, kappa2)
    
    # P.sim <- exp(outer(nu.vec, nu.vec, "+") -D)
    # A.sim <- 1*(P.sim >= 0.5)
 
    
    # making sure at least ~ 50 cliques are found. 
    # when scale > 8 we have to use an approximate clique search 
    # if(1 <= scale & scale < 4 & kappa1 <= 0 & kappa2 <= 0){
    #   
    #   clique.set <- clique_finder_2(A.sim, ell, 2500)
    # 
    # } else {
    # trim the network 
    
    #A.trim[d.sim <= quantile(d.sim, 0.2), ] = 0
    #A.trim[,d.sim <= quantile(d.sim, 0.2)] = 0
    # the trim is meant to speed up the clique search 
    
    #ell = 16
    #clique.set1 <- clique_finder(A.sim1, min_clique_size = ell,num_cliques_stop = 70, num_iters = 2000)
    # TODO: mention this variant on the clique_finding algorithm 
    
    # fast-ish methods of finding cliques 
    
    if(scale > 1){
      clique.set1 <- guided_clique_set(A.sim1,lpcm$cluster_labels, 
                                        min_clique_size = ell)
      clique.set2 <- guided_clique_set(A.sim2,lpcm$cluster_labels, 
                                        min_clique_size = ell)
    } else {
      clique.set1 <- cluster_clique_search_4(A.sim1,res = 1.2, min_clique_size = ell)
      clique.set2 <- cluster_clique_search_4(A.sim2,res = 1.2, min_clique_size = ell)
    }
    
    clique.set1 <- clique_split(clique.set1, min_clique_size = ell)
    clique.set2 <- clique_split(clique.set2, min_clique_size = ell)
    
    print(paste("Number of Cliques of size,",ell,":", length(clique.set1)))
    print(paste("Number of Cliques of size,",ell,":", length(clique.set2)))
    
    if(length(clique.set1) > 6 & length(clique.set2) > 6){
      if(length(clique.set1) > 60 ){
        clique.set1 <- clique.set1[1:60]
      }
      #tri.const = 1.5
      tri.const = 1.45
      estimates1 = estimate_curvature(A.sim1, clique.set1, tri.const = tri.const)
      
      
      #clique.set2 <- cluster_clique_search(A.sim2,res1 = 1,res2 = 1, min_clique_size = ell)
      #clique.set <- clique_split(clique.set, min_clique_size = ell)
      # }
      
      
      if(length(clique.set2) > 60 ){
        clique.set2 <- clique.set2[1:60]
      }
      estimates2 = estimate_curvature(A.sim2, clique.set2, tri.const = tri.const)
      
      
      
      
      kappa.vec <- c(estimates1$kappas, estimates2$kappas)
      kappa.vec <- as.numeric(kappa.vec)
      trim.kappa.vec <- scale_curvature(kappa.vec, c = curve.scale)
      index.vec <- c(rep(1,length(estimates1$kappas)),rep(2,length(estimates2$kappas)))
      if(length(unique(index.vec)) > 1){
        at.least.two.groups <- length(as.numeric(table(index.vec))) > 1 
        if(at.least.two.groups){
          trim.test.dat <- data.frame("loc" = index.vec, "trim.est" = trim.kappa.vec,"est" = kappa.vec)
          # normalized test
          norm.test <- tryCatch({
            kruskal.test(trim.est ~ loc, data = trim.test.dat) 
          }, error = function(error_condition) {
            return(NULL)
          }) 
          
          
        }  
        if(!is.null(norm.test$p.value)){
          print( norm.test$p.value)
          p.val.results[sim,scale.idx] = norm.test$p.value
        }
      }
 
    }
    
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
  
  
  
  kappa.est.1 <- estimates1$kappas
  kappa.est.2 <- estimates2$kappas
  
  if(length(kappa.est.1) == 0){
    kappa.est.1 == c()
  }
  if(length(kappa.est.2) == 0){
    kappa.est.2 == c()
  }
  plot.file <- paste0("plots/Multiview_CliqueSize_",ell,"_block_",id,".png")
  plot.title <- paste0("Multiview Curvature Clique Size: ",ell)
  
  
  labels <- c(rep(0,length(kappa.est.1)),rep(1,length(kappa.est.2)))
  bp.dat <- data.frame(matrix(c(labels, 
                                scale_curvature(kappa.est.1,curve.scale),
                                scale_curvature(kappa.est.2,curve.scale)), ncol = 2))
  colnames(bp.dat) <- c("group", "curvature")
  
  
  
  med.1 <- median(scale_curvature(kappa.est.1,curve.scale),na.rm = T)
  med.2 <- median(scale_curvature(kappa.est.2,curve.scale),na.rm = T)
  
  kappa.vec <- c(kappa.est.1,
                 kappa.est.2)
  med.group <- median(scale_curvature(kappa.vec,curve.scale),na.rm = T)
  
  plt <- ggplot(bp.dat, aes(y = curvature, x = group)) + 
    geom_jitter() + 
    geom_vline(xintercept = 0.5, col = "red") + 
    geom_segment(aes(x=-0.5,xend=0.5,y=med.1,yend=med.1), col = "blue", linetype = "dashed") + 
    geom_segment(aes(x=0.5,xend=1.5,y=med.2,yend=med.2), col = "blue", linetype = "dashed") + 
    geom_segment(aes(x=-0.5,xend=1.5,y=med.group,yend=med.group), col = "green", linetype = "dashed") + 
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

csv.file <- paste0("results/multiview_results","_block_",id,".csv")
write.csv(p.val.results, file = csv.file)

time.2 <- Sys.time()

print(paste("Time Difference:", time.2 - time.1))


