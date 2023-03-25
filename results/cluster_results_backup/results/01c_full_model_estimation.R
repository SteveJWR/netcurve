########################



# latent distribution 
#
# 
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


n.sims = 10
sim.idx <- 1:n.sims


set.seed(id)


kappa.set <- c(-2,-1,-0.5,0,0.5,1)
scale.set <- c(1/sqrt(2),1,2,4)


mu = -3
sd = 3

kappa.true = .5

sim.avg.variance <- 0.25**2
p = 3


res = 1

c1 = .5
c2 = 2
c3 = .25




graph.stat.names <- c("Graph size",
                      "Edge fraction", 
                      "Max Degree", 
                      "Mean Degree", 
                      "Distinct Cliques >= l", 
                      "Max Clique Size")

curve.scale <- 10 
time.1 <- Sys.time()
# todo: 
# find a corresponding 

for(kappa.idx in seq(length(kappa.set))){
  kappa = kappa.set[kappa.idx]
  if(kappa < 0){
    centers.radius = 2
  } else {
    centers.radius = 2.5
  }
  kappa.ests.results <- matrix(NA,nrow = n.sims, ncol = length(scale.set))
  p.val.results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
  for(scale.idx in seq(length(scale.set))){
    
    scale <- scale.set[scale.idx]
    n <- round(5000*scale)
    n.centers <- round(100*sqrt(scale))
    ell = round(8 + 4*log2(scale)) # min clique-size 
    approximate.variance <- sim.avg.variance        
    
    d.yz.min = 1.5
    if(ell < 8){
      d.yz.min = 1
    }
    #d.yz.max = -log(10/ell^2)
    d.yz.max = max(log(ell),log(ell^2/10)) 
    
    graph.stats <- matrix(NA, nrow = n.sims, ncol = length(graph.stat.names))
    colnames(graph.stats) <- graph.stat.names
    
    for(sim in sim.idx){
      print(paste("Simulation:", sim, "/",n.sims))
      
      PI <- as.numeric(rdirichlet(1, rep(2,n.centers)))
      
      
      
      cluster.model.variance = rgamma(n.centers, shape = approximate.variance)
      
      lpcm <- latent_position_cluster_model(n,n.centers, p, 
                                            centers.radius, 
                                            kappa, 
                                            cluster.model.variance, 
                                            PI = PI)
      
      Z <- lpcm$Z
      
      nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
      nu.vec <- nu.vec*(nu.vec < 0)
      nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 
      
      
      A.sim <- sim_ls_network_fast_2(nu.vec, Z, kappa)
      d.sim <- colSums(A.sim)
      
      if(scale < 4){
        
        g.large <- igraph::graph_from_adjacency_matrix(A.sim, mode = "undirected")
        
        #full.cliques <- largest_cliques(g.large)
        max.clique.size <- clique_num(g.large)
        # print(max.clique.size)
      } else {
        max.clique.size <- NA
      }
      
      #ell <- max.clique.size - 4
      
      if((max.clique.size < ell + 1) & scale < 4){
        next
      }
      print(paste("Max cliques size:",max.clique.size))
      
      
      # making sure at least ~ 50 cliques are found. 
      # when scale > 8 we have to use an approximate clique search 
      if(scale > 1){
        clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels, 
                                         min_clique_size = ell)
      } else {
        clique.set <- cluster_clique_search_4(A.sim,min_clique_size = ell,res = res)
      }
      
      clique.set <- clique_split(clique.set, min_clique_size = ell)
      
      print(paste("Number of Cliques of size,",ell,":", length(clique.set)))
      
      graph.stats[sim,1] <- n
      graph.stats[sim,2] <- sum(A.sim)/(length(A.sim) - n)
      graph.stats[sim,3] <- max(d.sim)
      graph.stats[sim,4] <- mean(d.sim)
      graph.stats[sim,5] <- length(clique.set)
      graph.stats[sim,6] <- max.clique.size
      
      #
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
      test.dat <- data.frame("loc" = x, "est" = kappa.vec)
      test <- kruskal.test(est ~ loc, data = test.dat) 
      kappa.ests.results[sim,scale.idx] = estimates$kappa.med
      p.val.results[sim,scale.idx] = test$p.value
    }
    
    file.graph.stats <- paste0("results/graph_stats_kappa_",kappa,"_scale_",round(scale,1),"_block_",id,".csv")
    write.csv(graph.stats, file = file.graph.stats)
  }
  file.kappa.ests <- paste0("results/estimates_kappa_",kappa,"_block_",id,".csv")
  file.p.vals <- paste0("results/p_vals_kappa_",kappa,"_block_",id,".csv")
  write.csv(kappa.ests.results, file = file.kappa.ests)
  write.csv(p.val.results, file = file.p.vals)
}

time.2 <- Sys.time()

print(paste("Time Difference:", round(time.2 - time.1,3)))



wilcox.test(kappa.set.1,kappa.set.3, paired=FALSE)




