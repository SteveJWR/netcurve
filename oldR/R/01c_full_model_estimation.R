########################


# Script to ensure consistency of curvature estimations and 
# ensure that the constant curvature test is valid
#
# 
rm(list = ls())

library(lolaR)
source("R/00_functions.R")
source("R/clique_finder.R")
rm(filter_indices)

### If Running on a cluster
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

# Blocking the simulations for running on the cluster 
kappa.idx <- floor(id/20) + 1
block <-  id %% 20 
if(block == 0){
  block = 20 
}
n.sims = 10
sim.idx <- 1:n.sims

#Setting Seed for Replicability, but having a different 
#seed for parallel branches 

set.seed(id)

# Values of kappa to iterate 
kappa.set <- c(-2,-1,-0.5,0,0.5,1)
# scale parameter for the size of the network 
scale.set <- c(1/2,1/sqrt(2),1,2,4)

mu = -3
sd = 3


# Simulation Parameters 
sim.avg.variance <- 0.25**2 # average variance of GMM 
p = 3 # Latent Dimension of the data
num.midpoints = 3 #TODO: Change this to J 

res = 1 # Used as a tuning parameter for the approximate clique search

tri.const.seq <- (seq(0, 1, length.out = 21)) + 1 # Tuning paremeter set

# Simulated graph statistics 
graph.stat.names <- c("Graph size",
                      "Edge fraction", 
                      "Max Degree", 
                      "Mean Degree", 
                      "Distinct Cliques >= l", 
                      "Max Clique Size", 
                      "Mean Degree Centrality")

curve.scale <- 10 # TODO: remove, this will be a tuning parameter for the old version.

time.1 <- Sys.time()


# for(kappa.idx in seq(length(kappa.set))){
#   
# }

kappa = kappa.set[kappa.idx]
if(kappa < 0){
  centers.radius = 2.5 # 2
} else {
  centers.radius = 2.5
}

centers.variance = 0.5**2

kappa.ests.results <- matrix(NA,nrow = n.sims*length(tri.const.seq), 
                             ncol = length(scale.set) + 1)
sl.kappa.est.results <- matrix(NA,nrow = n.sims, 
                               ncol = length(scale.set) + 1) 
p.val.results <- matrix(NA, nrow = n.sims*length(tri.const.seq),
                        ncol = length(scale.set) + 1)
normalized.p.val.results <- p.val.results


p.val.boot.results <- matrix(NA, nrow = n.sims,
                        ncol = length(scale.set))

get.largest.clique = F 
for(scale.idx in seq(length(scale.set))){
  
  scale <- scale.set[scale.idx]
  n <- round(5000*scale)
  
  n.centers <- round(100*sqrt(scale))
  ell = round(8 + 4*log2(scale)) # min clique-size 
  approximate.variance <- sim.avg.variance        
  
  # TODO: Ensure that this tuning parameter is > 1 as a default. 

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
    # lpcm <- latent_position_cluster_model_2(n,n.centers, p, kappa, 
    #                                         centers.variance =centers.variance,
    #                                         cluster.variance = approximate.variance, 
    #                                         PI = PI)
    Z <- lpcm$Z
    
    nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
    nu.vec <- nu.vec*(nu.vec < 0)
    nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 
    
    
    A.sim <- sim_ls_network_fast_2(nu.vec, Z, kappa)
    d.sim <- colSums(A.sim)
    
    # if(sim == 1 & plot.graph = T){
    #   A.noiso <- A.sim[d.sim > 0, d.sim > 0]
    #   #A.noiso <- A
    #   A.noiso <- A.noiso[1:2000,1:2000]
    #   d.noiso <- colSums(A.noiso)
    #   A.noiso <- A.noiso[d.noiso > 0, d.noiso > 0 ]
    #   g <- igraph::graph_from_adjacency_matrix(A.noiso, mode = "undirected")
    #   V(g)$labels = NA
    #   plot(g, vertex.size= 2,vertex.label=NA)
    #   
    # }
    
    # don't compute max clique
    if(get.largest.clique) { 
      if(scale < 4 & (scale < 2 | kappa >= 0)){
        g.large <- igraph::graph_from_adjacency_matrix(A.sim, mode = "undirected")
        #full.cliques <- largest_cliques(g.large)
        max.clique.size <- clique_num(g.large)
        # print(max.clique.size)
      } else {
        max.clique.size <- NA
      }
    } else {
      max.clique.size <- NA
    }
    
    
  
    print(paste("Max cliques size:",max.clique.size))
    
    
    # making sure at least ~ 50 cliques are found. 
    # when scale > 8 we have to use an approximate clique search 
    if(scale > 1 | kappa == -2){
      clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels, 
                                      min_clique_size = ell)
    } else {
      #clique.set <- cluster_clique_search_4(A.sim,min_clique_size = ell,res = res)
      clique.set <- lolaR::ClusterCliqueSearch(A.sim,min_clique_size = ell,res = res, verbose = T)
    }
    
    clique.set <- clique_split(clique.set, min_clique_size = ell)
    
    print(paste("Number of Cliques of size,",ell,":", length(clique.set)))
    
    graph.stats[sim,1] <- n #Size 
    graph.stats[sim,2] <- sum(A.sim)/(length(A.sim) - n) # edge density
    graph.stats[sim,3] <- max(d.sim) # maximum degree
    graph.stats[sim,4] <- mean(d.sim) # mean degree
    graph.stats[sim,5] <- length(clique.set) #number of cliques found
    graph.stats[sim,6] <- max.clique.size # largest clique
    
    g <- igraph::graph_from_adjacency_matrix(A.sim, mode = "undirected")
    graph.stats[sim,7] <- mean(igraph::centr_eigen(g)$vector) # average eigenvector centrality
    if(length(clique.set) > 60 ){
      clique.set <- clique.set[1:60]
    }
    
    
    estimates = estimate_curvature(A.sim, clique.set, no.refit = F, d.yz.min = d.yz.min,d.yz.max = d.yz.max)
    
    D.hat <- estimates$D
    kappa.vec <- c()
    group.vec <- c()
    tri.const <- 1.6
    for(k in seq(3)){
      y.opt <- estimates$midpoints[k,1]
      z.opt <- estimates$midpoints[k,2]
      m.opt <- estimates$midpoints[k,3]
      x.set <- filter_indices_2(D.hat, y.opt,
                                z.opt,m.opt,
                                tri.const = tri.const)
      kap.block <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)
      kappa.vec <- c(kappa.vec, kap.block)
      group.vec <- c(group.vec, rep(k,length(kap.block)))

    }
    #plot(group.vec, kappa.vec)

  

    
    norm.curve.test <- normalized_constant_curvature_test_seq(estimates, num.midpoints = num.midpoints,
                                                              tri.const.seq = tri.const.seq,
                                                              curve.scale = curve.scale,
                                                              heavy.tail.scale = T)

    D.hat <- estimates$D
    if(kappa > 0){
      sl.kappa.est.results[sim,scale.idx] <- Lubold_Estimate_Curvature_Sphere(D.hat)
    } else if(kappa < 0){
      
      sl.kappa.est.results[sim,scale.idx] <- Lubold_Estimate_Curvature_Hyperbolic(D.hat)
      
    } 
    
      
    for(tri.const.idx in seq(length(norm.curve.test))){
      tri.const = tri.const.seq[tri.const.idx]
      res.idx <- n.sims*(tri.const.idx - 1) + sim



      # curve.test <- constant_curvature_test(estimates, num.midpoints = 3,
      #                                       tri.const = tri.const)
      #
      best.estimates <- norm.curve.test[[tri.const.idx]]$estimates[norm.curve.test[[tri.const.idx]]$estimates[,1] == 1, 2]
      p.val <- norm.curve.test[[tri.const.idx]]$p.value
      norm.p.val <- norm.curve.test[[tri.const.idx]]$norm.p.value
      
      if(!is.null(norm.curve.test[[tri.const.idx]]$estimates)){
        kappa.med <- median(best.estimates, na.rm = T)
        kappa.ests.results[res.idx,scale.idx] = kappa.med
        p.val.results[res.idx,scale.idx] = p.val
        normalized.p.val.results[res.idx,scale.idx] = norm.p.val
      }
      kappa.ests.results[res.idx,length(scale.set) + 1] = tri.const
      p.val.results[res.idx,length(scale.set) + 1] = tri.const
      normalized.p.val.results[res.idx,length(scale.set) + 1] = tri.const
    }
  }

  
  
  file.graph.stats <- paste0("results/graph_stats_kappa_",kappa,"_scale_",round(scale,1),"_block_",block,".csv")
  write.csv(graph.stats, file = file.graph.stats)
}

file.kappa.ests <- paste0("results/estimates_kappa_",kappa,"_block_",block,".csv")
file.sl.kappa.ests <- paste0("results/sl_estimates_kappa_",kappa,"_block_",block,".csv")

file.p.vals <- paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
file.norm.p.vals <- paste0("results/norm_p_vals_kappa_",kappa,"_block_",block,".csv")

write.csv(sl.kappa.est.results, file = file.sl.kappa.ests)
write.csv(kappa.ests.results, file = file.kappa.ests)
write.csv(p.val.results, file = file.p.vals)
write.csv(normalized.p.val.results, file = file.norm.p.vals)

time.2 <- Sys.time()

print(paste("Time Difference:", round(time.2 - time.1,3)))



tri.const <- 1.65
ps <- normalized.p.val.results[normalized.p.val.results[,5] == tri.const, 3]
round(ps,3)
plot(density(ps), xlim = c(0,1))

sl.kappa.est.results




# hist(p.val.results[,1])
# hist(normalized.p.val.results[,1])
# 
# plot(density(normalized.p.val.results[,1], na.rm = T,bw = 0.1), xlim = c(0,1))
# mean(normalized.p.val.results[,1] < 0.05, na.rm = T)
# 
# c0.thresh = 1.4
# idx <- which(normalized.p.val.results[,2] > c0.thresh)
# mean(normalized.p.val.results[idx,1] < 0.05, na.rm = T)
# hist(normalized.p.val.results[idx,1])
# 
# idx <- which(p.val.results[,2] > c0.thresh)
# mean(p.val.results[idx,1] < 0.05, na.rm = T)
# hist(p.val.results[idx,1])

