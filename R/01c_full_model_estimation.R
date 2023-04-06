########################
#devtools::install_github("SteveJWR/lolaR", ref = "main")

# Script to ensure consistency of curvature estimations and
# ensure that the constant curvature test is valid
#
#
rm(list = ls())

library(lolaR)
source("R/00_functions.R")

#whether to write the files
write.files = T
write.graph.stats.files = T

### If Running on a cluster
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}
#Setting seed for reproducibility, but having a different
#seed for parallel branches
if(id == 15){
  set.seed(-id)
} else{
  set.seed(id)
}

# Simulation Parameters
# Values of kappa to iterate
kappa.set <- c(-2,-1,-0.5,0,0.5,1)

# Partitioning the simulations for running on the cluster

block <-  floor((id - 1)/length(kappa.set)) +1
kappa.idx <- (id - 1) %% length(kappa.set) + 1


n.sims = 10
sim.idx <- 1:n.sims




kappa = kappa.set[kappa.idx]
# Radius of the latent GMM
if(kappa < 0){
  centers.radius = 2.5 # 2
} else {
  centers.radius = 2.5
}
centers.variance = 0.5**2

# scale parameter for the size of the network
scale.set <- c(1/sqrt(2),1,2) #c(1/sqrt(2),1,2,4) #corresponds to ell = 6,8,12

# parameters for the distribution of the latent traits
mu = -3
sd = 3

sim.avg.variance <- 0.25**2 # average variance the variance of the clusters of the GMM
p = 3 # Latent Dimension of the data


# Method Tuning Parameters
num.midpoints = 3
tri.const = 1.4
tri.const.seq <- (seq(0, 1, length.out = 21)) + 1 # Tuning parameter set
max.num.cliques = 35
num.subsamples = 250
max.iter.estimate = 6
d.yz.min = 1

# Recorded simulated graph statistics
graph.stat.names <- c("Graph size",
                      "Edge fraction",
                      "Max Degree",
                      "Mean Degree",
                      "Distinct Cliques >= l",
                      "Max Clique Size",
                      "Mean Degree Centrality")



# matrix arrays for storing results
kappa.ests.results <- matrix(NA,nrow = n.sims*length(tri.const.seq),
                             ncol = length(scale.set) + 1)
# sl.kappa.est.results <- matrix(NA,nrow = n.sims,
#                                ncol = length(scale.set) + 1)
p.val.results <- matrix(NA, nrow = n.sims*length(tri.const.seq),
                        ncol = length(scale.set) + 1)
# normalized.p.val.results <- p.val.results

p.val.sub.results <- matrix(NA, nrow = n.sims,
                            ncol = length(scale.set))


# Additional Simulation Parameters
get.largest.clique = F

time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){

  scale <- scale.set[scale.idx]


  n <- round((3/4)*5000*scale)
  n.centers <- round((3/4)*100*sqrt(scale))

  ell = round(8 + 4*log2(scale)) # min clique-size
  approximate.variance <- sim.avg.variance




  #Rule of thumb for d.yz.max
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

    # Just search for cliques
    clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels,
                                    min_clique_size = ell)

    clique.set <- clique_split(clique.set, min_clique_size = ell)

    print(paste("Number of Cliques of size,",ell,":", length(clique.set)))

    graph.stats[sim,1] <- n # Size
    graph.stats[sim,2] <- sum(A.sim)/(length(A.sim) - n) # edge density
    graph.stats[sim,3] <- max(d.sim) # maximum degree
    graph.stats[sim,4] <- mean(d.sim) # mean degree
    graph.stats[sim,5] <- length(clique.set) # number of cliques found
    graph.stats[sim,6] <- max.clique.size # largest clique

    g <- igraph::graph_from_adjacency_matrix(A.sim, mode = "undirected")
    graph.stats[sim,7] <- mean(igraph::centr_eigen(g)$vector) # average eigenvector centrality


    clique.set <- clique.set
    if(length(clique.set) > max.num.cliques ){
      clique.set <- clique.set[1:max.num.cliques]
    }

    D.hat <- lolaR::EstimateD(A.sim, clique.set, max.iter = max.iter.estimate)

    # Select the reference set
    reference.set <- lolaR::selectReference(D.hat,
                                            J = num.midpoints,
                                            tri.const = tri.const,
                                            d.yz.min = d.yz.min,
                                            d.yz.max = d.yz.max)

    cc.test <- lolaR::SubSampleConstantCurvatureTest(A.sim,
                                                     clique.set,
                                                     reference.set,
                                                     B = num.subsamples)

    p.val.sub.results[sim,scale.idx] <- cc.test$`p.value`

    #p.values from multiple thresholds of tri.const
    #This uses the previously subsampled values of the matrix
    p.values.mult.thresholds <- SubSampleConstantCurvatureTestMultipleThresholds(D.hat = D.hat,
                                                                                 D.subsample = cc.test$`D.subs`,
                                                                                 reference.set = reference.set,
                                                                                 tri.const.seq = tri.const.seq,
                                                                                 verbose = F)

    for(tri.const.idx in seq(length(tri.const.seq))){

      tri.const.tmp = tri.const.seq[tri.const.idx]
      res.idx <- n.sims*(tri.const.idx - 1) + sim


      reference.set.tmp = reference.set
      y <- reference.set.tmp[[1]]$y
      z <- reference.set.tmp[[1]]$z
      m <- reference.set.tmp[[1]]$m

      x.set.tmp <- filter_indices(D.hat,y,z,m, tri.const = tri.const.tmp)
      reference.set.tmp[[1]]$xset = x.set.tmp
      if(length(x.set.tmp) > 0 ){
        estimates <- lolaR::EstimateKappaSet(D.hat,
                                             reference.set.tmp[[1]][["y"]],
                                             reference.set.tmp[[1]][["z"]],
                                             reference.set.tmp[[1]][["m"]],
                                             reference.set.tmp[[1]][["xset"]])

        kappa.med <- median(estimates, na.rm = T)
      } else {
        kappa.med <- NA
      }

      #best.estimates <- norm.curve.test[[tri.const.idx]]$estimates[norm.curve.test[[tri.const.idx]]$estimates[,1] == 1, 2]

      p.val <- p.values.mult.thresholds$p.value[tri.const.idx]


      kappa.ests.results[res.idx,scale.idx] = kappa.med
      p.val.results[res.idx,scale.idx] = p.val

      # if(!is.null(norm.curve.test[[tri.const.idx]]$estimates)){
      #   kappa.med <- median(best.estimates, na.rm = T)
      #   kappa.ests.results[res.idx,scale.idx] = kappa.med
      #   p.val.results[res.idx,scale.idx] = p.val
      #   #normalized.p.val.results[res.idx,scale.idx] = norm.p.val
      # }
      kappa.ests.results[res.idx,length(scale.set) + 1] = tri.const.tmp
      p.val.results[res.idx,length(scale.set) + 1] = tri.const.tmp
      #normalized.p.val.results[res.idx,length(scale.set) + 1] = tri.const.tmp
    }
  }
  #If tracking the graph statistics
  if(write.graph.stats.files){
    file.graph.stats <- paste0("results/graph_stats_kappa_",kappa,"_scale_",round(scale,1),"_block_",block,".csv")
    write.csv(graph.stats, file = file.graph.stats)
  }

}

if(write.files){
  file.kappa.ests <- paste0("results/estimates_kappa_",kappa,"_block_",block,".csv")
  # file.sl.kappa.ests <- paste0("results/sl_estimates_kappa_",kappa,"_block_",block,".csv")
  file.p.vals <- paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
  # file.norm.p.vals <- paste0("results/norm_p_vals_kappa_",kappa,"_block_",block,".csv")
  file.p.vals.sub <- paste0("results/p_vals_subsample_kappa_",kappa,"_block_",block,".csv")

  write.csv(kappa.ests.results, file = file.kappa.ests)
  # write.csv(sl.kappa.est.results, file = file.sl.kappa.ests)
  write.csv(p.val.results, file = file.p.vals)
  # write.csv(normalized.p.val.results, file = file.norm.p.vals)
  write.csv(p.val.sub.results, file = file.p.vals.sub)
}


time.2 <- Sys.time()

print(paste("Time Difference:", round(time.2 - time.1,3)))



