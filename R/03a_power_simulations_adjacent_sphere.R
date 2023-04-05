

rm(list = ls())


library(lolaR)
source("R/00_functions.R")

#whether to write the files
write.files = T

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}


if(id %in% c(1,2,15)){
  set.seed(-id)
} else {
  set.seed(id)
}
n.sims = 10

# Simulation Parameters
# probability of ending up on which side of the adjacent spheres
q = .5

# dimensions of the touching spheres
p1 = 2
p2 = 2

mu = -3
sd = 3

# curvature of each sphere
kappa1 = 0.8
kappa2 = 0.5


# For simulations, let's find the cliques based on the cluster membership


approximate.variance <- 0.25**2

centers.radius = 2.5


scale.set <- c(1/sqrt(2),1,2) #c(1/sqrt(2),1,2,4)
ell.set <- round(8 + 4*log2(scale.set)) # clique sizes based on the scale

# rule of thumb for number of connections
num.midpoints = 3
tri.const.seq = (seq(0, 1, length.out = 21)) + 1




# storage matrices for the results
p.val.results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
colnames(p.val.results) = paste0("CliqueSize_", ell.set)

# Method Tuning Parameters
d.yz.min <- 1
tri.const = 1.5 # constant for the filtering term
# Method Tuning Parameters
num.midpoints = 3
tri.const = 1.5
max.num.cliques = 35
num.subsamples = 250
max.iter.estimate = 3
d.yz.min = 1



time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){

  scale <- scale.set[scale.idx]

  n <- round(5000*scale)

  n.centers1 <- round(50*sqrt(scale))
  n.centers2 <- round(50*sqrt(scale))

  ell = round(8 + 4*log2(scale)) # min clique-size


  #rule of thumb for number of cliques
  d.yz.max = max(log(ell),log(ell^2/10))
  for(sim in seq(n.sims)){
    PI1 <- as.numeric(rdirichlet(1, rep(2,n.centers1)))
    PI2 <- as.numeric(rdirichlet(1, rep(2,n.centers2)))

    nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale)
    nu.vec <- nu.vec*(nu.vec < 0 )
    nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec

    # update the 00_functions with a description for the simulations etc.
    lpcm <- connected_spheres_lpcm(n, n.centers1, n.centers2,
                                   p1,p2,PI1, PI2, nu.vec,
                                   q, kappa1,kappa2,
                                   approximate.variance, max.rad = centers.radius)


    A.sim <- lpcm$A

    clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels,
                                    min_clique_size = ell)

    print(paste("Number of Cliques of size,",ell,":", length(clique.set)))
    if(length(clique.set) > 60 ){
      clique.set <- clique.set[1:60]
    }
    # ensuring a minimum number of cliques are present for the
    if(length(clique.set) > 6){

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
      if(!is.null(cc.test$`p.value`)){
        p.val.results[sim,scale.idx] = cc.test$`p.value`
      }
    }

    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
}

if(write.files){
  csv.file <- paste0("results/adjacent_spheres_results","_block_",id,".csv")
  write.csv(p.val.results, file = csv.file)

}

time.2 <- Sys.time()
print(paste("Time Difference:",round(time.2 - time.1,3)))









