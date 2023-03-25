
#load("3b_coverage.RData")
#save.image(file = "3b_coverage.RData")
rm(list = ls())

library(lolaR)
source("R/00_functions.R")
source("R/clique_finder.R")
source("R/SubsampleConstantCurvatureTest.R")

rm(filter_indices)
rm(optimal_midpoint_search)
# for block-diag 
library(Matrix)

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

sim.idx <- 1:n.sims

# Simulation Parameters
mu = -3
sd = 3

scale <- 1
scale.set <- c(1/sqrt(2),1,2,4)

ell.set <- round(8 + 4*log2(scale.set))


# Curvature of different viewed latent spaces 
kappa1 = 0.8
kappa2 = 0.5

p = 3

sim.avg.variance <- 0.25**2
centers.radius = 2.5




# Model tuning parameters
num.midpoints = 3
#tri.const.seq = (seq(0, 1, length.out = 21)) + 1
tri.const = 1.6
d.yz.min <- 1.5


# rule of thumb for number of connections




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
    
    A.sim1 <- sim_ls_network_fast_2(nu.vec, Z, kappa1)
    A.sim2 <- sim_ls_network_fast_2(nu.vec, Z, kappa2)
    
    clique.set1 <- guided_clique_set(A.sim1,lpcm$cluster_labels, 
                                     min_clique_size = ell)
    clique.set2 <- guided_clique_set(A.sim2,lpcm$cluster_labels, 
                                     min_clique_size = ell)
    
    clique.set1 <- clique_split(clique.set1, min_clique_size = ell)
    clique.set2 <- clique_split(clique.set2, min_clique_size = ell)
    
    print(paste("Number of Cliques of size,",ell,":", length(clique.set1)))
    print(paste("Number of Cliques of size,",ell,":", length(clique.set2)))
    
    if(length(clique.set1) > 6 & length(clique.set2) > 6){
      # minimizing numbers of cliques for computational speed 
      if(length(clique.set1) > 30 ){
        clique.set1 <- clique.set1[1:30]
      }
      
      if(length(clique.set2) > 30 ){
        clique.set2 <- clique.set2[1:30]
      }
      
      clique.set <- clique.set1
      for(k in 1:length(clique.set2)){
        clique.set[[k + length(clique.set1)]] <- clique.set2[[k]] + nrow(A.sim1)
      }
      # create a larger adjacency matrix 
      A.sim = bdiag(A.sim1, A.sim2)
      
      D.hat1 <- lolaR::EstimateD(A.sim1, clique.set1)
      D.hat2 <- lolaR::EstimateD(A.sim2, clique.set2)
      
      reference.set1 <- selectReference(D.hat1,
                                        J = num.midpoints, 
                                        tri.const = tri.const)
      midset1 <- c(reference.set1[[1]]$y,reference.set1[[1]]$z, reference.set1[[1]]$m) 
      D.hat1[midset1, midset1]
      reference.set2 <- selectReference(D.hat2,
                                        J = num.midpoints, 
                                        tri.const = tri.const)
      midset2 <- c(reference.set2[[1]]$y,reference.set2[[1]]$z, reference.set2[[1]]$m) 
      D.hat2[midset2, midset2]
      
      # picking a reference set point from each block. 
      reference.set = list(reference.set1[[1]], reference.set2[[1]])
      reference.set[[2]]$y = reference.set[[2]]$y + nrow(D.hat1)
      reference.set[[2]]$z = reference.set[[2]]$z + nrow(D.hat1)
      reference.set[[2]]$m = reference.set[[2]]$m + nrow(D.hat1)
      reference.set[[2]]$x.set = reference.set[[2]]$x.set + nrow(D.hat1)
      
      
      
      cc.test <- SubSampleConstantCurvatureTest(A.sim,
                                                clique.set,
                                                reference.set, 
                                                B = 200)

      
      if(!is.null(cc.test$p.value)){
        p.val.results[sim,scale.idx] = cc.test$p.value
      }
    }
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
}

csv.file <- paste0("results/multiview_results","_block_",id,".csv")
write.csv(p.val.results, file = csv.file)

time.2 <- Sys.time()

print(paste("Time Difference:", time.2 - time.1))


