
#load("3b_coverage.RData")
#save.image(file = "3b_coverage.RData")
rm(list = ls())

library(lolaR)
source("R/00_functions.R")
# for block_diag function
library(Matrix)

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

set.seed(id)
n.sims = 10

sim.idx <- 1:n.sims

# Simulation Parameters
mu = -3
sd = 3

scale.set <- c(1/sqrt(2),1,2) #c(1/sqrt(2),1,2,4)
ell.set <- round(8 + 4*log2(scale.set))


# Curvature of different viewed latent spaces
kappa1 = 0.8
kappa2 = 0.5

p = 3

approximate.variance <- 0.25**2
centers.radius = 2.5



# storage matrices for the results
p.val.results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
colnames(p.val.results) = paste0("CliqueSize_", ell.set)


# Method Tuning Parameters
d.yz.min <- 1
tri.const = 1.5 # constant for the filtering term
# Method Tuning Parameters
#num.midpoints = 3
tri.const = 1.5
max.num.cliques.per.view = 18 # number in each view.
min.num.cliques.per.view = 6
num.subsamples = 250
max.iter.estimate = 3
d.yz.min = 1


time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){

  scale <- scale.set[scale.idx]

  n <- round(5000*scale)

  n.centers <- round(100*sqrt(scale))

  ell = round(8 + 4*log2(scale)) # min clique-size



  #rule of thumb for number of cliques
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

    if(length(clique.set1) > min.num.cliques.per.view & length(clique.set2) > min.num.cliques.per.view){
      # minimizing numbers of cliques for computational speed
      if(length(clique.set1) > max.num.cliques.per.view ){
        clique.set1 <- clique.set1[1:max.num.cliques.per.view]
      }

      if(length(clique.set2) > max.num.cliques.per.view ){
        clique.set2 <- clique.set2[1:max.num.cliques.per.view]
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
                                        J = 1,
                                        tri.const = tri.const,
                                        d.yz.min = d.yz.min,
                                        d.yz.max = d.yz.max)
      midset1 <- c(reference.set1[[1]]$y,reference.set1[[1]]$z, reference.set1[[1]]$m)

      reference.set2 <- selectReference(D.hat2,
                                        J = 1,
                                        tri.const = tri.const,
                                        d.yz.min = d.yz.min,
                                        d.yz.max = d.yz.max)
      midset2 <- c(reference.set2[[1]]$y,reference.set2[[1]]$z, reference.set2[[1]]$m)


      # picking a reference set point from each block.
      reference.set = list(reference.set1[[1]], reference.set2[[1]])
      reference.set[[2]]$y = reference.set[[2]]$y + nrow(D.hat1)
      reference.set[[2]]$z = reference.set[[2]]$z + nrow(D.hat1)
      reference.set[[2]]$m = reference.set[[2]]$m + nrow(D.hat1)
      reference.set[[2]]$x.set = reference.set[[2]]$x.set + nrow(D.hat1)



      cc.test <- SubSampleConstantCurvatureTest(A.sim,
                                                clique.set,
                                                reference.set,
                                                B = num.subsamples)


      if(!is.null(cc.test$p.value)){
        p.val.results[sim,scale.idx] = cc.test$p.value
      }
    }
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
}

if(write.files){
  csv.file <- paste0("results/multiview_results","_block_",id,".csv")
  write.csv(p.val.results, file = csv.file)
}
time.2 <- Sys.time()

print(paste("Time Difference:", time.2 - time.1))


