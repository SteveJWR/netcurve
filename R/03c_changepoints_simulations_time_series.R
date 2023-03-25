

rm(list = ls())

library(lolaR)
source("R/00_functions.R")
source("R/clique_finder.R")
source("R/SubsampleConstantCurvatureTest.R")

rm(filter_indices)
rm(optimal_midpoint_search)

#install.packages("remotes")
#library(remotes)
#install_github("guillemr/robust-fpop")

library(robseg) # robust segmentation 


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

#n.sims = 200 
n.sims = 2



# Simulation Parameters
mu = -3
sd = 3
approximate.variance <- 0.25**2
centers.radius = 2.5
p = 3
rho = .5


# Changes in the curvature
Time.steps = 50
t.change1 = 15
t.change2 = 38
kappa1 = 1
kappa2 = 0.25
kappa3 = 1.3






scale <- 1

scale.set <- c(1/sqrt(2),1,2,4)
# scale 4 takes too long for this simulation
# consistency is already basically there earlier
scale.set <- c(1/2,1/sqrt(2),1,sqrt(2),2) 

ell.set <- round(8 + 4*log2(scale.set))

# Estimator tuning parameters
num.midpoints = 3
tri.const = 1.4
d.yz.min <- 1.5

# the problem is that this seems to be a hard geometry to search for cliques
#mu <- -4  # mu for spherical ish geometry 
results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
colnames(results) = paste0("CliqueSize_", ell.set)

time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){
  
  scale <- scale.set[scale.idx]
  
  n <- round(5000*scale)
  
  n.centers <- round(100*sqrt(scale))
  
  ell = round(8 + 4*log2(scale)) # min clique-size 
  
  d.yz.min = 1.0
  #d.yz.max = -log(10/ell^2)
  d.yz.max = max(log(ell),log(ell^2/10)) 
  for(sim in seq(n.sims)){
    
    PI <- as.numeric(rdirichlet(1, rep(2,n.centers)))
    #PI <- rep(1/n.centers, n.centers) 
    
    cluster.model.variance = rgamma(n.centers, shape = approximate.variance)
    #cluster.model.variance <- rep(approximate.variance, n.centers) 
    
    lpcm_T <- lpcm_spherical_rand_walk_dist(n,n.centers, p, 
                                            kappa = kappa1, 
                                            approximate.variance = approximate.variance, 
                                            PI = PI, 
                                            centers.radius = centers.radius, 
                                            Time.steps = Time.steps, rho = rho)
    
    
    
    kappa.seq <- list()
    
    nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
    nu.vec <- nu.vec*(nu.vec < 0 )
    nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 
    
    
    x <- c()
    kappa.seq <- c()
    kappa.true.seq <- c()
    for(time in seq(Time.steps)){
      cat(paste0(time,"/",Time.steps), end = "\r")
      if(time < t.change1){
        kappa.t = kappa1
      } else if(time < t.change2) {
        kappa.t = kappa2
      } else {
        kappa.t = kappa3
      }
      
      Dt <- pos_to_dist(lpcm_T$Z.set[[time]], kappa.t)
      
      
      At <- sim_ls_network_fast_2(nu.vec,lpcm_T$Z.set[[time]], kappa.t)
      
      
      #sim_ls_network_fast_2(nu.vec,Dt,kappa.t)
      clique.set.t <- guided_clique_set(At,lpcm_T$cluster_labels, 
                                        min_clique_size = ell)
      # idx <- lpcm_T$cluster_labels == 2
      # At[idx,idx]
      # print(mean(At))
      #length(clique.set.t)
      #print(mean(At))
      # skipping sections of few available points 
      if(length(clique.set.t) < 10){
        next 
      }
      estimates.t <- estimate_curvature(At, clique.set.t, tri.const = tri.const)
      
      if(length(estimates.t$kappas) > 0 ){
        kappa.seq <- c(kappa.seq, estimates.t$kappas)
        x <- c(x, rep(time, length(estimates.t$kappas)))
        kappa.true.seq <- c(kappa.true.seq, rep(kappa.t, length(estimates.t$kappas)))
      }
    }
    
    kappa.seq <- as.numeric(kappa.seq)
    y <- scale_curvature(kappa.seq, curve.scale)
    y.true <- scale_curvature(kappa.true.seq, curve.scale)
     
    idx.clean <- !is.na(y)
    x.clean <- x[idx.clean]
    y.clean <- y[idx.clean]
    y.true.clean <- y.true[idx.clean]
    dat <- data.frame("x" = x, "y" = y)
    
    est.sd <- mad(diff(y.clean)/sqrt(2))
    ## run dynamic programming
    res.l1 <- Rob_seg.std(x = y.clean,  
                          loss = "L1", 
                          lambda = 4*est.sd*log(length(y.clean)))
    y.hat.smt <- res.l1$smt
    loss <- mean((y.true.clean - y.hat.smt)^2)
    

    plot(y.true.clean, ylim = c(0,1.5))
    lines(y.hat.smt, col = "red")
    results[sim,scale.idx] = loss
    
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
    
  }
  
  ### number of measurements change each time. 
  # here we should see what happens 
  
  y.clean <- y[!is.na(y)]
  x.clean <- x[!is.na(y)]
  
  cpt.true <- c(min(which(x.clean == t.change1 )), 
                min(which(x.clean == t.change2 )))
  
  est.sd <- mad(diff(y.clean)/sqrt(2))
  ## run dynamic programming
  res.l1 <- Rob_seg.std(x = y.clean,  
                        loss = "L1", 
                        lambda = (4)*est.sd*log(length(y.clean)))
  ## estimated changepoints 
  cpt <- res.l1$t.est[-length(res.l1$t.est)]
  
}

csv.file <-  paste0("results/changepoint_results","_block_",id,".csv")
write.csv(results, file = csv.file)

time.2 <- Sys.time()

print(paste("Time Difference:", time.2 - time.1))
