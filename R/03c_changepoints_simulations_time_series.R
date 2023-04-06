#install.packages("remotes")
#library(remotes)
#install_github("guillemr/robust-fpop")

rm(list = ls())

library(lolaR)
source("R/00_functions.R")
library(robseg) # robust segmentation package for changepoints

#whether to write the files
write.files = T
plot.cpt = F

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

set.seed(id)

n.sims = 2

# Simulation Parameters
mu = -3
sd = 3
approximate.variance <- 0.25**2
centers.radius = 2.5
p = 3
rho = .5 #drift parameter for latent positions

# scale 4 takes too long for this simulation
# consistency is already basically there earlier
scale.set <- c(1/2,1/sqrt(2),1,sqrt(2),2)
ell.set <- round(8 + 4*log2(scale.set))

# Changes in the curvature
Time.steps = 50
t.change1 = 15
t.change2 = 38

# three curvature values
kappa1 = 1
kappa2 = 0.25
kappa3 = 1.3


# Method Tuning Parameters
tri.const = 1.3 # constant for the filtering term
# Method Tuning Parameters
#num.midpoints = 3
max.num.cliques.per.time = 18 # number in each view.
min.num.cliques.per.time = 6
#num.subsamples = 250
max.iter.estimate = 3
d.yz.min = 1

# Changepoint method parameters
curve.thresh = 10
change.reg = 2

# the problem is that this seems to be a hard geometry to search for cliques
#mu <- -4  # mu for spherical ish geometry
results <- matrix(NA, nrow = n.sims, ncol = length(scale.set))
results.med <- results
colnames(results) = paste0("CliqueSize_", ell.set)

time.1 <- Sys.time()
for(scale.idx in seq(length(scale.set))){

  scale <- scale.set[scale.idx]

  n <- round((3/4)*5000*scale)
  n.centers <- round((3/4)*100*sqrt(scale))
  ell = round(8 + 4*log2(scale)) # min clique-size

  #rule of thumb for number of cliques
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
    kappa.med.seq <- c()
    kappa.med.true.seq <- c()
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

      clique.set.t <- guided_clique_set(At,lpcm_T$cluster_labels,
                                        min_clique_size = ell)

      # skipping sections of few available points
      if(length(clique.set.t) > max.num.cliques.per.time){
        clique.set.t <- clique.set.t[1:max.num.cliques.per.time]
      }

      if(length(clique.set.t) < min.num.cliques.per.time){
        next
      }
      D.hat.t <- lolaR::EstimateD(At, clique.set.t, max.iter = max.iter.estimate)
      #print(paste0("Time: ", time, "/",Time.steps ))
      estimates.t <- lolaR::EstimateCurvature(D.hat.t,
                                              J = 1,
                                              tri.const = tri.const,
                                              d.yz.min = d.yz.min,
                                              d.yz.max = d.yz.max)

      if(length(estimates.t$kappas) > 0 ){
        kappa.seq <- c(kappa.seq, estimates.t$kappas)
        x <- c(x, rep(time, length(estimates.t$kappas)))
        kappa.true.seq <- c(kappa.true.seq, rep(kappa.t, length(estimates.t$kappas)))
        kappa.med.seq <- c(kappa.med.seq, median(estimates.t$kappas, na.rm = T))
        kappa.med.true.seq <- c(kappa.med.true.seq, kappa.t)
      }
    }

    kappa.seq <- as.numeric(kappa.seq)
    y <- kappa.med.seq #scale_curvature(kappa.seq, curve.scale)
    y.true <- kappa.med.true.seq #scale_curvature(kappa.true.seq, curve.scale)

    #y[y > curve.thresh] = curve.thresh
    #y[y < -curve.thresh] = -curve.thresh

    idx.clean <- !is.na(y)

    x.clean <- x[idx.clean]
    y.clean <- y[idx.clean]
    y.true.clean <- y.true[idx.clean]
    dat <- data.frame("x" = seq(length(y)), "y" = y)

    est.sd <- mad(diff(y.clean)/sqrt(2))
    ## run dynamic programming
    res.l1 <- Rob_seg.std(x = y.clean,
                          loss = "Outlier",
                          lambda = change.reg*est.sd*log(length(y.clean)),
                          lthreshold=3)
    y.hat.smt <- res.l1$smt
    loss <- mean(abs(y.true.clean - y.hat.smt))
    loss.med <- median(abs(y.true.clean - y.hat.smt))
    results[sim,scale.idx] = loss
    results.med[sim,scale.idx] = loss.med

    if(plot.cpt){
      plot(y.true.clean, ylim = c(-5,5))
      points(y.clean, col = "blue")
      lines(y.hat.smt, col = "red")
    }
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
}

if(write.files){
  csv.file <-  paste0("results/changepoint_results","_block_",id,".csv")
  write.csv(results, file = csv.file)
  csv.file <-  paste0("results/changepoint_results_median","_block_",id,".csv")
  write.csv(results.med, file = csv.file)
}


time.2 <- Sys.time()
print(paste("Time Difference:", time.2 - time.1))


