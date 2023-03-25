
#load("3b_coverage.RData")
#save.image(file = "3b_coverage.RData")
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
    
    if(length(clique.set1) > 60 ){
      clique.set1 <- clique.set1[1:60]
    }
    estimates1 = estimate_curvature(A.sim1, clique.set1)
    
    
    #clique.set2 <- cluster_clique_search(A.sim2,res1 = 1,res2 = 1, min_clique_size = ell)
    #clique.set <- clique_split(clique.set, min_clique_size = ell)
    # }
    
    
    if(length(clique.set2) > 60 ){
      clique.set2 <- clique.set2[1:60]
    }
    estimates2 = estimate_curvature(A.sim2, clique.set2)

    
    kappa.set.1 <- estimates1$kappas
    kappa.set.2 <- estimates2$kappas
    
    
    kappa.set.1 <- kappa.set.1[!is.na(kappa.set.1)]
    kappa.set.2 <- kappa.set.2[!is.na(kappa.set.2)]
    
    index <- c(rep(1,length(kappa.set.1)),
               rep(2,length(kappa.set.2)))
    kappa.vec <- c(kappa.set.1,
                   kappa.set.2)
    
    y = scale_curvature(kappa.vec,curve.scale)
    x = index
    
    # good, this test basically rejects every time
    
    test.dat <- data.frame("loc" = x, "est" = kappa.vec)
    test <- kruskal.test(est ~ loc, data = test.dat) 
    p.val.results[sim,scale.idx] = test$p.value
    cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
  }
  
  
  
  # Example plot of a change point problem: 
  # save one for each sample size 
  
  ### 
  
  ## simple plotting of changes and smoothed profile
  # plot(y.clean, pch=20, col="black")
  # lines(res.l1$smt, col="red", lwd=2)
  # lines(y.true.clean, lty=4, col="green", lwd=2)
  # abline(v=cpt, lty=2, col="red")
  # abline(v=cpt.true, lty=3, col="green")
  
  plot.file <- paste0("plots/Multiview_CliqueSize_",ell,"_block_",id,".png")
  plot.title <- paste0("Multiview Curvature Clique Size: ",ell)
  
  
  labels <- c(rep(0,length(kappa.set.1)),rep(1,length(kappa.set.2)))
  bp.dat <- data.frame(matrix(c(labels, 
                                scale_curvature(kappa.set.1,curve.scale),
                                scale_curvature(kappa.set.2,curve.scale)), ncol = 2))
  colnames(bp.dat) <- c("group", "curvature")
  
  
  
  med.1 <- median(scale_curvature(kappa.set.1,curve.scale),na.rm = T)
  med.2 <- median(scale_curvature(kappa.set.2,curve.scale),na.rm = T)
  
  kappa.vec <- c(kappa.set.1,
                 kappa.set.2)
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








# 
# 
# 
# 
# scale <- 1
# n <- round(5000*sqrt(scale))
# # n <- round(1000*sqrt(scale))
# 
# n.centers <- round(100*sqrt(scale))
# 
# # n.centers1 <- round(20*sqrt(scale))
# # n.centers2 <- round(20*sqrt(scale))
# # probability of ending up on which side of the sphere 
# 
# mu = -3
# sd = 3
# 
# kappa1 = 1
# kappa2 = 0.5
# 
# approximate.variance <- 0.25**2
# d.yz.min <- 1.5
# centers.radius = 2.2
# 
# 
# p = 3
# ell = round(8 + 4*log2(scale)) # min clique-size 
# c1 = .5
# c2 = 2
# c3 = .25
# d.yz.min = 1.5
# d.yz.max = -log(10/ell^2)
# d.yz.max <- max(log(ell),log(ell^2/10)) # rule of thumb for number of connections
# set.seed(1)
# 
# 
# 
# n.sims <- 200
# curve.scale <- 10 
# p.set <- rep(NA,n.sims)
# 
# # the problem is that this seems to be a hard geometry to search for cliques
# #mu <- -4  # mu for spherical ish geometry 
# for(sim in seq(n.sims)){
#   
#   PI <- as.numeric(rdirichlet(1, rep(2,n.centers)))
#   
#   cluster.model.variance = rgamma(n.centers, shape = approximate.variance)
#   
#   lpcm <- latent_position_cluster_model(n,n.centers, p, 
#                                         cluster.radius, 
#                                         kappa1, 
#                                         cluster.model.variance, 
#                                         PI = PI)
#                                      
#   Z <- lpcm$Z
#   
#   nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale) 
#   nu.vec <- nu.vec*(nu.vec < 0 )
#   nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec 
#   
#   D1 <- pos_to_dist(Z, kappa1)
#   D2 <- pos_to_dist(Z, kappa2)
#   
#   A.sim1 <- sim_ls_network(nu.vec, D1)
#   A.sim2 <- sim_ls_network(nu.vec, D2)
#   
#   # P.sim <- exp(outer(nu.vec, nu.vec, "+") -D)
#   # A.sim <- 1*(P.sim >= 0.5)
#   print(mean(A.sim1))
#   print(mean(A.sim2))
# 
#   # making sure at least ~ 50 cliques are found. 
#   # when scale > 8 we have to use an approximate clique search 
#   # if(1 <= scale & scale < 4 & kappa1 <= 0 & kappa2 <= 0){
#   #   
#   #   clique.set <- clique_finder_2(A.sim, ell, 2500)
#   # 
#   # } else {
#   # trim the network 
#   
#   #A.trim[d.sim <= quantile(d.sim, 0.2), ] = 0
#   #A.trim[,d.sim <= quantile(d.sim, 0.2)] = 0
#   # the trim is meant to speed up the clique search 
#   
#   #g.large = igraph::graph_from_adjacency_matrix(A.trim, mode = "undirected")
#   #full.cliques <- igraph::maximal.cliques(g.large, min = ell)
#   clique.set1 <- clique_finder(A.sim1, min_clique_size = ell,num_cliques_stop = 70, num_iters = 4000)
#   clique.set1 <- clique_split(clique.set1, min_clique_size = ell)
#   # }
#   clique.set2 <- clique_finder(A.sim2, min_clique_size = ell,num_cliques_stop = 70, num_iters = 4000)
#   clique.set2 <- clique_split(clique.set2, min_clique_size = ell)
#   
#   print(paste("Number of Cliques of size,",ell,":", length(clique.set1)))
#   print(paste("Number of Cliques of size,",ell,":", length(clique.set2)))
#   
#   nu.hats1 <- estimate_nus(A.sim1,clique.set1)
#   nu.hats2 <- estimate_nus(A.sim2,clique.set2)
#   
#   
#   # max.deg.idx <- which(colSums(A.sim) == max(colSums(A.sim)))
#   # ego.idx <- which(A.sim[max.deg.idx[1],] == 1)
#   # A.ego <- A.sim[ego.idx,ego.idx]
#   # g.ego <- igraph::graph_from_adjacency_matrix(A.ego, mode = "undirected")
#   # ego.cliques <- igraph::max_cliques(g.ego, min = ell)
#   # ego.largest.cliques <- largest_cliques(g.ego)
#   
#   K <- length(clique.set1)
#   clique.idx <- c()
#   subset.ids <- c()
#   fixed.effect.vec <- c()
#   for(k in seq(K)){
#     #cat(paste0(k,"  ", length(clique.set[[k]])), end = "\n")
#     clique.idx <- c(clique.idx,clique.set1[[k]])
#     subset.ids <- c(subset.ids, rep(k,length(clique.set1[[k]])))
#     fixed.effect.vec <- c(fixed.effect.vec,nu.hats1[[k]])
#   }
#   
#   A.sub <- A.sim1[clique.idx,clique.idx]
#   diag(A.sub) = 0 
#   
#   D0 = init_D0(A.sub,subset.ids,fixed.effect.vec)
#   D0 <- 0.8*D0
#   # this seems to be the fastest LCQP that I can use. 
#   D1 <- estimate_D_restricted_fast(A.sub,subset.ids,fixed.effect.vec,
#                                   D0,thresh = 10**(-3),
#                                   max.iter = 50, solver = "MOSEK",
#                                   verbose = T)
#   
#   
#   K <- length(clique.set2)
#   clique.idx <- c()
#   subset.ids <- c()
#   fixed.effect.vec <- c()
#   for(k in seq(K)){
#     #cat(paste0(k,"  ", length(clique.set[[k]])), end = "\n")
#     clique.idx <- c(clique.idx,clique.set2[[k]])
#     subset.ids <- c(subset.ids, rep(k,length(clique.set2[[k]])))
#     fixed.effect.vec <- c(fixed.effect.vec,nu.hats2[[k]])
#   }
#   
#   A.sub <- A.sim2[clique.idx,clique.idx]
#   diag(A.sub) = 0 
#   
#   D0 = init_D0(A.sub,subset.ids,fixed.effect.vec)
#   D0 <- 0.8*D0
#   # this seems to be the fastest LCQP that I can use. 
#   D2 <- estimate_D_restricted_fast(A.sub,subset.ids,fixed.effect.vec,
#                                    D0,thresh = 10**(-3),
#                                    max.iter = 50, solver = "MOSEK",
#                                    verbose = T)
#   
#   
#   mid.search1 <- optimal_midpoint_search(D1,top.k = 10, 
#                                         d.yz.min = d.yz.min, 
#                                         d.yz.max = d.yz.max)
#   
#   mid.search2 <- optimal_midpoint_search(D2,top.k = 10, 
#                                          d.yz.min = d.yz.min, 
#                                          d.yz.max = d.yz.max)
#   
#   print(mid.search1)
#   print(mid.search2)
#   
#   y.opt1 = mid.search1[1,1]
#   z.opt1 = mid.search1[1,2]
#   m.opt1 = mid.search1[1,3]
#   
#   
#   
#   
#   
#   x.set1 <- filter_indices(D1, 
#                           y.opt1,
#                           z.opt1,
#                           m.opt1, 
#                           c1 = c1,c2 = c2,c3 = c3)
#   
#   
#   kappa.set.1 <- estimate_kappa_set(D1,y.opt1,
#                                     z.opt1,
#                                     m.opt1,x.set1)
#   
#   
#   
#   y.opt2 = mid.search2[1,1]
#   z.opt2 = mid.search2[1,2]
#   m.opt2 = mid.search2[1,3]
#   
#   
#   
#   
#   
#   x.set2 <- filter_indices(D2, 
#                            y.opt2,
#                            z.opt2,
#                            m.opt2, 
#                            c1 = c1,c2 = c2,c3 = c3)
#   
#   
#   kappa.set.2 <- estimate_kappa_set(D2,y.opt2,
#                                     z.opt2,
#                                     m.opt2,x.set2)
#   
#   # median might not be enough on its own to detect changes 
#   # in curvature
#   
#   kappa.set.1 <- kappa.set.1[!is.na(kappa.set.1)]
#   kappa.set.2 <- kappa.set.2[!is.na(kappa.set.2)]
#   
#   index <- c(rep(1,length(kappa.set.1)),
#              rep(2,length(kappa.set.2)))
#   
#   kappa.vec <- c(kappa.set.1,
#                  kappa.set.2)
#   
#   y = scale_curvature(kappa.vec,curve.scale)
#   x = index
#   
#   # good, this test basically rejects every time
#   
#   pt <- perm_median_test(y, x)
#   p.set[sim] = pt$`p-value`
#   cat(paste("Simulation,", sim,"/",n.sims), end = "\r")
# }




