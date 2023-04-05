
###TODO: Delete File 


library(Matrix)
library(CVXR)
library(Rmosek)
library(MASS)
source("00_functions.R")
source("clique_finder.R")


tri.const <- 1.6
tri.const.seq = 1.6

tri.const <- 1.4
tri.const.seq = 1.4


d.yz.min = 1
d.yz.max = 6



files <- c("data/PhysicsColab/ca-AstroPh.txt",   ##Astrophysics
           "data/PhysicsColab/ca-CondMat.txt",   ##Condensed Matter Physics
           "data/PhysicsColab/ca-GrQc.txt",      ##General Relativity
           "data/PhysicsColab/ca-HepPh.txt",     ##High Energy Particle Physics
           "data/PhysicsColab/ca-HepTh.txt")     ##High Energy Particle Physics (Theory)
#file <- "data/git_web_ml/musae_git_edges.csv"
# edges

dataset.names <- c("Astro",
                   "CondMat",
                   "GR",
                   "HEPApp",
                   "HEPTheory")
K <- length(files)
ell.seq <- c(19,12,6,12,6)
# 

num.cliques.ell <- c()
clique.num <- c()
p.const.curve <- c()
curve.ests <- c()
estimates.set <- list()
midpoint.quality <- c()



k = 1 ## Which network


file = files[k]
ell = ell.seq[k]

d.yz.min = 1
d.yz.max = max(log(ell^2/6), log(ell),3)

E <- read.table(file)
n <- max(E) + 1

i = E[,1] + 1 # data is zero indexed 
j = E[,2] + 1
x = rep(1,nrow(E))
A <- sparseMatrix(i = i, j = j, x = x, dims = c(n,n))
if(!isSymmetric(A)){
  A <- A + t(A)
}

diag(A) = 0 # Some erroneous collaborations listed with self


# Remove Isolated Nodes
iso.nodes <- which(colSums(A) == 0)
no.iso.nodes <- which(colSums(A) > 0)
A <- A[no.iso.nodes,no.iso.nodes]
nn <- nrow(A)


# 
g.large <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")

clique.set <- igraph::maximal.cliques(g.large, min = ell)
clique.set <- clique_partition(clique.set)
clique.set <- clique_split(clique.set, min_clique_size = ell)
clique.set <- clique_set_connected_subgraph(A,clique.set) #reduces to largest connected component
print(paste("Number of Cliques of Size:", ell, ":::", length(clique.set)))
num.cliques.ell[k] <- length(clique.set) 
clique.num[k] <- clique_num(g.large)
plot_cliques(A,clique.set)

estimates <- estimate_curvature(A,clique.set, tri.const = tri.const, d.yz.min = d.yz.min, d.yz.max = d.yz.max)
D.hat <- estimates$D

dat.name <- dataset.names[k]
# save(E, file = paste0(dat.name,"_citation_edges.rda"))
# save(A, file = paste0(dat.name,"_adjacency_matrix.rda"))
# save(D.hat, file = paste0(dat.name,"_distance_estimate.rda"))

plot_cliques_curvature(A,estimates,clique.set, num.midpoints = 3, tri.const = tri.const)
# D0 <- estimates$D
# estimates <- estimate_curvature(A,clique.set, tri.const = 1.4, 
#                                 d.yz.min = 1, d.yz.max = log(ell), 
#                                 no.refit = T,D0)
#TODO: fix the discrepancy between above and below
norm.curve.test <- normalized_constant_curvature_test_seq(estimates, num.midpoints = 3,
                                                          tri.const.seq = tri.const.seq,
                                                          curve.scale = curve.scale)

curve.ests[k] <- estimates$kappa.med

p.const.curve[k] <- norm.curve.test[[1]]$p.value
estimates.set[[k]] <- estimates
mp.quality <- ((estimates$midpoints[1,4]^2 + estimates$midpoints[1,5]^2)/estimates$midpoints[1,6]^2)

midpoint.quality[k] <- mp.quality
print(p.const.curve)
print(curve.ests)
tmp <- norm.curve.test[[1]]$estimates

multi.loc.eds <- tmp %>% group_by(loc) %>% 
  summarise(meds = median(est))
print(multi.loc.eds)


J = 5  # Number of midpoints

tri.const <- 1.4
midpoints <- estimates$midpoints


lower.bounds <- c()
upper.bounds <- c()

lower.bounds.med <- c()
upper.bounds.med <- c()

locations <- c()
for(j in seq(J)){
  y = midpoints[j,1]
  z = midpoints[j,2]
  m = midpoints[j,3]
  
  x.set <- filter_indices_2(D.hat, y,
                            z,m, 
                            tri.const = tri.const)
  bounds <- estimateBounds(D.hat,y,z,m,x.set)
  
  lower.bounds <- c(lower.bounds, bounds$lower.bounds)
  upper.bounds <- c(upper.bounds, bounds$upper.bounds)
  
  lower.bounds.med <- c(lower.bounds.med, median(bounds$lower.bounds))
  upper.bounds.med <- c(upper.bounds.med, median(bounds$upper.bounds))
  
  locations <- c(locations, rep(j, length(bounds$upper.bounds)))
}

idx <- seq(length(locations))
dat <- data.frame("idx" = idx, 
                  "upper" = upper.bounds, 
                  "lower" = lower.bounds, 
                  "locations" = locations)
dat$locations <- as.factor(dat$locations)

dat.med <- data.frame("idx" = seq(J), 
                      "upper" = upper.bounds.med, 
                      "lower" = lower.bounds.med, 
                      "locations" = seq(J))
dat.med$locations <- as.factor(dat.med$locations)

ggplot(data = dat, aes(x = idx, y = upper, shape = locations, col = "upper")) +
  geom_errorbar(data = dat, aes(x = idx, ymin = lower, ymax = upper, col = locations)) + 
  ylim(-10,3) + 
  geom_hline(yintercept = kappa) + 
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 12") + 
  ylab("Curvature Estimate")

ggplot(data = dat.med, aes(x = idx, y = upper, shape = locations)) +
  geom_errorbar(data = dat.med, aes(x = idx, ymin = lower, ymax = upper, col = locations)) + 
  ylim(-300,3) + 
  geom_hline(yintercept = kappa) + 
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 12") + 
  ylab("Curvature Estimate")






##### Simulated Datasets 


scale.idx <- 2
kappa.idx <- 3


n.sims = 10
sim.idx <- 1:n.sims


kappa.set <- c(-2,-1,-0.5,0,0.5,1)
scale.set <- c(1/sqrt(2),1,2,4)


mu = -3
sd = 3

# midpoint objective constant > 0.5
#Cm <- 0.75

sim.avg.variance <- 0.25**2
p = 3
num.midpoints = 3

res = 1

c1 = .5
c2 = 2
c3 = .25


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



print(paste("Max cliques size:",max.clique.size))


# making sure at least ~ 50 cliques are found. 
# when scale > 8 we have to use an approximate clique search 
clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels, 
                                min_clique_size = ell)
clique.set <- clique_split(clique.set, min_clique_size = ell)

print(paste("Number of Cliques of size,",ell,":", length(clique.set)))

if(length(clique.set) > 60 ){
  clique.set <- clique.set[1:60]
}

estimates = estimate_curvature(A.sim, clique.set, no.refit = F, d.yz.min = d.yz.min,d.yz.max = d.yz.max)
D.hat <- estimates$D

J = 5  # Number of midpoints

tri.const <- 1.4
midpoints <- estimates$midpoints


lower.bounds <- c()
upper.bounds <- c()

lower.bounds.med <- c()
upper.bounds.med <- c()

locations <- c()
for(j in seq(J)){
  y = midpoints[j,1]
  z = midpoints[j,2]
  m = midpoints[j,3]
  
  x.set <- filter_indices_2(D.hat, y,
                            z,m, 
                            tri.const = tri.const)
  bounds <- estimateBounds(D.hat,y,z,m,x.set)
  
  lower.bounds <- c(lower.bounds, bounds$lower.bounds)
  upper.bounds <- c(upper.bounds, bounds$upper.bounds)
  
  lower.bounds.med <- c(lower.bounds.med, median(bounds$lower.bounds))
  upper.bounds.med <- c(upper.bounds.med, median(bounds$upper.bounds))
  
  locations <- c(locations, rep(j, length(bounds$upper.bounds)))
}

idx <- seq(length(locations))
dat <- data.frame("idx" = idx, 
                  "upper" = upper.bounds, 
                  "lower" = lower.bounds, 
                  "locations" = locations)
dat$locations <- as.factor(dat$locations)

dat.med <- data.frame("idx" = seq(J), 
                      "upper" = upper.bounds.med, 
                      "lower" = lower.bounds.med, 
                      "locations" = seq(J))
dat.med$locations <- as.factor(dat.med$locations)

ggplot(data = dat, aes(x = idx, y = upper, shape = locations, col = "upper")) +
  geom_errorbar(data = dat, aes(x = idx, ymin = lower, ymax = upper, col = locations)) + 
  ylim(-10,3) + 
  geom_hline(yintercept = kappa) + 
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 12") + 
  ylab("Curvature Estimate")

ggplot(data = dat.med, aes(x = idx, y = upper, shape = locations)) +
  geom_errorbar(data = dat.med, aes(x = idx, ymin = lower, ymax = upper, col = locations)) + 
  ylim(-20,3) + 
  geom_hline(yintercept = kappa) + 
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 6") + 
  ylab("Curvature Estimate")



ggplot(data = dat.med, aes(x = idx, y = upper - lower)) +
  geom_line() + 
  ylim(0,3) +
  ggtitle("Bias of Estimates Plot") + 
  ylab("Curvature Estimate")





