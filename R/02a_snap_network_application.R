
rm(list = ls())
library(Matrix)
library(CVXR)
library(Rmosek)
library(MASS)
library(lolaR)
source("R/00_functions.R")
source("R/clique_finder.R")


compute.medians = T
tri.const <- 1.2
#tri.const.seq = 1.2


d.yz.min = 1
d.yz.max.overall = 3
max.iter.estimate = 5
num.midpoints = 3
max.cliques = 80 #TODO: Set to 80
num.subsamples = 400 #TODO: Set to 100

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


for(k in seq(K)){
  file = files[k]
  ell = ell.seq[k]

  d.yz.min = 1
  d.yz.max = max(log(ell^2/6), log(ell),d.yz.max.overall)

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

  if(length(clique.set) > max.cliques){
    clique.set <- clique.set[1:max.cliques]
  }

  #estimates <- estimate_curvature(A,clique.set, tri.const = tri.const, d.yz.min = d.yz.min, d.yz.max = d.yz.max)
  #D.hat <- estimates$D
  D.hat <- lolaR::EstimateD(A, clique.set, max.iter = max.iter.estimate)

  # Select the reference set
  reference.set <- lolaR::selectReference(D.hat,
                                          J = num.midpoints,
                                          tri.const = tri.const,
                                          d.yz.min = d.yz.min,
                                          d.yz.max = d.yz.max)

  cc.test <- lolaR::SubSampleConstantCurvatureTest(A,
                                                   clique.set,
                                                   reference.set,
                                                   B = num.subsamples)

  dat.name <- dataset.names[k]
  # save(E, file = paste0(dat.name,"_citation_edges.rda"))
  # save(A, file = paste0(dat.name,"_adjacency_matrix.rda"))
  # save(D.hat, file = paste0(dat.name,"_distance_estimate.rda"))

  if(compute.medians){
    med1 <- median(lolaR::EstimateKappaSet(D.hat,
                                           reference.set[[1]]$y,
                                           reference.set[[1]]$z,
                                           reference.set[[1]]$m,
                                           reference.set[[1]]$xset), na.rm = T)
    med2 <- median(lolaR::EstimateKappaSet(D.hat,
                                           reference.set[[2]]$y,
                                           reference.set[[2]]$z,
                                           reference.set[[2]]$m,
                                           reference.set[[2]]$xset), na.rm = T)
    med3 <- median(lolaR::EstimateKappaSet(D.hat,
                                           reference.set[[3]]$y,
                                           reference.set[[3]]$z,
                                           reference.set[[3]]$m,
                                           reference.set[[3]]$xset), na.rm = T)
  }

  curve.estimate <- EstimateCurvature(D.hat, J = num.midpoints,
                                      tri.const = tri.const,
                                      d.yz.min = d.yz.min,
                                      d.yz.max = d.yz.max)

  plt <- plot_cliques_curvature(A,D.hat,curve.estimate$midpoints,
                                clique.set, num.midpoints = 3, tri.const = tri.const)

  # norm.curve.test <- normalized_constant_curvature_test_seq(estimates, num.midpoints = 3,
  #                                                           tri.const.seq = tri.const.seq,
  #                                                           curve.scale = curve.scale)

  curve.ests[k] <- curve.estimate$kappa.med

  p.const.curve[k] <- cc.test$p.value
  estimates.set[[k]] <- curve.estimate
  mp.quality <- ((curve.estimate$midpoints[1,4]^2 + curve.estimate$midpoints[1,5]^2)/curve.estimate$midpoints[1,6]^2)

  midpoint.quality[k] <- mp.quality
  #print(p.const.curve)
  #print(curve.ests)
  #tmp <- norm.curve.test[[1]]$estimates

  # multi.loc.eds <- tmp %>% group_by(loc) %>%
  #   summarise(meds = median(est))
  # print(multi.loc.eds)

}

print(round(p.const.curve, 3))
print(round(curve.ests, 3))
print(round(midpoint.quality, 6))

# we see that the astrophysics network does not appear to have constant curvature
# we plot these results and test whether curvature is different amongst the groups

#data.frame("loc" = index.vec, "trim.est" = trim.kappa.vec,"est" = kappa.vec)
# curve.test.dat <- matrix(NA, nrow = 0, ncol = 3)
# med.vec <- c()
# curve.scale = 10
# for(k in 1:5){
#   kappas <- estimates.set[[k]]$kappas
#   y <- scale_curvature(kappas, curve.scale)
#   y <- y[!is.na(y)]
#   scl <- mean(abs(y - median(y)))
#   y.norm <- (y - median(y))/scl + median(y)
#   dat.vec <- c(rep(k,length(y)),y.norm,y)
#   curve.test.dat.tmp <- matrix(dat.vec, ncol = 3)
#   curve.test.dat <- rbind(curve.test.dat,curve.test.dat.tmp)
#   med.vec[k] <- median(y)
# }
# colnames(curve.test.dat) <- c("group", "trim.est", "est")
# curve.test.dat <- as.data.frame(curve.test.dat)
# curve.test.dat.sub <- curve.test.dat[curve.test.dat$group != 1, ]
#
# cross.network.curvature.test <- kruskal.test(trim.est ~ group, data = curve.test.dat.sub)
#
# print(paste0("Networks Have Same Curvature p-value: ", round(cross.network.curvature.test$p.value, 4)))
#
# curve.test.dat.sub <- curve.test.dat[!curve.test.dat$group %in% c(1,3), ]
#
# cross.network.curvature.test <- kruskal.test(trim.est ~ group, data = curve.test.dat.sub)
#
# print(paste0("Networks Have Same Curvature p-value: ", round(cross.network.curvature.test$p.value, 4)))
#
#
# ggplot(curve.test.dat, aes(y = est, x = group)) +
#   geom_jitter() +
#   geom_vline(xintercept = 0.5, col = "red") +
#   geom_vline(xintercept = 1.5, col = "red") +
#   geom_vline(xintercept = 2.5, col = "red") +
#   geom_vline(xintercept = 3.5, col = "red") +
#   geom_vline(xintercept = 4.5, col = "red") +
#   geom_vline(xintercept = 5.5, col = "red") +
#   geom_segment(aes(x=0.5,xend=1.5,y=med.vec[1],yend=med.vec[1]), col = "blue") +
#   geom_segment(aes(x=1.5,xend=2.5,y=med.vec[2],yend=med.vec[2]), col = "blue") +
#   geom_segment(aes(x=2.5,xend=3.5,y=med.vec[3],yend=med.vec[3]), col = "blue") +
#   geom_segment(aes(x=3.5,xend=4.5,y=med.vec[4],yend=med.vec[4]), col = "blue") +
#   geom_segment(aes(x=4.5,xend=5.5,y=med.vec[5],yend=med.vec[5]), col = "blue") +
#   ggtitle("Curvature Estimates Across Collaboration Networks") +
#   ylab("Trimmed Curvature") +
#   xlab("Collaboration Network")

# num.cliques.ell
# clique.num










# plot.list <- list()
#
#
# for(k in seq(K)){
#   file = files[k]
#   ell = ell.seq[k]
#
#   d.yz.min = 1
#   d.yz.max = max(log(ell^2/6), log(ell),3)
#
#   E <- read.table(file)
#   n <- max(E) + 1
#
#   i = E[,1] + 1 # data is zero indexed
#   j = E[,2] + 1
#   x = rep(1,nrow(E))
#   A <- sparseMatrix(i = i, j = j, x = x, dims = c(n,n))
#   if(!isSymmetric(A)){
#     A <- A + t(A)
#   }
#
#   diag(A) = 0 # Some erroneous collaborations listed with self
#
#   # TODO: should we be removing these isolated Nodes?
#   # Most seem to be isolated.
#   iso.nodes <- which(colSums(A) == 0)
#   no.iso.nodes <- which(colSums(A) > 0)
#   A <- A[no.iso.nodes,no.iso.nodes]
#   nn <- nrow(A)
#
#
#   #
#   g.large <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
#
#   clique.set <- igraph::maximal.cliques(g.large, min = ell)
#   clique.set <- clique_partition(clique.set)
#   clique.set <- clique_split(clique.set, min_clique_size = ell)
#   clique.set <- clique_set_connected_subgraph(A,clique.set) #reduces to largest connected component
#   plot_cliques(A,clique.set)
#
# }




########### END ###########
# # what is going on with the different numbers of
# # cliques
# ell <- 8
# full.cliques <- igraph::maximal.cliques(g.large, min = ell)
#
# # clique searching doesn't even seem to really be a problem
#
#
# unique.full.cliques <- clique_partition(full.cliques,T)
#
#
# # clique.set <- clique_finder(A,min_clique_size=8,
# #                             num_cliques_stop=500,
# #                             num_iters=20000,
# #                             non_overlap_cliques = T)
#
#
# clique.set <- unique.full.cliques
#
# length(clique.set)
# unique.nodes <- unique(unlist(clique.set))
#
#
# nu.hats <- estimate_nus(A,clique.set)
#
#
#
# labels <- c()
# clique.ids <- c()
# for(lab in 1:length(clique.set)){
#   labels <- c(labels, rep(lab,length(clique.set[[lab]])))
#   clique.ids <- c(clique.ids,clique.set[[lab]])
# }
#
# A.sub <- A[clique.ids,clique.ids]
#
# d.clique <- colSums(A[,clique.ids])
#
# plot(hist(log(d)))
# plot(hist(log(d.clique)))
#
# diag(A.sub) = 0
#
#
# #D.hat <- estimate_D(A,clique.set,nu.hats)
#
# g <- graph_from_adjacency_matrix(A.sub,mode = "undirected")
#
#
# # doesn't seem to work well
#
#
#
#
# library(RColorBrewer)
# c28 <- c(
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "brown", "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "red", "white", "darkgrey"
# )
#
# V(g)$color <- c28[labels]
#
# plot(g,vertex.size= 6,vertex.label=NA)
#
#
#
#
# ######
# # clique.set[[10]] <- NULL
# # clique.set[[11]] <- NULL
# #
# # nu.hats[[10]] <- NULL
# # nu.hats[[11]] <- NULL
#
# K <- length(clique.set)
# clique.idx <- c()
# subset.ids <- c()
# fixed.effect.vec <- c()
# for(k in seq(K)){
#   cat(paste0(k,"  ", length(clique.set[[k]])), end = "\n")
#   clique.idx <- c(clique.idx,clique.set[[k]])
#   subset.ids <- c(subset.ids, rep(k,length(clique.set[[k]])))
#   fixed.effect.vec <- c(fixed.effect.vec,nu.hats[[k]])
# }
#
# for(k in seq(K)){
#   print(length(clique.set[[k]]))
#   print(length(nu.hats[[k]]))
# }
#
# A.cliques <- A[clique.idx,clique.idx]
#
# # for(k1 in seq(K)){
# #   for(k2 in seq(K)){
# #     if(k1 != k2){
# #       idx1 <- as.vector(clique.set[[k1]])
# #       idx2 <- as.vector(clique.set[[k2]])
# #
# #       A.cliques[which(subset.ids == k),which(subset.ids == k)] <- A[idx1,idx2]
# #     }
# #
# #   }
# #
# # }
#
# # based on the huge amount of time it takes to run.
# # I think I will need to write an optimization algorithm
# # from scratch
#
# cutoff.cliques <- 300
#
# subset.ids.small <- subset.ids[subset.ids <= cutoff.cliques]
# fixed.effect.vec.small <- fixed.effect.vec[subset.ids <= cutoff.cliques]
# A.cliques.small <- A.cliques[subset.ids <= cutoff.cliques,subset.ids <= cutoff.cliques]
#
# # library(MASS)
# # write.matrix(A.cliques.small,file="small_A.csv")
# # write.table(fixed.effect.vec.small, "small_nus.csv", row.names = F, col.names = F)
# # write.table(subset.ids.small, "small_cliques.csv", row.names = F, col.names = F)
#
# g.small <- graph_from_adjacency_matrix(A.cliques.small, mode = "undirected")
# V(g.small)$color <- c28[subset.ids.small]
# plot(g.small,vertex.size= 6,vertex.label=NA)
# D0 = init_D0(A.cliques.small,subset.ids.small,fixed.effect.vec.small)
# system.time(D1 <- estimate_D_restricted_fast(A.cliques.small,subset.ids.small,fixed.effect.vec.small,
#                                              D0,thresh = 10**(-3), max.iter = 50, solver = "MOSEK",
#                                              verbose = T))
# print(D1)
# max_triangle_deviation(D1)
# # source("/Users/Owner/mosek/9.3/tools/platform/osx64x86/rmosek/builder.R")
# # attachbuilder()
# # install.rmosek()
# #
# # require(Rmosek)
# # require(CVXR)
# # installed_solvers()
# #sort( sapply(ls(),function(x){object.size(get(x))}))
#
# #system.time(D0 <- estimate_D_restricted(A.cliques.small,subset.ids.small,fixed.effect.vec.small))
#
# #print(D0)
#
# #saveRDS(D.hat, "data/github_cliques_D.rds")
#
# D.hat <- D1
# #saveRDS(D.hat, "data/astrophysics_cliques_D.rds")
#
#
# #D.hat <- readRDS("data/astrophysics_cliques_D.rds")
# mid.search <- optimal_midpoint_search(D.hat,top.k = 10,
#                                       d.yz.min = d.yz.min,
#                                       d.yz.max = d.yz.max)
#
# print(mid.search)
#
#
# y.opt = mid.search[1,1]
# z.opt = mid.search[1,2]
# m.opt = mid.search[1,3]
#
#
# c1 = .5
# c2 = 2
# c3 = .5
#
# x.set <- filter_indices(D.hat, y.opt,
#                         z.opt,m.opt,
#                         c1 = c1,c2 = c2,c3 = c3)
#
#
#
#
# kappa.set.1 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)
#
#
# y.opt = mid.search[2,1]
# z.opt = mid.search[2,2]
# m.opt = mid.search[2,3]
#
# x.set <- filter_indices(D.hat, y.opt,
#                         z.opt,m.opt,
#                         c1 = c1,c2 = c2,c3 = c3)
#
#
# kappa.set.2 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)
#
#
# y.opt = mid.search[3,1]
# z.opt = mid.search[3,2]
# m.opt = mid.search[3,3]
#
# x.set <- filter_indices(D.hat, y.opt,
#                         z.opt,m.opt,
#                         c1 = c1,c2 = c2,c3 = c3)
#
#
# kappa.set.3 <- estimate_kappa_set(D.hat,y.opt,z.opt,m.opt,x.set)
#
#
# curve.scale = 10
#
#
# labels <- c(rep(0,length(kappa.set.1)),rep(1,length(kappa.set.2)),rep(2,length(kappa.set.3)))
# bp.dat <- data.frame(matrix(c(labels,
#                               scale_curvature(kappa.set.1,curve.scale),
#                               scale_curvature(kappa.set.2,curve.scale),
#                               scale_curvature(kappa.set.3,curve.scale)), ncol = 2))
# colnames(bp.dat) <- c("group", "curvature")
#
#
# ggplot(bp.dat, aes(y = curvature, group = group)) +
#   geom_boxplot()
#
#
# med.1 <- median(scale_curvature(kappa.set.1,curve.scale),na.rm = T)
# med.2 <- median(scale_curvature(kappa.set.2,curve.scale),na.rm = T)
# med.3 <- median(scale_curvature(kappa.set.3,curve.scale),na.rm = T)
#
# kappa.vec <- c(kappa.set.1,
#                kappa.set.2,
#                kappa.set.3)
# med.group <- median(scale_curvature(kappa.vec,curve.scale),na.rm = T)
#
# ggplot(bp.dat, aes(y = curvature, x = group)) +
#   geom_jitter() +
#   geom_vline(xintercept = 0.5, col = "red") +
#   geom_vline(xintercept = 1.5, col = "red") +
#   geom_segment(aes(x=-0.5,xend=0.5,y=med.1,yend=med.1), col = "blue") +
#   geom_segment(aes(x=0.5,xend=1.5,y=med.2,yend=med.2), col = "blue") +
#   geom_segment(aes(x=1.5,xend=2.5,y=med.3,yend=med.3), col = "blue") +
#   geom_segment(aes(x=-0.5,xend=2.5,y=med.group,yend=med.group), col = "green")
#
#
# kappa.set.1 <- kappa.set.1[!is.na(kappa.set.1)]
# kappa.set.2 <- kappa.set.2[!is.na(kappa.set.2)]
# kappa.set.3 <- kappa.set.3[!is.na(kappa.set.3)]
#
# index <- c(rep(1,length(kappa.set.1)),
#            rep(2,length(kappa.set.2)),
#            rep(3,length(kappa.set.3)))
# kappa.vec <- c(kappa.set.1,
#                kappa.set.2,
#                kappa.set.3)
# y = scale_curvature(kappa.vec,curve.scale)
# x = index
#
# perm_median_test(y, x)
#
# plot(density(y))
# boxplot(y)
#
# trim.kappa.vec <- scale_curvature(kappa.vec,curve.scale)
#
# # evaluate simulations over time,
# # Is is hard to detect global changes
#
# # begin the work of looking for networks of
# # constant curvature
#


