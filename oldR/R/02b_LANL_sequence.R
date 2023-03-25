

rm(list = ls())

source("00_functions.R")
source("clique_finder.R")
set.seed(1)

# https://github.com/guillemr/robust-fpop
library(robseg)


png.width = 1500  
png.height = 1500
png.res = 300


device.guide <- read.csv("data/LANL/4h/new_id_dictionary.csv")

curve.scale <- 10 
counter <- 1
time.counter <- 1
day.counter <- 1
daily.lag <- 4 
n <- max(device.guide$new_ids)
A.seq <- list()
# though this graph is directed, and weighted, we are only interested in the
# piece of geometry
A.seq.daily <- list()
for(period in 2:90) {
  data <- read.csv(paste("data/LANL/4h/netflow_re4h-d",formatC(period,width=2,flag="0"),".csv",sep=""))
  tbl = table(data$TimeFrame)
  # for(i in sort(unique(data$TimeFrame))){
  #   idx <- which(data$TimeFrame == i)
  #   data.sub <- data[idx,]
  #   
  #   A.sub <- sparseMatrix(i = data.sub$SrcDevice,
  #                         j = data.sub$DstDevice, 
  #                         x = rep(1,nrow(data.sub)),
  #                         dims = c(n,n))
  #   
  #   A.sub <- A.sub + t(A.sub)
  #   diag(A.sub) = 0
  #   A.sub[A.sub > 1] = 1
  #   
  #   A.seq[[time.counter]] = A.sub
  #   cat(paste("Period:",period," ::: Edge_fraction:", round(mean(A.sub),4)), end = "\r")
  #   time.counter <- time.counter + 1
  #   
  # }
  A <- sparseMatrix(i = data$SrcDevice,
                    j = data$DstDevice, 
                    x = rep(1,nrow(data)),
                    dims = c(n,n))
  
  A <- A + t(A)
  diag(A) = 0
  A[A > 1] = 1
  
  A.seq.daily[[day.counter]] = A
  cat(paste("Day:",day.counter," ::: Edge_fraction:", round(mean(A),4)), end = "\r")
  day.counter <- day.counter + 1
  #counter <- counter + i
}



A.seq.daily.lag <- list()
for(j in seq(length(A.seq.daily))){
  A.new <- A.seq.daily[[j]]
  for(k in seq(0,min(j - 1,2))){
    A.new <- A.new + A.seq.daily[[j - k]]
  }
  A.new[A.new > 1] = 1
  A.seq.daily.lag[[j]] = A.new
}



## sparsity is a problem, we deviate from the model class slightly
# 

#A.seq.daily.lag <- A.seq.daily
time.seq <- 1:length(A.seq.daily.lag)
ell <- 5 # 6
tri.const = sqrt(sqrt(2)) # 1.1
tri.const = 1.5
num.cliques.seq <- c()
kappa.seq <- c()
kappa.seq2 <- c()
kappa.med.seq <- c()
kappa.med.seq2 <- c()
x.idx <- c()
x.idx2 <- c()
missing.clique.times <- c()
set.seed(1)
# time is delayed three days 
for(time in time.seq){
  if(time == 57){
    original.rta <- length(kappa.seq) + 1
    original.rta2 <- length(kappa.seq2) + 1
  }
  if(time == 57 + daily.lag){
    approximate.rta <- length(kappa.seq) + 1
    approximate.rta2 <- length(kappa.seq2) + 1
  }
  A <- A.seq.daily.lag[[time]]
  g <- igraph::graph_from_adjacency_matrix(A,mode = "undirected")
  clique.set <- igraph::maximal.cliques(g, min = ell)
  clique.set <- clique_partition(clique.set,randomize = T)
  clique.set <- clique_set_connected_subgraph(A,clique.set) #reduces to largest connected component
  #plot_cliques(A,clique.set)
  if(length(clique.set) < 4){
    missing.clique.times <- c(missing.clique.times,time)
    next
  }
  #clique.set <- clique_finder_2(A, min_clique_size = ell)
  # misspecification of the model? 
  # large cliques in an otherwise sparse graph 
  # would suggest that this is a small difference
  # 
  
  # estimates = estimate_curvature(A, 
  #                                clique.set, 
  #                                verbose = F, 
  #                                tri.const = tri.const,
  #                                d.yz.min = .2, 
  #                                d.yz.max = 4,
  #                                rand.eff.0 = F)
  
  estimates2 = estimate_curvature(A, 
                                  clique.set, 
                                  verbose = F, 
                                  d.yz.min = 1,
                                  d.yz.max = 4,
                                  tri.const = tri.const,
                                  rand.eff.0 = T)
  

  
  
  #kappa.seq <- c(kappa.seq, estimates$kappas)
  kappa.seq2 <- c(kappa.seq2, estimates2$kappas)
  
  #kappa.med.seq <- c(kappa.med.seq,median(estimates$kappas,na.rm = T))
  kappa.med.seq2 <- c(kappa.med.seq2,median(estimates2$kappas,na.rm = T))
  
  #x.idx <- c(x.idx, rep(time,length(estimates$kappas)))
  x.idx2 <- c(x.idx2, rep(time,length(estimates2$kappas)))
  
  num.cliques.seq[time] = length(clique.set)
  print(length(clique.set))
  print(paste("Day:",time))
  print("")
  
  if(time %% 5 == 0){
    
    #plot(scale_curvature(unlist(kappa.seq), curve.scale))
    plot(scale_curvature(unlist(kappa.seq2), curve.scale))
  }

}

#y = unlist(kappa.seq)
y2 = unlist(kappa.seq2)
length(y2)


# 
# changepoint.data <- list("kappas" = y, "time.idx" = x.idx, 
#                          "rta.idx" = approximate.rta)
# changepoint.data2 <- list("kappas" = y2, "time.idx" = x.idx2,
#                           "rta.idx" = approximate.rta2, 
#                           "rta.idx.init" = original.rta2)
# file <- "results/LANL_netflow_estimates_ell_5_tri_1p6.rds"
# file.no.randeff <- "results/LANL_netflow_estimates_no_randeff_ell_5_tri_1p6.rds"
# saveRDS(changepoint.data,file)
# saveRDS(changepoint.data2,file.no.randeff)




###### scale sequence

file.no.randeff <- "results/LANL_netflow_estimates_no_randeff_ell_5_tri_1p6.rds"

changepoint.data2 <- readRDS(file.no.randeff)
y2 <- changepoint.data2$kappas
x.idx2 <- changepoint.data2$time.idx
rta.idx2 <- changepoint.data2$rta.idx
rta.idx.init2 <- changepoint.data2$rta.idx.init
approximate.rta2 <- rta.idx2
original.rta2 <- rta.idx.init2


y.clean2 <- scale_curvature(y2, curve.scale)
non.na.idx2 <- which(!is.na(y2))
y.clean2 <- y.clean2[non.na.idx2]
x.idx2.clean <- x.idx2[non.na.idx2]
cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)


### rescaling within every time block
y2.clean.rescaled <- c()
for(k in seq(max(x.idx2.clean))){
  
  y2.clean.block <- y.clean2[x.idx2.clean == k]
  y2.clean.rescale.block <- (y2.clean.block - med(y2.clean.block))/mean(abs(y2.clean.block - med(y2.clean.block))) + med(y2.clean.block)
  #y2.clean.rescale.block <- (y2.clean.block - med(y2.clean.block))/mad(y2.clean.block) + med(y2.clean.block)
  if(mean(abs(y2.clean.block - med(y2.clean.block))) == 0){
    y2.clean.rescale.block <- y2.clean.block
  }
  if(any(is.na(y2.clean.rescale.block))){
    break 
  }
  if(any(is.infinite(y2.clean.rescale.block))){
    break 
  }
  y2.clean.rescaled <- c(y2.clean.rescaled,y2.clean.rescale.block)
}




#approximate.rta <- 513
#approximate.rta2 <- 375
#cpt.true <- approximate.rta 
#plot(scale_curvature(y, curve.scale))
plot(scale_curvature(y2, curve.scale))


y.clean <- scale_curvature(y, curve.scale)
non.na.idx <- which(!is.na(y))
cpt.true <- which.min(non.na.idx <= approximate.rta)

# non.na.idx2 <- which(!is.na(y2))
# cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)

y.clean <- y.clean[non.na.idx]
est.sd <- mad(diff(y.clean)/sqrt(2))
est.sd <- 1
res.l1 <- Rob_seg.std(x = y.clean,  
                      loss = "L1", 
                      lambda = (2)*est.sd*log(length(y.clean)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y.clean, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=3, col="green")


# y.clean2 <- scale_curvature(y2, curve.scale)
# 
# non.na.idx2 <- which(!is.na(y2))
# cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)
# 
# y.clean2 <- y.clean2[non.na.idx2]
# est.sd <- mad(diff(y.clean2)/sqrt(2))
# 
# res.l1 <- Rob_seg.std(x = y.clean2,  
#                       loss = "L2", 
#                       lambda = (3)*est.sd*log(length(y.clean2)))
# ## estimated changepoints 
# cpt <- res.l1$t.est[-length(res.l1$t.est)]
# ## simple ploting of changes and smoothed profile
# plot(y.clean2, pch=20, col="black")
# lines(res.l1$smt, col="red", lwd=2)
# abline(v=cpt, lty=2, col="red")
# abline(v=cpt.true2, lty=3, col="green")



y.clean2 <- scale_curvature(y2, curve.scale)
non.na.idx2 <- which(!is.na(y2))
cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)
cpt.original.true2 <- which.min(non.na.idx2 <= original.rta2)

y.clean2 <- y.clean2[non.na.idx2]
est.sd <- mad(diff(y.clean2)/sqrt(2))
res.l1 <- Rob_seg.std(x = y.clean2,  
                      loss = "L1", 
                      lambda = (6)*est.sd*log(length(y.clean2)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y.clean2, pch=20, col="black")

lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true2, lty=4, col="green",lwd=2)
abline(v=cpt.original.true2, lty=4, col="blue",lwd=2)


non.na.idx2 <- which(!is.na(y2.clean.rescaled))
y2.clean.rescaled <- y2.clean.rescaled[non.na.idx2]
cpt.true2 <- approximate.rta2
cpt.original.true2 <- original.rta2
est.sd <- mad(diff(y2.clean.rescaled)/sqrt(2))
res.l1 <- Rob_seg.std(x = y2.clean.rescaled,  
                      loss = "L1", 
                      lambda = 11*est.sd*log(length(y2.clean.rescaled)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y2.clean.rescaled, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true2, lty=4, col="green",lwd=2)
abline(v=cpt.original.true2, lty=4, col="blue",lwd=2)



res.l1.low.smooth <- Rob_seg.std(x = y2.clean.rescaled,  
                                 loss = "L1", 
                                 lambda = 10*est.sd*log(length(y2.clean.rescaled)))

curve.rescaled.seq.dat <- data.frame("idx" = 1:length(y2.clean.rescaled), "curve" = y2.clean.rescaled, "estimate" = res.l1$smt, "estimate2" = res.l1.low.smooth$smt)


plt.curve.cpt <- ggplot(curve.rescaled.seq.dat, aes(x = idx, y = curve)) + 
  geom_point() + 
  geom_line(aes(x = idx, y = estimate), color = "red") + 
  geom_vline(xintercept = cpt.true2, color = "green") + 
  ylab("Curvature") + 
  xlab("idx (Ordered by Day)") + 
  ggtitle("LANL Curvature Changepoint Detection")

plt.curve.cpt2 <- ggplot(curve.rescaled.seq.dat, aes(x = idx, y = curve)) + 
  geom_point() + 
  geom_line(aes(x = idx, y = estimate2), color = "red") + 
  geom_vline(xintercept = cpt.true2, color = "green") + 
  ylab("Curvature") + 
  xlab("idx (ordered by day)") + 
  ggtitle("LANL Curvature Changepoint Detection (low lambda)")





png(filename = "plots/LANL_curve_changepoint_high_reg.png",
    width = 4*png.width, height = 4*png.height, res = 4*png.res)


plt.curve.cpt  

# Close the pdf file
dev.off() 


png(filename = "plots/LANL_curve_changepoint_low_reg.png",
    width = 4*png.width, height = 4*png.height, res = 4*png.res)


plt.curve.cpt2  

# Close the pdf file
dev.off() 



## comparison against edge and triangle densities 
cpt.true <- 57 + daily.lag

tri.seq <- c()
edge.dens.seq <- c()

# replace back daily lag 
for(time in seq(length(A.seq.daily))){
  print(time)
  A <- A.seq.daily[[time]]
  n = nrow(A)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  n.tri <- sum(igraph::count_triangles(g))/3
  n.edge <- sum(A)
  tri.seq[time] <- n.tri/choose(n,3)
  edge.dens.seq[time] <- n.edge/choose(n,2)
}


est.sd <- mad(diff(tri.seq)/sqrt(2))
res.l1.tri <- Rob_seg.std(x = tri.seq,  
                          loss = "L1", 
                          lambda = (1.5)*est.sd*log(length(tri.seq)))
## estimated changepoints 
res.l1.tri$smt
cpt <- res.l1.tri$t.est[-length(res.l1.tri$t.est)]
## simple ploting of changes and smoothed profile
plot(tri.seq, pch=20, col="black")
lines(res.l1.tri$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="green",lwd=2)


est.sd <- mad(diff(edge.dens.seq)/sqrt(2))
res.l1.edge <- Rob_seg.std(x = edge.dens.seq,  
                           loss = "L1", 
                           lambda = (0.9)*est.sd*log(length(edge.dens.seq)))
## estimated changepoints 
res.l1.edge$smt
cpt <- res.l1.edge$t.est[-length(res.l1.edge$t.est)]
## simple ploting of changes and smoothed profile
plot(edge.dens.seq, pch=20, col="black")
lines(res.l1.edge$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="green",lwd=2)



graph.stats.dat <- data.frame("day" = 1:length(edge.dens.seq), 
                              "edge.dens" = edge.dens.seq, 
                              "tri.dens" = tri.seq,
                              "pred.edge" = res.l1.edge$smt, 
                              "pred.tri" = res.l1.tri$smt)


plt.edge <- ggplot(graph.stats.dat, aes(x = day, y = edge.dens)) + 
  geom_point() + 
  geom_line(aes(x = day, y = pred.edge), color = "red") + 
  geom_vline(xintercept = cpt.true, color = "green") + 
  ylab("Edge Density") + 
  xlab("Day") + 
  ggtitle("LANL Edge Density Changepoint Detection")


png(filename = "plots/LANL_edge_changepoint.png",
    width = 4*png.width, height = 4*png.height, res = 4*png.res)


plt.edge  

# Close the pdf file
dev.off() 


plt.tri <- ggplot(graph.stats.dat, aes(x = day, y = tri.dens)) + 
  geom_point() + 
  geom_line(aes(x = day, y = pred.tri), color = "red") + 
  geom_vline(xintercept = cpt.true, color = "green") + 
  ylab("Triangle Density") + 
  xlab("Day") + 
  ggtitle("LANL Triangle Density Changepoint Detection")


png(filename = "plots/LANL_triangle_changepoint.png",
    width = 4*png.width, height = 4*png.height, res = 4*png.res)


plt.tri  

# Close the pdf file
dev.off() 





# Plots minimizing time after event as a function of number of checks. 

beta.seq <- exp(seq(log(0.001),log(13),length.out = 3000))
#beta.seq <- seq(0.001,13,length.out = 1000)
#plot(beta.seq)
curve.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
tri.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
edge.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
event.day <- 57

for(beta.idx in seq(beta.seq)){
  beta = beta.seq[beta.idx]
  
  est.sd <- mad(diff(y2.clean.rescaled)/sqrt(2))
  res.l1 <- Rob_seg.std(x = y2.clean.rescaled,  
                        loss = "L1", 
                        lambda = beta*est.sd*log(length(y2.clean.rescaled)))
  ## estimated changepoints 
  cpt <- res.l1$t.est[-length(res.l1$t.est)]
  cpt.days <- x.idx2.clean[cpt]
  
  
  day.lag <- min(cpt.days[cpt.days >= event.day]) - event.day
  # for plotting purposes
  if(day.lag > 100){
    day.lag = 100
  }
  detect.rate <- length(unique(cpt.days))
  
  curve.dat[beta.idx,] <- c(day.lag,detect.rate)
  
  
  est.sd <- mad(diff(tri.seq)/sqrt(2))
  res.l1.tri <- Rob_seg.std(x = tri.seq,  
                            loss = "L1", 
                            lambda = beta*est.sd*log(length(tri.seq)))
  res.l1.tri$smt
  cpt.days <- res.l1.tri$t.est[-length(res.l1.tri$t.est)]
  
  
  day.lag <- min(cpt.days[cpt.days >= event.day]) - event.day
  # for plotting purposes
  if(day.lag > 100){
    day.lag = 100
  }
  detect.rate <- length(cpt.days)
  
  tri.dat[beta.idx,] <- c(day.lag,detect.rate)
  
  
  
  
  
  
  
  est.sd <- mad(diff(edge.dens.seq)/sqrt(2))
  res.l1.edge <- Rob_seg.std(x = edge.dens.seq,  
                             loss = "L1", 
                             lambda = beta*est.sd*log(length(edge.dens.seq)))
  res.l1.edge$smt
  cpt.days <- res.l1.edge$t.est[-length(res.l1.edge$t.est)]
  
  day.lag <- min(cpt.days[cpt.days >= event.day]) - event.day
  # for plotting purposes
  if(day.lag > 100){
    day.lag = 100
  }
  detect.rate <- length(cpt.days)
  
  edge.dat[beta.idx,] <- c(day.lag,detect.rate)
  
  
}


plot(tri.dat[,2],tri.dat[,1], col = "blue", ylim = c(0,100), type = "l")
lines(edge.dat[,2],edge.dat[,1], col = "red")
lines(curve.dat[,2],curve.dat[,1], col = "green")


plot.dat <- data.frame(matrix(NA,nrow = 0, ncol = 3))

plot.dat <- rbind(plot.dat,cbind(rep("curvature", nrow(curve.dat)), curve.dat))
plot.dat <- rbind(plot.dat,cbind(rep("edges", nrow(edge.dat)), edge.dat))
plot.dat <- rbind(plot.dat,cbind(rep("triangles", nrow(tri.dat)), tri.dat))


colnames(plot.dat) = c("Method", "TimeDelay","Budget")
plot.dat$Budget <- as.numeric(plot.dat$Budget)
plot.dat$TimeDelay <- as.numeric(plot.dat$TimeDelay)

ggplot(data = plot.dat, aes(x = Budget, y = TimeDelay, col = Method)) + 
  geom_line() + 
  ggtitle("Alarm Rate And Detection Delay of Red Team Attack") + 
  xlab("Alarm Rate") + 
  ylim(c(0,25))


#### END  ####
# 
# est.sd <- mad(diff(y2.clean.rescaled)/sqrt(2))
# res.l1 <- Rob_seg.std(x = y2.clean.rescaled,  
#                       loss = "L1", 
#                       lambda = (1.5)*est.sd*log(length(y.clean2)))
# ## estimated changepoints 
# cpt <- res.l1$t.est[-length(res.l1$t.est)]
# ## simple ploting of changes and smoothed profile
# plot(y2.clean.rescaled, pch=20, col="black")
# lines(res.l1$smt, col="red", lwd=2)
# abline(v=cpt, lty=2, col="red")
# abline(v=cpt.true2, lty=4, col="green",lwd=2)
# 
# 
# y.clean2 <- scale_curvature(y2, curve.scale)
# non.na.idx2 <- which(!is.na(y2))
# cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)
# 
# y.clean2 <- y.clean2[non.na.idx2]
# est.sd <- mad(diff(y.clean2)/sqrt(2))
# res.l1 <- Rob_seg.std(x = y.clean2,  
#                       loss = "Outlier", 
#                       lambda = (.1)*est.sd*log(length(y.clean2)),
#                       lthreshold=.7)
# ## estimated changepoints 
# cpt <- res.l1$t.est[-length(res.l1$t.est)]
# ## simple ploting of changes and smoothed profile
# plot(y.clean2, pch=20, col="black")
# lines(res.l1$smt, col="red", lwd=2)
# abline(v=cpt, lty=2, col="red")
# abline(v=cpt.true2, lty=3, col="green")
# 
# 
# y.clean2 <- scale_curvature(y2, curve.scale)
# non.na.idx2 <- which(!is.na(y2))
# cpt.true2 <- which.min(non.na.idx2 <= approximate.rta2)
# 
# y.clean2 <- y.clean2[non.na.idx2]
# est.sd <- mad(diff(y.clean2)/sqrt(2))
# res.l1 <- Rob_seg.std(x = y.clean2,  
#                       loss = "Huber", 
#                       lambda = (4.5)*est.sd*log(length(y.clean2)),
#                       lthreshold=sqrt(2))
# ## estimated changepoints 
# cpt <- res.l1$t.est[-length(res.l1$t.est)]
# ## simple ploting of changes and smoothed profile
# plot(y.clean2, pch=20, col="black")
# lines(res.l1$smt, col="red", lwd=2)
# abline(v=cpt, lty=2, col="red")
# abline(v=cpt.true2, lty=3, col="green")





y.med.clean <- scale_curvature(unlist(kappa.med.seq), curve.scale)
non.na.idx <- which(!is.na(y.med.clean))
cpt.true <- 57

y.clean <- y.clean[non.na.idx]
est.sd <- mad(diff(y.med.clean)/sqrt(2))
res.l1 <- Rob_seg.std(x = y.med.clean,  
                      loss = "L1", 
                      lambda = (.4)*est.sd*log(length(y.clean)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y.med.clean, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=3, col="green")

y.med.clean2 <- scale_curvature(unlist(kappa.med.seq2), curve.scale)
non.na.idx <- which(!is.na(y.med.clean2))
cpt.true2 <- 57

est.sd <- mad(diff(y.med.clean2)/sqrt(2))
res.l1 <- Rob_seg.std(x = y.med.clean,  
                      loss = "L1", 
                      lambda = (.7)*est.sd*log(length(y.clean)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y.med.clean2, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true2, lty=3, col="green")


med.seq2 <- c()

for(i in 1:89){
  idx <- which(x.idx2 == i)
  med.seq2[i] <- median(scale_curvature(unlist(kappa.seq2[idx]),curve.scale),na.rm = T)
}


no.na.idx <- !is.na(med.seq2)
med.idx2 <- 1:89
med.original.cpt2 <- which.max(med.idx2[no.na.idx] >= 57)
med.cpt2 <- which.max(med.idx2[no.na.idx] >= 57 + daily.lag)
med.seq2 <- med.seq2[no.na.idx]




est.sd <- mad(diff(med.seq2)/sqrt(2))
res.l1 <- Rob_seg.std(x = med.seq2,  
                      loss = "L1", 
                      lambda = (.8)*est.sd*log(length(med.seq2)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(med.seq2, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=med.original.cpt2, lty=3, col="blue")
abline(v=med.cpt2, lty=3, col="green")





plot(tri.seq)
plot(edge.dens.seq)






est.sd <- mad(diff(edge.dens.seq)/sqrt(2))
res.l1 <- Rob_seg.std(x = edge.dens.seq,  
                      loss = "L1", 
                      lambda = (.8)*est.sd*log(length(edge.dens.seq)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(edge.dens.seq, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="green",lwd=2)


est.sd <- mad(diff(tri.seq)/sqrt(2))
res.l1 <- Rob_seg.std(x = tri.seq,  
                      loss = "L1", 
                      lambda = (1.0)*est.sd*log(length(tri.seq)))
## estimated changepoints 
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(tri.seq, pch=20, col="black")
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="green",lwd=2)












