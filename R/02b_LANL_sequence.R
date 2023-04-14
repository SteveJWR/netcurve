

rm(list = ls())

source("R/00_functions.R")
source("R/clique_finder.R")
set.seed(1)

# download robust changepoint detection algorithm
# https://github.com/guillemr/robust-fpop
library(robseg)
library(ggpubr)

save.plots = T
fig.height = 2000
fig.width = 2000
fig.res = 350

# saved parameters for most costly computational step
reestimate.curvature = F
save.curve.seq = F

device.guide <- read.csv("data/LANL/4h/new_id_dictionary.csv")

#Tuning parameters for dealing with outliers
curve.scale <- 10

day.counter <- 1 #
daily.lag <- 4 # days to average over for connections for cliques
n <- max(device.guide$new_ids)
A.seq <- list()
# though this graph is directed, and weighted, we are only interested in the geometry
A.seq.daily <- list()
for(period in 2:90) {
  # read in the 4h interval cleaned data
  data <- read.csv(paste("data/LANL/4h/netflow_re4h-d",formatC(period,width=2,flag="0"),".csv",sep=""))
  tbl = table(data$TimeFrame)

  A <- sparseMatrix(i = data$SrcDevice,
                    j = data$DstDevice,
                    x = rep(1,nrow(data)),
                    dims = c(n,n))

  A <- A + t(A) # symmetrizing the adjacency matrix
  diag(A) = 0
  A[A > 1] = 1 # triming the double counts

  A.seq.daily[[day.counter]] = A
  cat(paste("Day:",day.counter," ::: Edge_fraction:", round(mean(A),4)), end = "\r")
  day.counter <- day.counter + 1

}


if(reestimate.curvature){
  # adjusting the adjacency matrices so that there is an edge if any signal sent in the last
  # number of days according to the daily lag
  # This is so that we can find cliques on particular days.
  A.seq.daily.lag <- list()
  for(j in seq(length(A.seq.daily))){
    A.new <- A.seq.daily[[j]]
    for(k in seq(0,min(j - 1,daily.lag))){
      A.new <- A.new + A.seq.daily[[j - k]]
    }
    A.new[A.new > 1] = 1
    A.seq.daily.lag[[j]] = A.new
    cat(paste0(j,"/",length(A.seq.daily)), end = "\r")

    ## sparsity is a problem,
    # we assume the random effect are zero

    #A.seq.daily.lag <- A.seq.daily
    # l
    time.seq <- 1:length(A.seq.daily.lag)
    ell <- 5 # 6
    tri.const = sqrt(sqrt(2)) # 1.1
    tri.const = 1.3
    num.cliques.seq <- c()
    kappa.seq <- c()
    kappa.seq2 <- c()
    kappa.med.seq <- c()
    kappa.med.seq2 <- c()
    x.idx <- c()
    x.idx2 <- c()
    missing.clique.times <- c()

    max.iter.estimate = 1

    d.yz.max = log(ell)

    d.yz.min = 1

    set.seed(1)
    # time is delayed due to rolling average
    time.seq <- seq(1,89)
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
      if(length(clique.set) < 6){
        missing.clique.times <- c(missing.clique.times,time)
        next
      }


      D.hat <- lolaR::EstimateD(A, clique.set, max.iter = max.iter.estimate, rand.eff.0 = T)
      #print(paste0("Time: ", time, "/",Time.steps ))
      estimates <- lolaR::EstimateCurvature(D.hat,
                                            J = 1,
                                            tri.const = tri.const,
                                            d.yz.min = d.yz.min,
                                            d.yz.max = d.yz.max)


      kappa.seq2 <- c(kappa.seq2, estimates$kappas)


      kappa.med.seq2 <- c(kappa.med.seq2,median(estimates$kappas,na.rm = T))


      x.idx2 <- c(x.idx2, rep(time,length(estimates$kappas)))

      num.cliques.seq[time] = length(clique.set)
      print(length(clique.set))
      print(paste("Day:",time))

      if(time %% 5 == 0){
        plot(scale_curvature(unlist(kappa.seq2), curve.scale))
        plot(scale_curvature(unlist(kappa.med.seq2), curve.scale))
      }

    }

    if(save.curve.seq){
      saveRDS(kappa.seq2, file = "data/KappaSequence.rds")
      saveRDS(kappa.med.seq2, file = "data/KappaMedianSequence.rds")
      saveRDS(x.idx2, file = "data/DaySequence.rds")
    }

  }
}






kappa.seq2 <- readRDS(file = "data/KappaSequence.rds")
kappa.med.seq2 <- readRDS(file = "data/KappaMedianSequence.rds")
x.idx2 <- readRDS(file = "data/DaySequence.rds")


y <- scale_curvature(kappa.med.seq2, curve.scale)

cpt.true = 57 + daily.lag
lambda.reg = 1.5
est.sd <- mad(diff(y)/sqrt(2))

res.l1 <- Rob_seg.std(x = y,
                      loss = "Outlier",
                      lambda = lambda.reg*est.sd*log(length(y)),
                      lthreshold = 3)
## estimated changepoints
cpt <- res.l1$t.est[-length(res.l1$t.est)]
## simple ploting of changes and smoothed profile
plot(y, pch=20, col="black", ylim = c(-20,5))
lines(res.l1$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=3, col="green")
abline(v=57, lty=3, col="blue")



lambda.reg = .8
est.sd <- mad(diff(y)/sqrt(2))

res.l1.low.smooth <- Rob_seg.std(x = y,
                                 loss = "Outlier",
                                 lambda = lambda.reg*est.sd*log(length(y)),
                                 lthreshold = 3)
## estimated changepoints
cpt <- res.l1.low.smooth$t.est[-length(res.l1.low.smooth$t.est)]
## simple ploting of changes and smoothed profile
plot(y, pch=20, col="black", ylim = c(-20,5))
lines(res.l1.low.smooth$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=3, col="green")
abline(v=57, lty=3, col="blue")


curve.rescaled.seq.dat <- data.frame("idx" = 1:length(y), "curve" = y, "estimate" = res.l1$smt, "estimate2" = res.l1.low.smooth$smt)


plt.curve.cpt <- ggplot(curve.rescaled.seq.dat, aes(x = idx, y = curve)) +
  geom_point() +
  geom_line(aes(x = idx, y = estimate), color = "red") +
  #geom_vline(xintercept = cpt.true, color = "blue") +
  geom_vline(xintercept = 57, color = "purple") +
  ylab("Transformed Curvature") +
  xlab("Day") +
  theme_bw() +
  ggtitle("LANL Curvature Changepoint Detection")

plt.curve.cpt2 <- ggplot(curve.rescaled.seq.dat, aes(x = idx, y = curve)) +
  geom_point() +
  geom_line(aes(x = idx, y = estimate2), color = "red") +
  #geom_vline(xintercept = cpt.true, color = "blue") +
  geom_vline(xintercept = 57, color = "purple") +
  ylab("Transformed Curvature") +
  xlab("Day") +
  theme_bw() +
  ggtitle("LANL Curvature Changepoint Detection (low lambda)")



#dev.off()
if(save.plots){
  plt.curve.cpt %>%
    ggexport(filename = "plots/LANL_curve_changepoint_high_reg.png", width = fig.width, height = fig.height, res = fig.res)
  plt.curve.cpt2 %>%
    ggexport(filename = "plots/LANL_curve_changepoint_low_reg.png", width = fig.width, height = fig.height, res = fig.res)

}




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


lambda.red.tri = 10**(-7)
est.sd <- mad(diff(tri.seq)/sqrt(2))
res.l1.tri <- Rob_seg.std(x = tri.seq,
                          loss = "Outlier",
                          lambda = lambda.red.tri*est.sd*log(length(tri.seq)),
                          lthreshold = 3)
## estimated changepoints
res.l1.tri$smt
cpt <- res.l1.tri$t.est[-length(res.l1.tri$t.est)]
## simple ploting of changes and smoothed profile
plot(tri.seq, pch=20, col="black")
lines(res.l1.tri$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="purple",lwd=2)


lambda.red.edge = 5*10**(-4)
est.sd <- mad(diff(edge.dens.seq)/sqrt(2))
res.l1.edge <- Rob_seg.std(x = edge.dens.seq,
                           loss = "Outlier",
                           lambda = lambda.red.edge*est.sd*log(length(edge.dens.seq)),
                           lthreshold = 3)
## estimated changepoints
res.l1.edge$smt
cpt <- res.l1.edge$t.est[-length(res.l1.edge$t.est)]
## simple ploting of changes and smoothed profile
plot(edge.dens.seq, pch=20, col="black")
lines(res.l1.edge$smt, col="red", lwd=2)
abline(v=cpt, lty=2, col="red")
abline(v=cpt.true, lty=4, col="purple",lwd=2)



graph.stats.dat <- data.frame("day" = 1:length(edge.dens.seq),
                              "edge.dens" = edge.dens.seq,
                              "tri.dens" = tri.seq,
                              "pred.edge" = res.l1.edge$smt,
                              "pred.tri" = res.l1.tri$smt)


plt.edge <- ggplot(graph.stats.dat, aes(x = day, y = edge.dens)) +
  geom_point() +
  geom_line(aes(x = day, y = pred.edge), color = "red") +
  geom_vline(xintercept = cpt.true, color = "purple") +
  ylab("Edge Density") +
  xlab("Day") +
  theme_bw() +
  ggtitle("LANL Edge Density Changepoint Detection")



if(save.plots){
  plt.edge %>%
    ggexport(filename = "plots/LANL_edge_changepoint.png", width = fig.width, height = fig.height, res = fig.res)

}



plt.tri <- ggplot(graph.stats.dat, aes(x = day, y = tri.dens)) +
  geom_point() +
  geom_line(aes(x = day, y = pred.tri), color = "red") +
  geom_vline(xintercept = cpt.true, color = "purple") +
  ylab("Triangle Density") +
  xlab("Day") +
  theme_bw() +
  ggtitle("LANL Triangle Density Changepoint Detection")




if(save.plots){
  plt.tri %>%
    ggexport(filename = "plots/LANL_triangle_changepoint.png", width = fig.width, height = fig.height, res = fig.res)

}






# Plots minimizing time after event as a function of number of checks.

beta.seq <- exp(seq(log(0.00001),log(13),length.out = 3000))
#beta.seq <- seq(0.001,13,length.out = 1000)
#plot(beta.seq)
curve.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
tri.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
edge.dat <- matrix(NA,nrow = length(beta.seq), ncol = 2)
event.day <- 57

for(beta.idx in seq(beta.seq)){
  beta = beta.seq[beta.idx]

  est.sd <- mad(diff(y)/sqrt(2))
  res.l1 <- Rob_seg.std(x = y,
                        loss = "Outlier",
                        lambda = beta*est.sd*log(length(y)),
                        lthreshold = 3)
  ## estimated changepoints
  cpt.days <- res.l1$t.est[-length(res.l1$t.est)]


  #minimum lag until day after detection
  day.lag <- min(cpt.days[cpt.days >= event.day]) - event.day
  # for plotting purposes
  if(day.lag > 100){
    day.lag = 100
  }
  detect.rate <- length(unique(cpt.days))

  curve.dat[beta.idx,] <- c(day.lag,detect.rate)


  est.sd.tri <- mad(diff(tri.seq)/sqrt(2))
  res.l1.tri <- Rob_seg.std(x = tri.seq,
                            loss = "Outlier",
                            lambda = 10**(-5)*beta*est.sd.tri*log(length(tri.seq)),
                            lthreshold = 3)
  res.l1.tri$smt
  cpt.days <- res.l1.tri$t.est[-length(res.l1.tri$t.est)]


  day.lag <- min(cpt.days[cpt.days >= event.day]) - event.day
  # for plotting purposes
  if(day.lag > 100){
    day.lag = 100
  }
  detect.rate <- length(cpt.days)

  tri.dat[beta.idx,] <- c(day.lag,detect.rate)







  est.sd.edge <- mad(diff(edge.dens.seq)/sqrt(2))
  res.l1.edge <- Rob_seg.std(x = edge.dens.seq,
                             loss = "Outlier",
                             lambda = 10**(-4)*beta*est.sd.edge*log(length(edge.dens.seq)),
                             lthreshold = 3)
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

alarm.plot <- ggplot(data = plot.dat, aes(x = Budget, y = TimeDelay, col = Method)) +
                      geom_line() +
                      ggtitle("Alarm Rate And Detection Delay of Red Team Attack") +
                      xlab("Alarm Rate") +
                      theme_bw() +
                      coord_cartesian(ylim=c(0,25))


if(save.plots){
  alarm.plot %>%
    ggexport(filename = "plots/LANL_alarm_rate.png", width = fig.width, height = fig.height, res = fig.res)
}

#### END  ####
