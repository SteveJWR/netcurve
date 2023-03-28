


rm(list = ls())

source("00_functions.R")
source("clique_finder.R")
set.seed(1)
library(latex2exp)
library(scales)
### Convergence Results 

png.width = 1500
png.height = 1500
png.res = 300

kappa.set <- c(-2,-1,-0.5,0,0.5,1)
#estimates 
ell.set <- c(6,8,12,16)

plot.dat <- matrix(NA, 0,3)
tri.const.filter <- 1.4

for(kappa in kappa.set){
  
  full.estimates <- matrix(NA, nrow = 0, ncol = 4)
  for(block in 1:20){
    if(kappa == 0.5 & block == 9){
      next
    }
    file <-  paste0("results/estimates_kappa_",kappa,"_block_",block,".csv")
    if(file.exists(file)){
      est.block <- read.csv(file)
      rownames(est.block) = est.block[,1]
      est.block <- est.block[,2:6]
      est.block <- est.block[est.block[,5] == tri.const.filter, ]
      full.estimates <- rbind(full.estimates, est.block[,1:4])
    }
    
    ell.seq <- rep(ell.set, each = nrow(full.estimates))
    
    dat.tmp = matrix(c(ell.seq, as.numeric(unlist(full.estimates))), ncol = 2, nrow = 4*nrow(full.estimates))
    dat.tmp <- cbind(dat.tmp, rep(kappa, nrow(full.estimates)))
    dat.tmp[,2] <- dat.tmp[,2] - kappa
  }
  plot.dat <- rbind(plot.dat, dat.tmp)
}

plot.dat <- as.data.frame(plot.dat)




colnames(plot.dat) = c("CliqueSize","Bias", "Curvature")
plot.dat$Curvature = as.factor(plot.dat$Curvature)
plot.summary.data <- plot.dat %>% group_by(CliqueSize, Curvature) %>% 
  summarize(q9 = quantile(Bias, 0.9, na.rm = T), 
            q1  = quantile(Bias, 0.1, na.rm = T), 
            med = median(Bias, na.rm = T))


plt <- ggplot(data = plot.summary.data, aes(x = CliqueSize, y = med, color = Curvature)) + 
  geom_line() + 
  geom_errorbar(aes(ymax = q9, ymin = q1),position=position_dodge(width=0.5)) + 
  geom_point(position=position_dodge(width=0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  #ggtitle(paste0(TeX("Consistency of Estimator: $C_{\\Delta}=$"), tri.const)) + 
  ggtitle(TeX(sprintf(r'(Consistency of Estimator: $C_{\Delta}=%f$)', tri.const.filter))) +
  ylab("Estimate Deviation")+ 
  xlab("Clique Size") + 
  scale_x_continuous(breaks = ell.set,labels = ell.set)

 
plt


# png(filename = "plots/GMM_consistancy.png",
#     width = png.width, height = png.height, res = png.res)
# 
# 
# plt  
# 
# # Close the pdf file
# dev.off() 




### Coverage Results 

kappa.set <- c(-2,-1,-0.5,0,0.5,1)
#estimates 
ell.set <- c(6,8,12,16)

p.val.plot.dat <- matrix(NA, 0,3)
tri.const.filter <- 1.0
for(kappa in kappa.set){
  
  p.val.full <- matrix(NA, nrow = 0, ncol = 4)
  for(block in 1:10){
    
    #file <-  paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
    file <-  paste0("results/norm_p_vals_kappa_",kappa,"_block_",block,".csv")
    if(file.exists(file)){
      p.val.block <- read.csv(file)
      rownames(p.val.block) = p.val.block[,1]
      p.val.block <- p.val.block[,2:6]
      p.val.block <- p.val.block[p.val.block[,5] == tri.const.filter, ]
      p.val.full <- rbind(p.val.full, p.val.block[,1:4])
      ell.seq <- rep(ell.set, each = nrow(p.val.full))
      dat.tmp = matrix(c(ell.seq, as.numeric(unlist(p.val.full))), ncol = 2, nrow = 4*nrow(p.val.full))
      dat.tmp <- cbind(dat.tmp, rep(kappa, nrow(p.val.full)))
    }
    
  }
  p.val.plot.dat <- rbind(p.val.plot.dat, dat.tmp)
}

p.val.plot.dat <- as.data.frame(p.val.plot.dat)




colnames(p.val.plot.dat) = c("CliqueSize","pval", "Curvature")
p.val.plot.dat$Curvature = as.factor(p.val.plot.dat$Curvature)
p.val.plot.summary.data <- p.val.plot.dat %>% group_by(CliqueSize, Curvature) %>% 
  summarize(fpr1 = mean(pval < 0.1, na.rm = T),
            fpr05 = mean(pval < 0.05, na.rm = T),
            fpr01  = mean(pval < 0.01, na.rm = T),
            fpr005  = mean(pval < 0.005, na.rm = T),
            var.fpr1 = (1/(sum(pval > 0.0, na.rm = T)))*mean(pval < 0.1, na.rm = T)*(1 - mean(pval < 0.1, na.rm = T)),
            var.fpr05 = (1/(sum(pval > 0.0, na.rm = T)))*mean(pval < 0.05, na.rm = T)*(1 - mean(pval < 0.05, na.rm = T)),
            var.fpr01  = (1/(sum(pval > 0.0, na.rm = T)))*mean(pval < 0.01, na.rm = T)*(1 - mean(pval < 0.01, na.rm = T)),
            var.fpr005  = (1/(sum(pval > 0.0, na.rm = T)))*mean(pval < 0.005, na.rm = T)*(1 - mean(pval < 0.005, na.rm = T)))    


ggplot(data = p.val.plot.dat, 
       aes(x = CliqueSize, y = pval, color = Curvature)) + 
  geom_point() + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  ggtitle("Coverage Of Test") + 
  ylab("Estimate Deviation")


plt <- ggplot(data = p.val.plot.summary.data, 
       aes(x = CliqueSize, y = fpr05, color = Curvature)) + 
  geom_line() + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  # geom_errorbar(aes(ymax = fpr05 + 2*sqrt(var.fpr05), 
  #                   ymin = fpr05 - 2*sqrt(var.fpr05))) +  
  ggtitle(TeX(sprintf(r'(FPR of Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  ylab("False Positive Rate") +
  xlab("Clique Size") + 
  ylim(c(0,1)) + 
  scale_x_continuous(breaks = ell.set,labels = ell.set)
plt

png(filename = "plots/fpr_C1p0.png",
    width = png.width, height = png.height, res = png.res)
plt
# Close the pdf file
dev.off()

tmp <- p.val.plot.dat[p.val.plot.dat$CliqueSize == 12 & p.val.plot.dat$Curvature == -1,  2]
plot(hist(tmp))

x = seq(0,1,length.out = length(tmp))
qqplot(x,tmp)

# png(filename = "plots/GMM_coverage.png",
#     width = png.width, height = png.height, res = png.res)
# 
# 
# plt  
# 
# # Close the pdf file
# dev.off() 




#### Adjacent Spheres Power #####
# we could update this later, add more simulations? 

ell.set <- c(6,8,12,16)

ad.sphere.p.val.full <- matrix(NA, nrow = 0, ncol = 4)

tri.const.filter = 1.6
for(block in 1:100){
  
  file <-  paste0("results/adjacent_spheres_results_block_",block,".csv")
  if(file.exists(file)){
    p.val.block <- read.csv(file)
    rownames(p.val.block) = p.val.block[,1]
    p.val.block <- p.val.block[,2:5]
    ad.sphere.p.val.full <- rbind(ad.sphere.p.val.full, p.val.block[,1:4])
  }
}

alpha <- 0.05
n.vec <- colSums(!is.na(ad.sphere.p.val.full))
pow.vec <- colSums(ad.sphere.p.val.full <= alpha, na.rm = T)/n.vec
sd.vec <- sqrt((pow.vec)*(1 - pow.vec)/n.vec)

plot.dat <- data.frame("Clique_Size" = ell.set, 
                       "Power" = pow.vec, 
                       "SD" = sd.vec)

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Power)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = Power - 2*SD, ymax = Power + 2*SD)) + 
  ylim(0,.7) + 
  scale_x_continuous(breaks = pretty_breaks()) + 
  ggtitle(TeX(sprintf(r'(Adjacent Spheres Power C.C. Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size")


plt 

png(filename = "plots/ad_spheres_power.png",
    width = png.width, height = png.height, res = png.res)


plt

# Close the pdf file
dev.off()


#### Multiview Spheres Power #####
# we could update this later, add more simulations? 

ell.set <- c(6,8,12,16)

multi.p.val.full <- matrix(NA, nrow = 0, ncol = 4)

tri.const.filter = 1.6
for(block in 1:100){
  
  file <-  paste0("results/multiview_results_block_",block,".csv")
  if(file.exists(file)){
    p.val.block <- read.csv(file)
    rownames(p.val.block) = p.val.block[,1]
    p.val.block <- p.val.block[,2:5]
    multi.p.val.full <- rbind(multi.p.val.full, p.val.block[,1:4])
    
  }
  
}


alpha <- 0.05
n.vec <- colSums(!is.na(multi.p.val.full))
pow.vec <- colSums(multi.p.val.full <= alpha, na.rm = T)/n.vec
sd.vec <- sqrt((pow.vec)*(1 - pow.vec)/n.vec)

plot.dat <- data.frame("Clique_Size" = ell.set, 
                       "Power" = pow.vec, 
                       "SD" = sd.vec)

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Power)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = Power - 2*SD, ymax = Power + 2*SD)) + 
  ylim(0,.9) + 
  scale_x_continuous(breaks = pretty_breaks()) + 
  ggtitle(TeX(sprintf(r'(Multi View Power of Constant Curvature Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size")


plt 

png(filename = "plots/multiview_power.png",
    width = png.width, height = png.height, res = png.res)


plt

# Close the pdf file
dev.off()



#### Changepoint Consistency #####
# Here we are ensuring consistency of the estimated 
# Changepoints 


ell.set <- c(4,6,8,10,12)

changepoint.rmse <- matrix(NA, nrow = 0, ncol = 5)

tri.const.filter = 1.6
for(block in 1:100){
  
  file <-  paste0("results/changepoint_results_block_",block,".csv")
  if(file.exists(file)){
    rmse.block <- read.csv(file)
    rownames(rmse.block) = rmse.block[,1]
    rmse.block <- rmse.block[,2:6]
    changepoint.rmse <- rbind(changepoint.rmse, rmse.block[,1:5])
    
  }
  
}


n.vec <- colSums(!is.na(changepoint.rmse))
mean.vec <- colMeans(changepoint.rmse, na.rm = T)
median.vec <- mean.vec
median.vec[1] <- quantile(changepoint.rmse[,1], 0.5)
median.vec[2] <- quantile(changepoint.rmse[,2], 0.5)
median.vec[3] <- quantile(changepoint.rmse[,3], 0.5)
median.vec[4] <- quantile(changepoint.rmse[,4], 0.5)
median.vec[5] <- quantile(changepoint.rmse[,5], 0.5)

sd.vec <- colSDs(changepoint.rmse)/sqrt(n.vec)

plot.dat <- data.frame("Clique_Size" = ell.set, 
                       "Mean_RMSE" = mean.vec, 
                       "Median_RMSE" = median.vec,
                       "SD" = sd.vec)

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Median_RMSE), group = "Median", color = "black") +
  geom_line() + 
  geom_errorbar(aes(ymin = Median_RMSE - 2*SD, ymax = Median_RMSE + 2*SD)) + 
  ylim(-0.1,.9) + 
  #geom_line(aes(x = Clique_Size, y = Mean_RMSE), group = "Mean", color = "red") + 
  scale_x_continuous(breaks = pretty_breaks()) + 
  ggtitle(TeX(sprintf(r'(Consistency of Changepoints: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size") + 
  ylab("Median RMSE")



plt 

png(filename = "plots/Changepoint_RMSE.png",
    width = png.width, height = png.height, res = png.res)


plt

# Close the pdf file
dev.off()






###### Graph Statistics 


kappa.set <- c(-2,-1,-0.5,0,0.5,1)
#estimates 
scale.set <- c(0.7,1,2,4)


list.of.graph.stats.means <- list()
list.of.graph.stats.sds <- list()
for(kappa.idx in seq(length(kappa.set))){
  kappa <- kappa.set[kappa.idx]
  full.graph.stats.means <- matrix(NA, nrow = 6, ncol = 4)
  full.graph.stats.sds <- matrix(NA, nrow = 6, ncol = 4)
  for(scale.idx in seq(length(scale.set))){
    scale = scale.set[scale.idx]
    scale.full.graph.stats <- matrix(NA, nrow = 0, ncol = 7)
    for(block in 1:20){
      file <-  paste0("results/graph_stats_kappa_",kappa,"_scale_",scale,"_block_",block,".csv")
      if(file.exists(file)){
        stat.block <- read.csv(file)
        rownames(stat.block) = stat.block[,1]
        stat.block <- stat.block[,3:8]
        scale.full.graph.stats <- rbind(scale.full.graph.stats, cbind(stat.block, rep(scale,nrow(stat.block))))
        
      }
    }
    colnames(scale.full.graph.stats)[7] = "scale"
    full.graph.stats.means[,scale.idx] <- colMeans(scale.full.graph.stats[,1:6], na.rm = T)
    full.graph.stats.sds[,scale.idx] <- colSDs(scale.full.graph.stats[,1:6])
  }
  list.of.graph.stats.means[[kappa.idx]] <- full.graph.stats.means
  list.of.graph.stats.sds[[kappa.idx]] <- full.graph.stats.sds
}

library(stargazer)
#kappa = -2
tmp <- round(list.of.graph.stats.means[[1]], 3)
tmp2 <- round(list.of.graph.stats.sds[[1]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)

#kappa = -1
tmp <- round(list.of.graph.stats.means[[2]], 3)
tmp2 <- round(list.of.graph.stats.sds[[2]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)


#kappa = -0.5
tmp <- round(list.of.graph.stats.means[[3]], 3)
tmp2 <- round(list.of.graph.stats.sds[[3]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)


#kappa = 0
tmp <- round(list.of.graph.stats.means[[4]], 3)
tmp2 <- round(list.of.graph.stats.sds[[4]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)



#kappa = 0.5
tmp <- round(list.of.graph.stats.means[[5]], 3)
tmp2 <- round(list.of.graph.stats.sds[[5]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)


#kappa = 1
tmp <- round(list.of.graph.stats.means[[6]], 3)
tmp2 <- round(list.of.graph.stats.sds[[6]], 3)


tmp.out <- matrix("", nrow = 6,ncol = 4)
for(i in 1:6){
  for(j in 1:4){
    tmp.out[i,j] = paste0(as.character(tmp[i,j])," (", as.character(tmp2[i,j]), ")")
  }
}

colnames(tmp.out) <- scale.set
rownames(tmp.out) <- colnames(stat.block)
stargazer(tmp.out, summary = F)



