


rm(list = ls())

source("R/00_functions.R")
#source("clique_finder.R")
set.seed(1)
library(latex2exp)
library(scales)
library(ggpubr)
### Convergence Results


save.plots = T
fig.height = 1500
fig.width = 1500
fig.res = 350


kappa.set <- c(-2,-1,-0.5,0,0.5,1)
#estimates
ell.set <- c(6,8,12,16) #c(6,8,12,16)
n.l = length(ell.set)

bad.ids = c() #TODO when some issues occur

plot.dat <- matrix(NA, 0,3)
tri.const.filter <- 1.55

trim.neginf = -Inf

for(tri.const.filter in seq(1,2, length.out = 21)) {
  plot.dat <- matrix(NA, 0,3)
  for(kappa in kappa.set){
    full.estimates <- matrix(NA, nrow = 0, ncol = n.l + 1)
    for(block in 1:125){
      id = which(kappa == kappa.set) + length(kappa.set)*(block - 1)
      if(id %in% bad.ids){
        next
      }
      file <-  paste0("results/estimates_kappa_",kappa,"_block_",block,".csv")
      if(file.exists(file)){
        est.block <- read.csv(file)
        if(trim.neginf){
          est.block[est.block < trim.neginf] = NA
        }
        if(ncol(est.block) != n.l + 2){
          next
        }
        rownames(est.block) = est.block[,1]
        est.block <- est.block[,2:(n.l + 2)]
        if(trim.neginf){
          est.block[est.block < trim.neginf] = NA
        }
        est.block <- est.block[est.block[,n.l + 1] == tri.const.filter, ]
        full.estimates <- rbind(full.estimates, est.block[,1:n.l])
      }

      ell.seq <- rep(ell.set, each = nrow(full.estimates))

      dat.tmp = matrix(c(ell.seq, as.numeric(unlist(full.estimates))), ncol = 2, nrow = n.l*nrow(full.estimates))
      dat.tmp <- cbind(dat.tmp, rep(kappa, nrow(full.estimates)))
      dat.tmp[,2] <- dat.tmp[,2] - kappa

    }
    plot.dat <- rbind(plot.dat, dat.tmp)
  }

  plot.dat <- as.data.frame(plot.dat)

  colnames(plot.dat) = c("CliqueSize","Bias", "Curvature")
  plot.dat$Curvature = as.factor(plot.dat$Curvature)
  plot.summary.data <- plot.dat %>% group_by(CliqueSize, Curvature) %>%
    summarize(q9 = quantile(Bias, 0.95, na.rm = T),
              q1  = quantile(Bias, 0.05, na.rm = T),
              med = median(Bias, na.rm = T),
              mean = mean(Bias, na.rm = T))

  # x.subset <- plot.dat %>% filter(CliqueSize == 12, Curvature == 0) %>% select(Bias)
  # plot(density(x.subset$Bias, na.rm = T,n = 5000), xlim = c(-2,2))

  plt <- ggplot(data = plot.summary.data, aes(x = CliqueSize, y = med, color = Curvature)) +
    geom_line() +
    geom_errorbar(aes(ymax = q9, ymin = q1),position=position_dodge(width=0.5)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #ggtitle(paste0(TeX("Consistency of Estimator: $C_{\\Delta}=$"), tri.const)) +
    ggtitle(TeX(sprintf(r'(Consistency of Estimator: $C_{\Delta}=%f$)', tri.const.filter))) +
    ylab("Estimate Deviation")+
    xlab("Clique Size") +
    theme_bw() +
    scale_x_continuous(breaks = ell.set,labels = ell.set)

  #dev.off()
  plot(plt)
  #Sys.sleep(3)
  #dev.off()
}


# remove the giant outliers.
trim.outlier = 100
num.total.points = 0
num.removed.points = 0

tri.const.filter = 1.5
plot.dat <- matrix(NA, 0,3)
for(kappa in kappa.set){
  full.estimates <- matrix(NA, nrow = 0, ncol = 3)
  for(block in 1:125){
    id = which(kappa == kappa.set) + length(kappa.set)*(block - 1)
    if(id %in% bad.ids){
      next
    }
    file <-  paste0("results/estimates_kappa_",kappa,"_block_",block,".csv")
    if(file.exists(file)){

      est.block <- read.csv(file)
      if(ncol(est.block) != n.l + 2){
        next
      }
      rownames(est.block) = est.block[,1]
      est.block <- est.block[,2:(n.l + 2)]
      est.block <- est.block[est.block[,n.l + 1] == tri.const.filter, ]
      full.estimates <- rbind(full.estimates, est.block[,1:n.l])


    }

    ell.seq <- rep(ell.set, each = nrow(full.estimates))

    dat.tmp = matrix(c(ell.seq, as.numeric(unlist(full.estimates))), ncol = 2, nrow = n.l*nrow(full.estimates))
    dat.tmp <- cbind(dat.tmp, rep(kappa, nrow(full.estimates)))
    dat.tmp[,2] <- dat.tmp[,2] - kappa
  }
  plot.dat <- rbind(plot.dat, dat.tmp)
}

plot.dat <- as.data.frame(plot.dat)

colnames(plot.dat) = c("CliqueSize","Bias", "Curvature")
plot.dat$nooutlier = 1*(abs(plot.dat$Bias) < trim.outlier)
plot.dat$Curvature = as.factor(plot.dat$Curvature)
plot.summary.data <- plot.dat %>% group_by(CliqueSize, Curvature) %>%
  summarize(q9 = quantile(Bias, 0.95, na.rm = T),
            q1  = quantile(Bias, 0.05, na.rm = T),
            med = median(Bias, na.rm = T),
            mean = mean(Bias*nooutlier, na.rm = T),
            mad = mad(Bias, na.rm = T))

# x.subset <- plot.dat %>% filter(CliqueSize == 12, Curvature == 0) %>% select(Bias)
# plot(density(x.subset$Bias, na.rm = T,n = 5000), xlim = c(-2,2))

plt <- ggplot(data = plot.summary.data, aes(x = CliqueSize, y = med, color = Curvature)) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymax = q9, ymin = q1),position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #ggtitle(paste0(TeX("Consistency of Estimator: $C_{\\Delta}=$"), tri.const)) +
  ggtitle(TeX(sprintf(r'(Consistency of Estimator: $C_{\Delta}=%f$)', tri.const.filter))) +
  ylab("Estimate Deviation")+
  xlab("Clique Size") +
  theme_bw() +
  scale_x_continuous(breaks = ell.set,labels = ell.set)

plt <- ggplot(data = plot.summary.data, aes(x = CliqueSize, y = mean, color = Curvature)) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymax = q9, ymin = q1),position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #ggtitle(paste0(TeX("Consistency of Estimator: $C_{\\Delta}=$"), tri.const)) +
  ggtitle(TeX(sprintf(r'(Consistency of Estimator: $C_{\Delta}=%f$)', tri.const.filter))) +
  ylab("Estimate Deviation")+
  xlab("Clique Size") +
  theme_bw() +
  scale_x_continuous(breaks = ell.set,labels = ell.set)

dev.off()
plt
print(paste0("Fraction of filtered Estimates: ", round(mean(1 - plot.dat$nooutlier, na.rm = T), 5)))

if(save.plots){
  plt %>%
    ggexport(filename = "plots/GMM_consistency.png", width = fig.width, height = fig.height, res = fig.res)
}




# png(filename = "plots/GMM_consistency.png",
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
ell.set <- c(6,8,12)#c(6,8,12,16)

p.val.plot.dat <- matrix(NA, 0,3)
tri.const.filter <- 1.2
for(tri.const.filter in seq(1,2, length.out = 21)) {
  p.val.plot.dat <- matrix(NA, 0,3)
  for(kappa in kappa.set){

    p.val.full <- matrix(NA, nrow = 0, ncol = 4)
    for(block in 1:125){
      id = which(kappa == kappa.set) + length(kappa.set)*(block - 1)
      if(id %in% bad.ids){
        next
      }
      #file <-  paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
      file <-  paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
      if(file.exists(file)){
        p.val.block <- read.csv(file)
        rownames(p.val.block) = p.val.block[,1]
        p.val.block <- p.val.block[,2:5]
        p.val.block <- p.val.block[p.val.block[,4] == tri.const.filter, ]
        p.val.full <- rbind(p.val.full, p.val.block[,1:3])
        ell.seq <- rep(ell.set, each = nrow(p.val.full))
        dat.tmp = matrix(c(ell.seq, as.numeric(unlist(p.val.full))), ncol = 2, nrow = 3*nrow(p.val.full))
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
    geom_line(position=position_dodge(width=0.5)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    # geom_errorbar(aes(ymax = fpr05 + 2*sqrt(var.fpr05),
    #                   ymin = fpr05 - 2*sqrt(var.fpr05))) +
    ggtitle(TeX(sprintf(r'(FPR of Test: $C_{\Delta}=%f$)', tri.const.filter))) +
    ylab("False Positive Rate") +
    xlab("Clique Size") +
    ylim(c(0,1)) +
    theme_bw() +
    scale_x_continuous(breaks = ell.set,labels = ell.set)

  plot(plt)
}



tri.const.filter = 1.2
p.val.plot.dat <- matrix(NA, 0,3)
for(kappa in kappa.set){

  p.val.full <- matrix(NA, nrow = 0, ncol = 4)
  for(block in 1:125){
    id = which(kappa == kappa.set) + length(kappa.set)*(block - 1)
    if(id %in% bad.ids){
      next
    }
    #file <-  paste0("results/p_vals_kappa_",kappa,"_block_",block,".csv")
    file <-  paste0("results/true_model_R3_rho_0p5/p_vals_kappa_",kappa,"_block_",block,".csv")
    if(file.exists(file)){
      p.val.block <- read.csv(file)
      rownames(p.val.block) = p.val.block[,1]
      p.val.block <- p.val.block[,2:5]
      p.val.block <- p.val.block[p.val.block[,4] == tri.const.filter, ]
      p.val.full <- rbind(p.val.full, p.val.block[,1:3])
      ell.seq <- rep(ell.set, each = nrow(p.val.full))
      dat.tmp = matrix(c(ell.seq, as.numeric(unlist(p.val.full))), ncol = 2, nrow = 3*nrow(p.val.full))
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
            var.fpr1 = (1/(sum(pval >= 0.0, na.rm = T)))*mean(pval < 0.1, na.rm = T)*(1 - mean(pval < 0.1, na.rm = T)),
            var.fpr05 = (1/(sum(pval >= 0.0, na.rm = T)))*mean(pval < 0.05, na.rm = T)*(1 - mean(pval < 0.05, na.rm = T)),
            var.fpr01  = (1/(sum(pval >= 0.0, na.rm = T)))*mean(pval < 0.01, na.rm = T)*(1 - mean(pval < 0.01, na.rm = T)),
            var.fpr005  = (1/(sum(pval >= 0.0, na.rm = T)))*mean(pval < 0.005, na.rm = T)*(1 - mean(pval < 0.005, na.rm = T)))

ggplot(data = p.val.plot.dat,
       aes(x = CliqueSize, y = pval, color = Curvature)) +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  ggtitle("Coverage Of Test") +
  ylab("Estimate Deviation")


plt <- ggplot(data = p.val.plot.summary.data,
              aes(x = CliqueSize, y = fpr05, color = Curvature)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_errorbar(aes(ymax = fpr05 + 2*sqrt(var.fpr05),
                    ymin = fpr05 - 2*sqrt(var.fpr05)),position=position_dodge(width=0.5)) +
  ggtitle(TeX(sprintf(r'(FPR of Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  ylab("False Positive Rate") +
  xlab("Clique Size") +
  ylim(c(0,.5)) +
  theme_bw() +
  scale_x_continuous(breaks = ell.set,labels = ell.set)

plot(plt)


# tmp <- p.val.plot.dat %>% filter(Curvature == 0, CliqueSize == 12) %>%  select(pval)

if(save.plots){
  plt %>%
    ggexport(filename = "plots/GMM_fpr_cc_test.png", width = fig.width, height = fig.height, res = fig.res)
}




#### Adjacent Spheres Power #####
# we could update this later, add more simulations?

ell.set <- c(6,8,12)#c(6,8,12,16)

ad.sphere.p.val.full <- matrix(NA, nrow = 0, ncol = 3)

tri.const.filter = 1.2
for(block in 1:20){

  file <-  paste0("results/adjacent_spheres_results_block_",block,".csv")
  if(file.exists(file)){
    p.val.block <- read.csv(file)
    rownames(p.val.block) = p.val.block[,1]
    p.val.block <- p.val.block[,2:4]
    ad.sphere.p.val.full <- rbind(ad.sphere.p.val.full, p.val.block[,1:3])
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
  theme_bw() +
  ggtitle(TeX(sprintf(r'(Adjacent Spheres Power C.C. Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size")


plt

if(save.plots){
  plt %>%
    ggexport(filename = "plots/ad_spheres_power.png", width = fig.width, height = fig.height, res = fig.res)
}


#### Multiview Spheres Power #####
# we could update this later, add more simulations?

ell.set <- c(6,8,12)

multi.p.val.full <- matrix(NA, nrow = 0, ncol = 3)

tri.const.filter = 1.2
for(block in 1:20){

  file <-  paste0("results/multiview_results_block_",block,".csv")
  if(file.exists(file)){
    p.val.block <- read.csv(file)
    rownames(p.val.block) = p.val.block[,1]
    p.val.block <- p.val.block[,2:4]
    multi.p.val.full <- rbind(multi.p.val.full, p.val.block[,1:3])

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
  theme_bw() +
  ggtitle(TeX(sprintf(r'(Multi View Power of C.C. Test: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size")

dev.off()
plt

if(save.plots){
  plt %>%
    ggexport(filename = "plots/multiview_power.png", width = fig.width, height = fig.height, res = fig.res)
}




#### Changepoint Consistency #####
# Here we are ensuring consistency of the estimated
# Changepoints


ell.set <- c(4,6,8,10,12)

changepoint.mad <- matrix(NA, nrow = 0, ncol = 5)

tri.const.filter = 1.4
for(block in 1:100){

  file <-  paste0("results/changepoint_results_block_",block,".csv")
  if(file.exists(file)){
    mad.block <- read.csv(file)
    rownames(mad.block) = mad.block[,1]
    mad.block <- mad.block[,2:6]
    changepoint.mad <- rbind(changepoint.mad, mad.block[,1:5])

  }

}


n.vec <- colSums(!is.na(changepoint.mad))
mean.vec <- colMeans(changepoint.mad, na.rm = T)
median.vec <- mean.vec
median.vec[1] <- quantile(changepoint.mad[,1], 0.5)
median.vec[2] <- quantile(changepoint.mad[,2], 0.5)
median.vec[3] <- quantile(changepoint.mad[,3], 0.5)
median.vec[4] <- quantile(changepoint.mad[,4], 0.5)
median.vec[5] <- quantile(changepoint.mad[,5], 0.5)

sd.vec <- colSDs(changepoint.mad)/sqrt(n.vec)

plot.dat <- data.frame("Clique_Size" = ell.set,
                       "Mean_MAD" = mean.vec,
                       "Median_MAD" = median.vec,
                       "SD" = sd.vec)

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Median_MAD), group = "Median", color = "black") +
  geom_line() +
  geom_errorbar(aes(ymin = Median_MAD - 2*SD, ymax = Median_MAD + 2*SD)) +
  ylim(-0.1,.9) +
  #geom_line(aes(x = Clique_Size, y = Mean_MAD), group = "Mean", color = "red") +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggtitle(TeX(sprintf(r'(Consistency of Changepoints: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size") +
  ylab("Median MAD")

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Mean_MAD), group = "Median", color = "black") +
  geom_line() +
  geom_errorbar(aes(ymin = Mean_MAD - 2*SD, ymax = Mean_MAD + 2*SD)) +
  ylim(-0.1,.9) +
  #geom_line(aes(x = Clique_Size, y = Mean_MAD), group = "Mean", color = "red") +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggtitle(TeX(sprintf(r'(Consistency of Changepoints: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size") +
  theme_bw() +
  ylab("Mean AD")

dev.off()
plt


if(save.plots){
  plt %>%
    ggexport(filename = "plots/changepoint_MeanAD.png", width = fig.width, height = fig.height, res = fig.res)
}




ell.set <- c(4,6,8,10,12)

changepoint.mad <- matrix(NA, nrow = 0, ncol = 5)

tri.const.filter = 1.4
for(block in 1:100){

  file <-  paste0("results/changepoint_results_median_block_",block,".csv")
  if(file.exists(file)){
    mad.block <- read.csv(file)
    rownames(mad.block) = mad.block[,1]
    mad.block <- mad.block[,2:6]
    changepoint.mad <- rbind(changepoint.mad, mad.block[,1:5])

  }

}


n.vec <- colSums(!is.na(changepoint.mad))
mean.vec <- colMeans(changepoint.mad, na.rm = T)
median.vec <- mean.vec
median.vec[1] <- quantile(changepoint.mad[,1], 0.5)
median.vec[2] <- quantile(changepoint.mad[,2], 0.5)
median.vec[3] <- quantile(changepoint.mad[,3], 0.5)
median.vec[4] <- quantile(changepoint.mad[,4], 0.5)
median.vec[5] <- quantile(changepoint.mad[,5], 0.5)

sd.vec <- colSDs(changepoint.mad)/sqrt(n.vec)

plot.dat <- data.frame("Clique_Size" = ell.set,
                       "Mean_MAD" = mean.vec,
                       "Median_MAD" = median.vec,
                       "SD" = sd.vec)

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Median_MAD), group = "Median", color = "black") +
  geom_line() +
  geom_errorbar(aes(ymin = Median_MAD - 2*SD, ymax = Median_MAD + 2*SD)) +
  ylim(-0.1,.9) +
  #geom_line(aes(x = Clique_Size, y = Mean_MAD), group = "Mean", color = "red") +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggtitle(TeX(sprintf(r'(Consistency of Changepoints: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size") +
  ylab("Median AD")

plt <- ggplot(plot.dat, aes(x = Clique_Size, y = Mean_MAD), group = "Median", color = "black") +
  geom_line() +
  geom_errorbar(aes(ymin = Mean_MAD - 2*SD, ymax = Mean_MAD + 2*SD)) +
  ylim(-0.1,.9) +
  #geom_line(aes(x = Clique_Size, y = Mean_MAD), group = "Mean", color = "red") +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggtitle(TeX(sprintf(r'(Consistency of Changepoints: $C_{\Delta}=%f$)', tri.const.filter))) +
  xlab("Clique Size") +
  theme_bw() +
  ylab("Median MAD")

dev.off()
plt


if(save.plots){
  plt %>%
    ggexport(filename = "plots/changepoint_MedianAD.png", width = fig.width, height = fig.height, res = fig.res)
}




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




