source("00_functions.R")

library(ggplot2)
library(dplyr)
library(tidyr)

png.width = 1500
png.height = 1500
png.res = 300

#Colorblind Friendly Palette 
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#CC79A7", "#F0E442", "#0072B2", "#D55E00")

## Estimating Equations Plot 
grid.size <- 10000
kappa.seq <- seq(-8,3, length.out = grid.size)

dxy <- 1
dxz <- 1
dyz <- 1

dxmM2 <- d_xm(-2,dxy,dxz,dyz)
dxmM1 <- d_xm(-1,dxy,dxz,dyz)
dxm0 <- d_xm(0,dxy,dxz,dyz)
dxm1 <- d_xm(1,dxy,dxz,dyz)
dxm2 <- d_xm(2,dxy,dxz,dyz)


plot.dat <- data.frame(matrix(NA, nrow = grid.size, ncol = 7))

# Spherical Schoenberg
D <- matrix(0,nrow = 4, ncol = 4)
D[1,2] = dxy
D[1,3] = dxz 
D[1,4] = dxm1
D[2,3] = dyz
D[2,4] = dyz/2
D[3,4] = dyz/2

D <- D + t(D)
for(kappa.idx in seq(grid.size)){
  kap <- kappa.seq[kappa.idx]
  
  plot.dat[kappa.idx,2] <- g_ee(kap,dxy,dxz,dyz,dxmM2)
  plot.dat[kappa.idx,3] <- g_ee(kap,dxy,dxz,dyz,dxmM1)
  plot.dat[kappa.idx,4] <- g_ee(kap,dxy,dxz,dyz,dxm0)
  plot.dat[kappa.idx,5] <- g_ee(kap,dxy,dxz,dyz,dxm1)
  plot.dat[kappa.idx,6] <- g_ee(kap,dxy,dxz,dyz,dxm2)
  
  C = cos(D * sqrt(kap + 0i))
  C1 <- C/(kap + 10^(-20))
  eig <- eigen(C1, symmetric = T)
  if(kap < -0.001){
    plot.dat[kappa.idx,7] <- 0.5*abs(eig$values[2])
  } else if(kap > 0.001){
    plot.dat[kappa.idx,7] <- 100*abs(eig$values[4])
  }
  
  
}

#plot.dat[,7] <- (0.1)*plot.dat[,7]
#plot.dat[which(plot.dat[,1] >= 0),7] <- 1000*plot.dat[which(plot.dat[,1] >= 0),7] 
colnames(plot.dat) <- c("kappa", 
                        "-2",
                        "-1",
                        "0",
                        "1",
                        "2", 
                        "Schoenberg: 1")
plot.dat[,1] <- kappa.seq

plot.dat.long <- plot.dat %>%
  pivot_longer(!kappa, names_to = "Curvature", values_to = "ee")


plot.dat.long$Curvature <- factor(plot.dat.long$Curvature, ordered = TRUE, 
                                  levels = c("-2",
                                             "-1",
                                             "0",
                                             "1",
                                             "2", 
                                             "Schoenberg: 1"))


plot.dat.long <- plot.dat.long[plot.dat.long$Curvature != "Schoenberg: 1", ]

plt <- ggplot(plot.dat.long, aes(x = kappa, y = ee, color = Curvature)) + 
  geom_line() + 
  ggtitle("Implicit Estimating Function") + 
  xlab("Kappa") + 
  ylab("g Value") + 
  ylim(-1,0.5) + 
  scale_color_manual(values = cbp1[1:5])
plt

scale.factor <- 1.5
png(filename = "plots/estimating_equations_example.png",
    width = scale.factor*png.width, height = scale.factor*png.height, res = (scale.factor)^2*png.res)


plt

# Close the pdf file
dev.off()




# 
library(rgl)
library(expm)
if (requireNamespace("MASS")) {
  Sigma <- matrix(c(10, -3, 0, -3, 5, 0, 0, 0, 1), 3, 3)
  Mean <- rep(0,3)
  x <- MASS::mvrnorm(500, Mean, Sigma)
  x <- x/sqrt(rowSums(x^2))
  x <- x %*% sqrtm(Sigma)
  #x <- x %*% Sigma
  open3d()
  
  plot3d(x, box = FALSE, 
         xlab = "x", 
         ylab = "y", 
         zlab = "z",
         xlim = c(-3,3), 
         ylim = c(-3,3), 
         zlim = c(-3,3))
  
  plot3d( ellipse3d((1/8)*Sigma, centre = Mean), 
          col = "green", alpha = 0.5, add = TRUE)
  rgl.postscript('plots/3dplotEllipse.eps', fmt = 'eps')
}  


# Plot the estimate and joint 90% confidence region for the displacement and cylinder
# count linear coefficients in the mtcars dataset



if (requireNamespace("MASS")) {
  Sigma <- matrix(c(5, 0, 0, 0, 5, 0, 0, 0, 5), 3, 3)
  Mean <- rep(0,3)
  x <- MASS::mvrnorm(500, Mean, Sigma)
  x <- x/sqrt(rowSums(x^2))
  x <- x %*% sqrtm(Sigma)
  #x <- x %*% Sigma
  open3d()
  
  plot3d(x, box = FALSE, 
         xlab = "x", 
         ylab = "y", 
         zlab = "z",
         xlim = c(-3,3), 
         ylim = c(-3,3), 
         zlim = c(-3,3))
  
  plot3d( ellipse3d((1/8)*Sigma, centre = Mean), col = "blue", alpha = 0.5, add = TRUE)
  rgl.postscript('plots/3dplotSphere.eps', fmt = 'eps')
}  









