
#library(mvtnorm)
# A warning, sometimes mvtnorm just breaks

#library(truncnorm)
library(igraph)
library(dplyr)
#library(Rfast)
library(ggplot2)
library(CVXR)
#
library(truncnorm)
library(Matrix)
#require(mvtnorm)
library(igraph)
library(parallel)
library(Rfast)
library(poweRlaw)
library(LaplacesDemon)
library(sClust)
library(RColorBrewer)
library(colorspace)

sec <- function(x){
  1/cos(x)
}

expit <- function(x){
  1/(1 + exp(-x))
}

# fast computation of the elementary symmetric polynomial
# Computing the elementary symmetric polynomial is hard to do directly
# this simplified version makes it much more feasible
# if t gets larger than say 400, then this can quickly become numerically unstable
es_poly <- function(x){
  t = length(x)
  prev.polynomial = c(1)
  for(s in 1:t){
    next.polynomial = c(0,prev.polynomial)*x[s] + c(prev.polynomial,0)
    prev.polynomial = next.polynomial
  }
  return(next.polynomial)
}


# estimating equation for kappa
g_ee <- function(kappa,dxy,dxz,dyz,dxm){
  # using a single value to construct the estimating function use 0i to compute using complex numbers.
  thresh = 10**(-13)
  if(abs(kappa)  > thresh){
    g = (2*cos(dxm*sqrt(kappa + 0i)) - sec(dyz/2*sqrt(kappa + 0i))*(cos(dxy*sqrt(kappa + 0i)) + cos(dxz*sqrt(kappa + 0i))))/(kappa)
    g = Re(g)
  } else if(abs(kappa)  <= thresh){
    g = ((dxy**2 + dxz**2)/2) - dyz**2/4 - dxm**2
  }
  return(g)
}

# true distance as a function of dxy, dxz, dyz
d_xm <- function(kappa,dxy,dxz,dyz){
  out.com <- (1/sqrt(kappa + 0i))*acos((1/2)*sec(dyz/2*sqrt(kappa + 0i))*(cos(dxy*sqrt(kappa + 0i)) + cos(dxz*sqrt(kappa + 0i))))
  if(kappa > 0){
    out <- Re(out.com)
  } else if(kappa < 0){
    out <- Re(out.com)
  } else {
    out <- sqrt((dxy**2 + dxz**2)/2 - dyz**2/4)
  }
  return(out)
}



# derivative of the estimating function.
g_grad_kappa <- function(kappa,dxy,dxz,dyz,dxm, manual.grad = F, eps = 10**(-9)){

  # log-scale computing for small kappa
  if(manual.grad){
    gp = (lolaR::gEF(kappa + eps,dxy,dxz,dyz,dxm) - lolaR::gEF(kappa + eps,dxy,dxz,dyz,dxm))/eps
  } else{
    if(kappa != 0){
      gp = (-1)*(1/(4*kappa^2))*(8*cos(dxm*sqrt(kappa + 0i)) - 4*cos(dxz*sqrt(kappa + 0i))*sec((dyz*sqrt(kappa + 0i))/2) +
                                   4*dxm*sqrt(kappa + 0i)*sin(dxm*sqrt(kappa + 0i)) - 2*dxy*sqrt(kappa + 0i)*sec((dyz*sqrt(kappa + 0i))/2)*sin(dxy*sqrt(kappa + 0i)) -
                                   2*dxz*sqrt(kappa + 0i)*sec((dyz*sqrt(kappa + 0i))/2)*sin(dxz*sqrt(kappa + 0i)) +
                                   dyz*sqrt(kappa + 0i)*cos(dxz*sqrt(kappa + 0i))*sec((dyz*sqrt(kappa + 0i))/2)*tan((dyz*sqrt(kappa + 0i))/2) +
                                   cos(dxy*sqrt(kappa + 0i))*sec(dyz*sqrt(kappa + 0i)/2)*(-4 + dyz*sqrt(kappa + 0i)*tan(dyz*sqrt(kappa + 0i)/2)))
      gp = Re(gp)
    } else if(kappa == 0){
      gp = (-1)*(1/192)*(-16*dxm^4 + 8*dxy^4 + 8*dxz^4 - 12*dxy^2*dyz^2 - 12*dxz^2*dyz^2 + 5*dyz^4)
    }
  }
  return(gp)
}

# gradient with respect to d
g_grad_d <- function(kappa, dxy, dxz, dyz, dxm, manual.negative = F, eps = 10**(-8)){
  if(kappa > 0){
    g1 <- (1/(sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*sin(dxy*sqrt(kappa + 0i))
    g2 <- (1/(sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*sin(dxz*sqrt(kappa + 0i))
    g3 <- (1/(2*sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*tan(dyz/2*sqrt(kappa + 0i))*(cos(dxy*sqrt(kappa + 0i)) + cos(dxz*sqrt(kappa + 0i)))
    g4 <- -(2/(sqrt(kappa + 0i)))*sin(dxm*sqrt(kappa + 0i))
    g.vec <- c(g1,g2,g3,g4)
  }
  else if(kappa == 0){
    g1 <- dxy + 0i
    g2 <- dxz + 0i
    g3 <- dyz/2 + 0i
    g4 <- -2*dxm + 0i
    g.vec <- c(g1,g2,g3,g4)
  } else {
    g1 <- (1/(sqrt( -kappa + 0i)))*sech(dyz/2*sqrt(-kappa + 0i))*sinh(dxy*sqrt(-kappa + 0i))
    g2 <- (1/(sqrt(-kappa + 0i)))*sech(dyz/2*sqrt(-kappa + 0i))*sinh(dxz*sqrt(-kappa + 0i))
    g3 <- -(1/(2*sqrt(-kappa + 0i)))*sech(dyz/2*sqrt(-kappa + 0i))*tanh(dyz/2*sqrt(-kappa + 0i))*(cosh(dxy*sqrt(-kappa + 0i)) + cosh(dxz*sqrt(-kappa + 0i)))
    g4 <- (2/(sqrt(-kappa + 0i)))*sinh(dxm*sqrt(-kappa + 0i))

    if(manual.negative){
      g1 <- (lolaR::gEF(kappa, dxy + eps, dxz, dyz, dxm) - lolaR::gEF(kappa, dxy, dxz, dyz, dxm))/eps
      g2 <- (lolaR::gEF(kappa, dxy, dxz + eps, dyz, dxm) - lolaR::gEF(kappa, dxy, dxz, dyz, dxm))/eps
      g3 <- (lolaR::gEF(kappa, dxy, dxz, dyz + eps, dxm) - lolaR::gEF(kappa, dxy, dxz, dyz, dxm))/eps
      g4 <- (lolaR::gEF(kappa, dxy, dxz, dyz, dxm + eps) - lolaR::gEF(kappa, dxy, dxz, dyz, dxm))/eps
    }

    g.vec <- c(g1,g2,g3,g4)
  }
  return(g.vec)
}

sech <- function(x){
  1/cosh(x)
}

#Upper bound on curvature Estimate
# g_u <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
#   out <- g_ee(kappa,dxy,dxz,dyz,dxm + d_xm(kappa,dym,dzm,dyz))
# }
#
# g_l <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
#   out <- g_ee(kappa,dxy,dxz,dyz,dxm - d_xm(kappa,dym,dzm,dyz))
# }
#
#
#
#
# kappa_u <- function(dxy, dxz, dyz, dxm, dym, dzm,
#                     kappa.prec = 10**(-5),
#                     min.curvature = -1000){
#
#   # first check for whether the midpoint estimate is already too far for any curvature:
#
#   if(dxm < min(dxy,dxz) - (1/2)*dyz){
#     warning("Triangle Inequality is not satisfied")
#     return(-Inf)
#   }
#
#
#   # Picking good initialization for the grid search.
#   max.curvature = (pi/max(c(dxy,dxz,dyz,dxm,dym,dzm)))**2
#
#   kappa.upper <- max.curvature
#   kappa.lower <- min.curvature
#
#   g.u.upper <- g_u(kappa.upper, dxy, dxz, dyz, dxm, dym, dzm)
#   g.u.lower <- g_u(kappa.lower, dxy, dxz, dyz, dxm, dym, dzm)
#
#
#   if(g.u.upper < 0 ) {
#     return(max.curvature)
#   } else if(g.u.lower > 0) {
#     return(min.curvature)
#   } else {
#     kappa.gap <- kappa.upper - kappa.lower
#     while(kappa.gap > kappa.prec){
#       kappa.mid <-  mean(c(kappa.upper, kappa.lower))
#       g.u.mid <- g_u(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
#       #print(kappa.mid)
#       # Handles a few numerical errors
#       tries <- 0
#       while(is.nan(g.u.mid) & tries < 10){
#         #print(tries)
#
#         kappa.mid <- runif(1,kappa.lower, kappa.upper)
#         #print(kappa.mid)
#         g.u.mid <- g_u(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
#         tries <- tries + 1
#       }
#       if(tries >= 10){
#         return(max.curvature)
#       }
#
#       if(g.u.mid > 0){
#         kappa.upper <- kappa.mid
#       } else if(g.u.mid == 0) {
#         break
#       } else {
#         kappa.lower <- kappa.mid
#       }
#       kappa.gap <- kappa.upper - kappa.lower
#     }
#     return(kappa.mid)
#   }
# }
#
# kappa_l <- function(dxy, dxz, dyz, dxm, dym, dzm,
#                     kappa.prec = 10**(-5),
#                     min.curvature = -1000){
#
#   # first check for whether the midpoint estimate is already too far for any curvature:
#
#   if(dxm < min(dxy,dxz) - (1/2)*dyz){
#     warning("Triangle Inequality is not satisfied")
#     return(-Inf)
#   }
#
#
#   # Picking good initialization for the grid search.
#   max.curvature = (pi/max(c(dxy,dxz,dyz,dxm,dym,dzm)))**2
#
#   kappa.upper <- max.curvature
#   kappa.lower <- min.curvature
#
#   g.l.upper <- g_l(kappa.upper, dxy, dxz, dyz, dxm, dym, dzm)
#   g.l.lower <- g_l(kappa.lower, dxy, dxz, dyz, dxm, dym, dzm)
#
#
#   if(g.l.upper < 0 ) {
#     return(max.curvature)
#   } else if(g.l.lower > 0) {
#     return(min.curvature)
#   } else {
#     kappa.gap <- kappa.upper - kappa.lower
#     while(kappa.gap > kappa.prec){
#       kappa.mid <- mean(c(kappa.upper, kappa.lower))
#
#       g.l.mid <- g_l(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
#
#       # Handles a few numerical errors
#       tries <- 0
#       while(is.nan(g.l.mid) & tries < 10){
#         kappa.mid <- runif(1,kappa.lower, kappa.upper)
#         g.l.mid <- g_l(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
#         tries <- tries + 1
#       }
#       if(tries >= 10){
#         return(min.curvature)
#       }
#
#       if(g.l.mid > 0){
#         kappa.upper <- kappa.mid
#       } else if(g.l.mid == 0) {
#         break
#       } else {
#         kappa.lower <- kappa.mid
#       }
#       kappa.gap <- kappa.upper - kappa.lower
#     }
#     return(kappa.mid)
#   }
# }
#
# estimateBounds <- function(D,y,z,m,x.set){
#   kappa.us <- c()
#   kappa.ls <- c()
#
#   for(x in x.set){
#
#     dxy <- D[x,y]
#     dxz <- D[x,z]
#     dyz <- D[y,z]
#     dxm <- D[x,m]
#     dym <- D[y,m]
#     dzm <- D[z,m]
#     set.seed(2)
#     kap.u <- kappa_u(dxy, dxz, dyz, dxm, dym, dzm,
#                      kappa.prec = 10**(-5),
#                      min.curvature = -1000)
#
#     kap.l <- kappa_l(dxy, dxz, dyz, dxm, dym, dzm,
#                      kappa.prec = 10**(-5),
#                      min.curvature = -1000)
#
#     kappa.us <- c(kappa.us,kap.u)
#     kappa.ls <- c(kappa.ls,kap.l)
#   }
#   return(list("upper.bounds" = kappa.us,
#               "lower.bounds" = kappa.ls))
# }
#
#
# ### Subsample
# SubSampleConstantCurvatureTest <- function(A,clique.set,reference.set,
#                                            subsample.rate = 1,B = 1000){
#
#
#
#   A.subset = A
#   K = length(clique.set)
#   for(b in seq(B)){
#     for(k in seq(K)){
#       ell = length(clique.set[[k]])
#       clique.subsample[[k]] <- sample(clique.set[[k]],
#                                       size = ell - subsample.rate,
#                                       replace = F)
#
#     }
#   }
# }
#
#
#
#
# # g1 <- c()
# # g2 <- c()
# #
# # kappa.seq <- seq(-5,0.5,length.out = 200)
# # for(kap in kappa.seq){
# #   g1 <- c(g1,g_u(kap, dxy, dxz, dyz, dxm, dym, dzm))
# #   g2 <- c(g2,g_l(kap, dxy, dxz, dyz, dxm, dym, dzm))
# # }
#
# # plot(kappa.seq, g1, type = "l", col = "red")
# # lines(kappa.seq, g2, col = "blue")
# #
# # plot(g1,g2, type = "l")
# # Newton method for estimating kappa.hat
# # y = 1
# # z = 2
# # m = 3
# # x = 4
# # dyz <- D[y,z]
# # dxy <- D[x,y]
# # dxz <- D[x,z]
# # dxm <- D[x,m]
# #
# # estimate_kappa(dxy,dxz,dyz,dxm)
#
# # dyz <- 1
# # dxy <- 1
# # dxz <- 1
# # dxm <- 1/2 + 0.01
# estimate_kappa <- function(dxy,dxz,dyz,dxm,
#                            kappa.init = 0, thresh = 10^(-6),
#                            max.iter = 10, max.curve = Inf,
#                            ee.thresh = 0.1, kappa.reset.step = 10){
#
#   # first check for whether the midpoint estimate is already too far for any curvature:
#   min.curvature = -5000
#   neg.init <- dxm <= d_xm(0,dxy,dxz,dyz)
#   max.curvature = (pi/max(c(dxy,dxz,dyz)))**2
#   kap.seq <- seq(0,max.curvature, 0.01)
#   d.xm.vec <- c()
#   for(kap in kap.seq){
#     suppressWarnings(d.xm.vec <- c(d.xm.vec, d_xm(kap,dxy,dxz,dyz)))
#   }
#   max.dxm = max(d.xm.vec, na.rm = T)
#
#   if(!neg.init){
#     idx <- which.min(abs(d.xm.vec - dxm) )
#     kappa.init <- kap.seq[idx]
#   }
#
#   kap.seq <- seq(min.curvature,0, 1)
#   d.xm.vec <- c()
#   for(kap in kap.seq){
#     suppressWarnings(d.xm.vec <- c(d.xm.vec, d_xm(kap,dxy,dxz,dyz)))
#   }
#   min.dxm = min(d.xm.vec, na.rm = T)
#
#   if(neg.init){
#     idx <- which.min(abs(d.xm.vec - dxm) )
#     kappa.init <- kap.seq[idx]
#   }
#
#   # cases for impossible midpoint distances
#   if(any(is.na(c(dxy,dxz,dyz)))){
#     return(NA)
#   } else if(dxm > max.dxm){
#     return(Inf)
#   } else if(dxm < min.dxm){
#     return(-Inf)
#   } else {
#     max.curvature = (pi/max(c(dxy,dxz,dyz,dxm)))**2
#
#     kappa.prev = kappa.init
#     diff = Inf
#
#     max.curvature = (pi/max(c(dxy,dxz,dyz,dxm)))**2 # maximal value in estimating equation
#
#     iter = 0
#     reset.last = F
#     was.reset = F
#     while(diff > thresh & iter < max.iter & abs(kappa.prev) < max.curve){
#       kappa.next = kappa.prev - g_ee(kappa.prev,dxy,dxz,dyz,dxm)/g_grad_kappa(kappa.prev,dxy,dxz,dyz,dxm)
#       diff = abs(kappa.next - kappa.prev)
#       # resetting step if jumping outside the feasible range:
#       # TODO: Check with Tyler if this is a reasonable optimization strategy
#       kappa.prev = kappa.next
#       # resets the curvature
#
#       # what happens when we had reset and then still run into the same problem
#       if(reset.last & was.reset){
#         iter = max.iter + 1
#         break
#       }
#       # restart looking for the new midpoint.
#       if(reset.last ){
#         kappa.prev = max.curvature - 10**(-1)
#         reset.last = F
#         edge.sign = sign(g_ee(kappa.prev,dxy,dxz,dyz,dxm))
#         kappa.reset = -kappa.reset.step
#         while((sign(g_ee(kappa.reset,dxy,dxz,dyz,dxm))) == edge.sign){
#           kappa.reset = kappa.reset-kappa.reset.step
#           if(kappa.reset < - 1000){
#             iter = max.iter + 1
#             break
#           }
#         }
#         kappa.prev = kappa.reset # reset the placement
#         iter = 0
#         was.reset = T
#       }
#
#       if(kappa.prev > max.curvature){
#         kappa.prev = max.curvature - 10**(-1)
#         reset.last = T # note that we went outside of the
#       }
#       iter = iter + 1
#       #print(kappa.prev)
#     }
#     kappa.est = kappa.next
#     if(iter > max.iter | abs(kappa.est) > max.curve | abs(g_ee(kappa.prev,dxy,dxz,dyz,dxm)) > ee.thresh) {
#       kappa.est = NA
#     }
#   }
#   return(kappa.est)
# }
#


# kappa.set <- seq(-1,1,0.01)
# ee <- c()
# for(kap in kappa.set){
#   ee <- c(ee,g_ee(kap,dxy,dxz,dyz,dxm))
# }
# plot(kappa.set,ee)

# kappa_variance <- function(kappa.est, d.vec, d.var, rho.correction = 0){
#   Sigma.hat <- diag(d.var)
#   #Sigma.hat[2,4] <- rho.correction
#   Sigma.hat[3,4] <- rho.correction
#   #Sigma.hat[4,2] <- Sigma.hat[2,4]
#   Sigma.hat[4,3] <- Sigma.hat[3,4]
#
#   dxy = d.vec[1]
#   dxz = d.vec[2]
#   dyz = d.vec[3]
#   dxm = d.vec[4]
#   dgd <- g_grad_d(kappa.est, dxy, dxz, dyz, dxm)
#   dgk <- g_grad_d(kappa.est, dxy, dxz, dyz, dxm)
#   kap.var.hat <- (t(dgd) %*% Sigma.hat %*% dgd)/(dgk**2)
#   return(kap.var.hat)
#
# }

#simulate latent space positions.
# flatness compresses the isotropic map into a approximately a 2d disk
sim_latent_pos <- function(sigma,p,kappa, n, flatness = 1){
  if(p > 2){
    if(kappa == 0){
      sigma.flat = c(rep(flatness,2)*sigma,rep(sigma, p - 2))
      Z = mvtnorm::rmvnorm(n, sigma = sigma.flat)
    } else if (kappa < 0){
      sigma.flat = c(rep(flatness,2)*sigma,rep(sigma, p - 2))
      Z = mvtnorm::rmvnorm(n, sigma = sigma.flat)
      Z.0 = apply(Z, 1, function(z){
        out = sum(z^2) + 1
      })
      Z = cbind(sqrt(Z.0),Z)
    } else if (kappa > 0){ #ends up being uniform in this parameterization
      sigma.flat = c(rep(flatness,2)*sigma,rep(sigma, p - 2))
      Z = mvtnorm::rmvnorm(n, sigma = sigma.flat)
      Z.norm = apply(Z, 1, function(z){
        out = sum(z^2)
      })
      # should we or should we not embed in a literally larger sphere
      #Z = (1/sqrt(kappa))*Z/sqrt(Z.norm)
      Z = Z/sqrt(Z.norm)
    }
  } else {
    if(kappa == 0){
      Z = mvtnorm::rmvnorm(n, sigma = diag(rep(sigma, p)))
    } else if (kappa < 0){
      Z = mvtnorm::rmvnorm(n, sigma = diag(rep(sigma, p)))
      Z.0 = apply(Z, 1, function(z){
        out = sum(z^2) + 1
      })
      Z = cbind(sqrt(Z.0),Z)
    } else if (kappa > 0){ #ends up being uniform in this parameterization
      Z = mvtnorm::rmvnorm(n, sigma = diag(rep(sigma, p + 1)))
      Z.norm = apply(Z, 1, function(z){
        out = sum(z^2)
      })
      # should we or should we not embed in a literally larger sphere
      #Z = (1/sqrt(kappa))*Z/sqrt(Z.norm)
      Z = Z/sqrt(Z.norm)
    }
  }

  return(Z)
}


## TODO: Correct this
# sampling from a uniform ball in various geometries
sim_latent_uniform_ball <- function(n,p,kappa,radius, flatness = 1){
  # equal angle sampler
  # if(perp){
  #   X <- mvtnorm::rmvnorm(n, mean = rep(0,p-1), sigma = diag(rep(1,p -1)))
  #   X <- cbind(rep(0,n),X)
  # } else {
  #   X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = diag(rep(1,p)))
  # }
  sigma = 1
  if(p > 2){
    sigma.flat = c(rep(flatness,2)*sigma,rep(sigma, p - 2))
    X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = diag(rep(1,p)))
    X <- cbind(rep(0,n),X)
  }


  X = X/sqrt(rowSums(X**2)) # direction X
  R = radial_rejection_sampler(n,p,kappa,radius)
  if(kappa == 0){
    Z = R*X
  } else if(kappa > 0) {
    Z0 = cos(R*sqrt(kappa))
    scale_factor = sqrt(1 - Z0**2)
    Z = cbind(Z0,scale_factor*X)
  } else if(kappa < 0) {
    Z0 = cosh(R*sqrt(abs(kappa)))
    scale_factor = sqrt(Z0**2 - 1)
    Z = cbind(Z0,scale_factor*X)
  }
  return(Z)
}

sim_projected_uniform_ball <- function(n,p,kappa,radius){
  # equal angle sampler
  X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = diag(rep(1,p)))
  X = X/sqrt(rowSums(X**2)) # direction X
  # setting sampled radius to be fixed
  if(kappa > 0 ){
    max.radius <- pi/sqrt(kappa)
    if(radius > max.radius){
      R = radial_rejection_sampler(n,p,0,max.radius)
    } else {
      R = radial_rejection_sampler(n,p,0,radius)
    }
  } else {
    R = radial_rejection_sampler(n,p,0,radius)
  }

  if(kappa == 0){
    Z = R*X
  } else if(kappa > 0) {
    Z0 = cos(R*sqrt(kappa))
    scale_factor = sqrt(1 - Z0**2)
    Z = cbind(Z0,scale_factor*X)
  } else if(kappa < 0) {
    Z0 = cosh(R*sqrt(abs(kappa)))
    scale_factor = sqrt(Z0**2 - 1)
    Z = cbind(Z0,scale_factor*X)
  }
  return(Z)
}

sim_projected_conic_distribution <- function(n,p,kappa,radius, flatness = 1){
  # equal angle sampler
  if(p > 2){
    sigma.flat = c(rep(flatness,2),rep(1, p - 2))
  } else {
    sigma.flat = rep(1,p)
  }

  X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = diag(sigma.flat))
  X = X/sqrt(rowSums(X**2)) # direction X
  # setting sampled radius to be fixed
  if(kappa > 0 ){
    max.radius <- pi/sqrt(kappa)
    if(radius > max.radius){
      R = runif(n,min = -max.radius/2, max = max.radius/2) + runif(n,min = -max.radius/2, max = max.radius/2)
      R = abs(R)
      #R = radial_rejection_sampler(n,p,0,max.radius)
    } else {
      #R = radial_rejection_sampler(n,p,0,radius)
      R = runif(n,min = -radius/2, max = radius/2) + runif(n,min = -radius/2, max = radius/2)
      R = abs(R)
    }
  } else {
    #R = radial_rejection_sampler(n,p,0,radius)
    R = runif(n,min = -radius/2, max = radius/2) + runif(n,min = -radius/2, max = radius/2)
    R = abs(R)
  }

  if(kappa == 0){
    Z = R*X
  } else if(kappa > 0) {
    Z0 = cos(R*sqrt(kappa))
    scale_factor = sqrt(1 - Z0**2)
    Z = cbind(Z0,scale_factor*X)
  } else if(kappa < 0) {
    Z0 = cosh(R*sqrt(abs(kappa)))
    scale_factor = sqrt(Z0**2 - 1)
    Z = cbind(Z0,scale_factor*X)
  }
  return(Z)
}

# write down how I sampled from each of these

#rejection sampler for radius
radial_rejection_sampler <- function(n,p,kappa,radius){
  if(kappa == 0){
    M = (radius)**p
    block_sample_size = round(5*M + 20)
    R.sample <- c()
    while(length(R.sample) < n){
      U <- runif(block_sample_size)
      Rtrial <- runif(block_sample_size,0,radius)
      accept <- U <= (1/M)*(Rtrial)**p
      R.sample <- c(R.sample, Rtrial[accept])
    }
    R.sample <- R.sample[1:n]
  } else if(kappa > 0) {
    if(sqrt(kappa)*radius > pi){
      stop("radius cannot exceed pi/sqrt(kappa)")
    }
    if(sqrt(kappa)*radius > pi/2){
      M = (sin(pi/2)/sqrt(kappa))**p
    } else {
      M = (sin(sqrt(kappa)*radius)/sqrt(kappa))**p
    }

    block_sample_size = round(5*M + 20)
    R.sample <- c()
    while(length(R.sample) < n){
      U <- runif(block_sample_size)
      Rtrial <- runif(block_sample_size,0,radius)
      accept <- U <= (1/M)*(sin(sqrt(kappa)*Rtrial)/sqrt(kappa))**p
      R.sample <- c(R.sample, Rtrial[accept])
    }
    R.sample <- R.sample[1:n]
  } else if(kappa < 0) {
    M = (sinh(sqrt(-kappa)*radius)/sqrt(-kappa))**p
    block_sample_size = round(5*M + 20)
    R.sample <- c()
    while(length(R.sample) < n){
      U <- runif(block_sample_size)
      Rtrial <- runif(block_sample_size,0,radius)
      accept <- U <= (1/M)*(sinh(sqrt(-kappa)*Rtrial)/sqrt(-kappa))**p
      R.sample <- c(R.sample, Rtrial[accept])
    }
    R.sample <- R.sample[1:n]
  }
  return(R.sample)

}


# pos_to_dist <- function(Z, kappa){
#   n = nrow(Z)
#   id.pairs = expand.grid(1:n,1:n)
#   id.pairs = id.pairs[id.pairs[,1] > id.pairs[,2],]
#   n.pairs = nrow(id.pairs)
#   D = matrix(data = 0, nrow = n, ncol = n)
#   if(kappa == 0){
#     dists = mclapply(1:n.pairs, function(x){
#       z = id.pairs[x,]
#       d = sqrt(sum((Z[z[[1]],] - Z[z[[2]],])**2)) # Assignment outside of the loop
#       #D[z[[1]],z[[2]]] <<- sqrt(sum((Z[z[[1]],] - Z[z[[2]],])**2)) # Assignment outside of the loop
#       return(d)
#     })
#   } else if (kappa < 0){
#     dists = mclapply(1:n.pairs, function(x){
#       z = id.pairs[x,]
#       # bilinear form
#       B.vec = Z[z[[1]],] * Z[z[[2]],]
#       B = 2*B.vec[1] - sum(B.vec)
#       d = (1/sqrt(abs(kappa)))*acosh(B)
#       #D[z[[1]],z[[2]]] <<- d # Assignment outside of the loop
#       return(d)
#     })
#   } else if (kappa > 0){
#     dists = mclapply(1:n.pairs, function(x){
#       z = id.pairs[x,]
#       # bilinear form
#       B.vec = Z[z[[1]],] * Z[z[[2]],]
#       B = sum(B.vec)
#       d = (1/sqrt(abs(kappa)))*acos(B)
#       return(d) # Assignment outside of the loop
#     })
#   }
#   for(i in 1:n.pairs){
#     D[id.pairs[i,1],id.pairs[i,2]] <- dists[[i]]
#   }
#   D = D + t(D)
#   return(D)
# }



# Alternative version parallelized
pos_to_dist_pair <- function(Zi, Zj, kappa){
  n = nrow(Zi)
  p = ncol(Zi)

  if(kappa == 0){
    #bilinear form
    Bmat = Zi %*% t(Zj)
    norms2_i = rowSums(Zi^2)
    norms2_j = rowSums(Zj^2)
    #D2 = log(exp(norms_i) %*% t(exp(norms_j))) - 2*Bmat
    D2 = outer(norms2_i,norms2_j,"+") - 2*Bmat
    #D2 = outer(Bmat,Bmat,"+")
    # rounding errors may occur
    D2[D2 < 0] = 0
    D = sqrt(D2)
  } else if (kappa < 0){

    sigmat = diag(c(1,rep(-1,(p - 1))))
    Bmat = Zi %*% sigmat %*% t(Zj)
    # rounding errors may occur
    Bmat[Bmat < 1] = 1
    D = (1/sqrt(abs(kappa)))*acosh(Bmat)
  } else if (kappa > 0){
    Bmat = Zi %*% t(Zj)
    #norms = sqrt(diag(Bmat))
    # ensures embedding on the unit sphere
    #Bmat = (norms^(-1)) %*% t((norms^(-1))) * Bmat
    # rounding errors may occur
    Bmat[Bmat > 1] = 1
    D = (1/sqrt(abs(kappa)))*acos(Bmat)
  }
  return(D)
}


pos_to_dist <- function(Z, kappa){
  n = nrow(Z)
  p = ncol(Z)
  D = matrix(data = 0, nrow = n, ncol = n)
  if(kappa == 0){
    #bilinear form
    Bmat = Z %*% t(Z)
    norms2 = rowSums(Z^2)
    #D2 = log(exp(norms) %*% t(exp(norms))) - 2*Bmat
    D2 = outer(norms2,norms2,"+") - 2*Bmat
    #D2 = outer(Bmat,Bmat,"+")
    # rounding errors may occur
    D2[D2 < 0] = 0
    D = sqrt(D2)
  } else if (kappa < 0){

    sigmat = diag(c(1,rep(-1,(p - 1))))
    Bmat = Z %*% sigmat %*% t(Z)
    # rounding errors may occur
    Bmat[Bmat < 1] = 1
    D = (1/sqrt(abs(kappa)))*acosh(Bmat)
  } else if (kappa > 0){
    Bmat = Z %*% t(Z)
    norms = sqrt(diag(Bmat))
    # ensures embedding on the unit sphere
    Bmat = (norms^(-1)) %*% t((norms^(-1))) * Bmat
    # rounding errors may occur
    Bmat[Bmat > 1] = 1
    D = (1/sqrt(abs(kappa)))*acos(Bmat)
  }
  return(D)
}



# simulate adjacency matrix
sim_ls_network <- function(rand.eff, D){
  # check that dimensions of rand.eff and D match, and check the positivity conditions
  nu.mat = exp(rand.eff) %*% t(exp(rand.eff))
  diag(nu.mat) =  0
  P = nu.mat*exp(-D)
  P.set = P[lower.tri(P)]

  A = matrix(data = 0, nrow = nrow(P), ncol = ncol(P))
  U = runif(n = length(P[lower.tri(P)]))
  A[lower.tri(A)] = ifelse(U <= P[lower.tri(P)], 1, 0)
  A = A + t(A)
  return(A)
}

sim_diffusion_network <- function(rand.eff, D, sig = 1){
  # check that dimensions of rand.eff and D match, and check the positivity conditions
  nu.mat = exp(rand.eff) %*% t(exp(rand.eff))
  diag(nu.mat) =  0
  P = nu.mat*exp(-D^2/sig^2)
  P.set = P[lower.tri(P)]

  A = matrix(data = 0, nrow = nrow(P), ncol = ncol(P))
  U = runif(n = length(P[lower.tri(P)]))
  A[lower.tri(A)] = ifelse(U <= P[lower.tri(P)], 1, 0)
  A = A + t(A)
  return(A)
}


# Z = Z.true
# kappa = kappa.true
# rand.eff = nu.vec

sim_ls_network_fast <- function(rand.eff,Z,kappa){
  n = length(rand.eff)
  connect.i <- c()
  connect.j <- c()
  p = ncol(Z)

  i = 1
  j = 2
  z.j <- Z[j,]
  z.i <- Z[i,]
  if(kappa == 0){
    #bilinear form
    #Bmat = z.j %*% t(z.i)
    #norms = diag(Bmat)
    norm.vec <- t(z.i) - z.j
    d.vec <- sqrt(sum(norm.vec**2))
  } else if (kappa < 0){
    sigmat = diag(c(1,rep(-1,(p - 1))))
    Bmat =   t(z.i*diag(sigmat)) %*% z.j
    # rounding errors may occur
    Bmat[Bmat < 1] = 1
    d.vec = (1/sqrt(abs(kappa)))*acosh(Bmat)
  } else if (kappa > 0){
    Bmat = t(z.j) %*% z.i
    #norms = sqrt(diag(Bmat))
    # ensures embedding on the unit sphere
    #Bmat = (norms^(-1)) %*% t((norms^(-1))) * Bmat
    # rounding errors may occur
    Bmat[Bmat > 1] = 1
    d.vec = (1/sqrt(abs(kappa)))*acos(Bmat)
  }
  nu.sum.vec <- rand.eff[i] + rand.eff[j]

  p.vec <- exp(nu.sum.vec - d.vec)
  n.sim <- length(p.vec)
  U <- runif(n.sim)
  a.vec <- 1*(U <= p.vec)
  if(a.vec == 1){
    connect.i <- c(connect.i,i)
    connect.j <- c(connect.j,j)
  }

  for(j in 3:(n)){
    if(j %% 1000 == 0){
      cat(paste("Mat Row:", j,"/", n), end = "\r")
    }
    i.seq <- seq(j - 1)
    z.j <- Z[j,]
    z.i <- Z[i.seq,]
    if(kappa == 0){
      #bilinear form
      #Bmat = z.j %*% t(z.i)
      #norms = diag(Bmat)
      norm.vec <- t(z.i) - z.j
      d.vec <- sqrt(colSums(norm.vec**2))
    } else if (kappa < 0){
      sigmat = diag(c(1,rep(-1,(p - 1))))
      Bmat = z.j %*% sigmat %*% t(z.i)
      # rounding errors may occur
      Bmat[Bmat < 1] = 1
      d.vec = (1/sqrt(abs(kappa)))*acosh(Bmat)
    } else if (kappa > 0){
      Bmat = z.j %*% t(z.i)
      #norms = sqrt(diag(Bmat))
      # ensures embedding on the unit sphere
      #Bmat = (norms^(-1)) %*% t((norms^(-1))) * Bmat
      # rounding errors may occur
      Bmat[Bmat > 1] = 1
      d.vec = (1/sqrt(abs(kappa)))*acos(Bmat)
    }
    d.vec <- as.numeric(d.vec)
    nu.sum.vec <- rand.eff[i.seq] + rand.eff[j]

    p.vec <- exp(nu.sum.vec - d.vec)
    n.sim <- length(p.vec)
    U <- runif(n.sim)
    a.vec <- 1*(U <= p.vec)
    i.seq.connected <- i.seq[which(a.vec == 1)]
    connect.i <- c(connect.i,i.seq.connected)
    connect.j <- c(connect.j, rep(j,length(i.seq.connected)))
  }
  A <- sparseMatrix(i = connect.i, j = connect.j,
                    x = rep(1,length(connect.i)), dims = c(n,n))
  A <- A + t(A)
  return(A)
}


# rand.eff <- nu.vec
# Z <- Z.true
# kappa <- kappa.true
sim_ls_network_fast_2 <- function(rand.eff,Z,kappa, max.n = 5000){
  n = length(rand.eff)
  connect.i <- c()
  connect.j <- c()
  p = ncol(Z)
  n.blocks <- ceiling(n/max.n)
  A <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(n,n))
  for(i in seq(n.blocks)){
    for(j in seq(n.blocks)){
      if(i < j){
        cat(paste0("block pair: (", i,",",j,")/","(", n.blocks,",",n.blocks,")"), end = "\r")
        i.set <- (max.n*(i - 1) + 1):(max.n*(i - 1) + max.n)
        j.set <- (max.n*(j - 1) + 1):(max.n*(j - 1) + max.n)
        i.set <- i.set[i.set <= n]
        j.set <- j.set[j.set <= n]

        Zi <- Z[i.set,]
        Zj <- Z[j.set,]
        D.sub <- pos_to_dist_pair(Zi, Zj, kappa)
        nui <- rand.eff[i.set]
        nuj <- rand.eff[j.set]

        nu.sub <- outer(nui,nuj,"+")
        p.sub <- exp(nu.sub - D.sub)
        U <- runif(length(p.sub))
        input = 1*(U <= p.sub)
        #A[i.set,j.set] = input
        A.sub <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(n,n))
        A.sub[i.set,j.set] = input
        A <- A + A.sub
      } else if(i == j){
        cat(paste0("block pair: (", i,",",j,")/","(", n.blocks,",",n.blocks,")"), end = "\r")
        i.set <- (max.n*(i - 1) + 1):(max.n*(i - 1) + max.n)
        j.set <- (max.n*(j - 1) + 1):(max.n*(j - 1) + max.n)
        i.set <- i.set[i.set <= n]
        j.set <- j.set[j.set <= n]

        Zi <- Z[i.set,]
        Zj <- Z[j.set,]
        D.sub <- pos_to_dist_pair(Zi, Zj, kappa)
        nui <- rand.eff[i.set]
        nuj <- rand.eff[j.set]

        nu.sub <- outer(nui,nuj,"+")
        p.sub <- exp(nu.sub - D.sub)
        U <- runif(length(p.sub))
        input <- 1*(U <= p.sub)

        diag(input) = 0
        input[lower.tri(input)] = 0

        A.sub <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(n,n))
        A.sub[i.set,j.set] = input
        A <- A + A.sub
      }
    }
  }
  # A[lower.tri(A)] = 0
  # diag(A) = 0
  A <- A + t(A)
  return(A)
}

# another option for making this faster,
# block versions of the whole thing.


# function for separating cliques into a non-overlapping set
# distinct_cliques <- function(clique.mat){
#   idx = 1:nrow(clique.mat)
#   idx.short = c()
#   nodes.counted = c()
#   for(i in idx){
#     if(any(clique.mat[i,] %in% nodes.counted)){
#
#     } else {
#       idx.short = c(idx.short, i)
#       nodes.counted = c(nodes.counted, clique.mat[i,])
#     }
#   }
#   distinct.cliques = clique.mat[idx.short,]
#   return(distinct.cliques)
# }
#
# # Currently pre-defined functions for the non-overlapping clique search
# # Could be made more slick in the future
#
# non_overlapping_clique_search <- function(A, clique.size){
#   l = clique.size
#   g = graph_from_adjacency_matrix(A, mode = "undirected")
#   clique.list = cliques(g, min=l, max=l)
#   clique.mat = matrix(unlist(clique.list), ncol = l, byrow = T)
#   distinct.cliques = distinct_cliques(clique.mat)
#   return(distinct.cliques)
# }
#
# # Is this going to be a problem, the fact we are choosing the second largest clique set.
# demaximalize_cliques <- function(distinct.cliques){
#   t <- ncol(distinct.cliques)
#   K <- nrow(distinct.cliques)
#   nm.cliques <- sapply(1:K, function(z){
#     sub.idx <- sample(1:t,(t - 1), replace = F)
#     out <- distinct.cliques[z,sub.idx]
#     return(out)
#   })
#   out <- t(nm.cliques)
#   return(out)
# }
#
# # estimate the fixed effects of a particular clique nodes.
# clique_fixed_effect_ard <- function(A,clique.idx){
#   Ysums = colSums(A[-clique.idx,clique.idx])
#   nus = c(log(Ysums/Ysums[1])) # estimated fixed effects with shift
#   return(nus)
# }
#
#
# # y must be a binary vector
# # vectors must be of the same length
# likelihood_vector <- function(y, fixed.effects) {
#   l = length(y)
#   s = sum(y)
#   ones.idx = (y == 1)
#   zeros.idx = (y == 0)
#   V1 = exp(fixed.effects)[ones.idx]
#   V0 = exp(fixed.effects)[zeros.idx]
#
#   if (length(V1) != 0){
#     c.coef = es_poly(V1)[s+1]
#   } else {
#     c.coef = 1
#   }
#
#   if(length(V0) == 0){
#     vec.tmp = 1
#     alt.series = 1
#   } else {
#     vec.tmp = es_poly(V0)
#     alt.series = (-1)*(-1)**(1:(length(vec.tmp)))
#   }
#   lik.vec = alt.series*c.coef*vec.tmp
#   if(length(lik.vec) != (l + 1)){
#     lik.vec = c(rep(0,1 + l -length(lik.vec)),lik.vec)
#   }
#   return(lik.vec)
# }
#
# # check if m and vec are the same length
# mgf_lik_score <- function(m,vec){
#   score = vec/(sum(m * vec))
#   return(score)
# }
#
# mgf_lik_hessian <- function(m,vec){
#   hess = -(vec %*% t(vec))/((sum(m * vec))**2)
#   return(hess)
# }
#
# mgf_dual <- function(nu,M,b,w){
#   # w a vector of weights corresponding to the observations
#   g = -sum(b * nu) + sum(w * log(t(M) %*% nu))
#   return(g)
# }
#
# mgf_dual_grad <- function(nu,M,b, w){
#   q = as.vector(w/(t(M) %*% nu))
#   dg = -b + M %*% q
#   return(dg)
# }
#
# mgf_dual_hess <- function(nu,M,b, w){
#   H <- matrix(data = NA, nrow = 2, ncol = 2)
#   for(i in 1:2){
#     for(j in 1:2){
#       q = as.vector(w/(t(M) %*% nu)^2)
#       H[i,j] = -sum(M[i,]*M[j,]*q)
#     }
#   }
#   return(H)
# }

# Function: compute empirical distribution function
# Input: N, (support size {0,...,N}); integer
#        x, sample of observations; vector of integers
# Output: Empirical distribution function for x
compute_edf <- function(x,N){
  p.hat <- c()
  for(y in 0:(N)){
    prop <- mean(x == y)
    p.hat <- c(p.hat, prop)
  }
  return(p.hat)
}

# mgf_alpha_estimation_dual <- function(A,clique.idx, fixed.effects, thresh = 10**(-6)){
#   l = length(clique.idx)
#   Y.samp = A[-clique.idx,clique.idx]
#   n.samp = nrow(Y.samp)
#   S.counts <- rowSums(Y.samp)
#   ell = ncol(Y.samp)
#   Sn <- compute_edf(S.counts,ell)
#
#   # compute the vectors V for the pattern vectors.
#   Y.patterns = expand.grid(rep(list(0:1),l))
#   S.key <- rowSums(Y.patterns)
#   lik.patterns = matrix(NA, nrow = nrow(Y.patterns), ncol = (ncol(Y.patterns) + 1))
#
#   lapply(1:nrow(Y.patterns), function(z){
#     pat = Y.patterns[z,]
#     lik.patterns[z,] <<- likelihood_vector(pat,fixed.effects)
#   })
#
#   V <- matrix(NA, nrow = (ell + 1), ncol = (ell + 1))
#
#   for(s in 0:ell){
#     if(s %in% c(0,ell)){
#       s.pattern = lik.patterns[S.key == s, ]
#     } else {
#       s.pattern = colSums(lik.patterns[S.key == s, ])
#     }
#
#     V[s + 1, ] = s.pattern
#   }
#   r = colSums(V)
#   # Closed form solution:
#   V.inv = solve(V)
#   q = t(V.inv) %*% c(1,rep(0,l))
#   z.opt = Sn/q
#   m.opt = V.inv %*% z.opt
#
#   # setting all moments to be 1, is this a reasonable initialization scheme,
#   # this should imply a proper cdf since this would imply a point mass at 0, it is indeed a feasible
#   # result, the only requirement is to convert it to a dual variable.
#   # or we could just find a dual variable.
#   return(as.vector(m.opt))
# }
#
# mgf_alpha_estimation <- function(A,clique.idx, fixed.effects, thresh = 10**(-6)){
#   l = length(clique.idx)
#   Y.samp = A[-clique.idx,clique.idx]
#   n.samp = nrow(Y.samp)
#
#   Y.patterns = expand.grid(rep(list(0:1),l))
#   count.patterns = rep(0,nrow(Y.patterns))
#   lapply(1:nrow(Y.patterns), function(z){
#     pat = Y.patterns[z,]
#
#     s = lapply(1:n.samp,function(x){
#       out = all(Y.samp[x,] == pat)
#     })
#     count.patterns[z] <<- sum(unlist(s))
#   })
#   n.eff <- sum(count.patterns) # effective sample size
#   # likelihood vector patterns corresponding to the outcome patterns
#   lik.patterns = matrix(NA, nrow = nrow(Y.patterns), ncol = (ncol(Y.patterns) + 1))
#
#   # use the loop for fast global assignment. TODO:  Maybe its faster to do something else.
#   tmp <- lapply(1:nrow(Y.patterns), function(z){
#     pat = Y.patterns[z,]
#     lik.patterns[z,] <<- likelihood_vector(pat,fixed.effects)
#   })
#   rm(tmp)
#
#   # This is just a particular choice of initialization that seemed to work reasonably
#   m.last = (0.9)^(1:l) # keeping all moments to be 1,
#   m.diff = Inf
#   while(m.diff > thresh){
#     m.last.long = c(1,m.last)
#
#     # computing entirety of score function.
#     score.array = sapply(1:nrow(Y.patterns), function(z){
#       vec = lik.patterns[z,]
#
#       weight = count.patterns[z]
#       out = weight*mgf_lik_score(m.last.long,vec)
#       return(out)
#     })
#     score.vec = rowSums(score.array)/n.eff
#     # computing entirety of hessian function.
#     hess.array = sapply(1:nrow(Y.patterns), function(z){
#       vec = lik.patterns[z,]
#       weight = count.patterns[z]
#       out = weight*mgf_lik_hessian(m.last.long,vec)
#       return(out)
#     })
#     hess.vec = rowSums(hess.array)/n.eff
#     hess = matrix(hess.vec, nrow = (l + 1), ncol = (l + 1), byrow = T)
#     score.sub = score.vec[-1]
#     hess.sub = hess[-1,-1]
#     inv.hess = solve(hess.sub)
#     # replace inv.hess with a constant for gradient descent
#     #
#     m.next = m.last - inv.hess %*% score.sub
#     m.diff = sqrt(sum((m.next - m.last)**2))
#     m.last = m.next
#   }
#   print("Return exponential moments: ")
#   return(as.vector(m.next))
# }



# Add c.hat to the estimate for the fixed effects to
scale_estimate_base <- function(moments,lead.poly = 2){
  l = length(moments)
  # effective moments
  eff.moments <- moments[moments > 0]
  l.eff <- length(eff.moments)
  c.hat = log(moments[l.eff]/moments[l.eff - 1]) + lead.poly*(log(1 - 1/l.eff))
  if(l.eff < l){
    warning(paste0("Highest Moment estimate is 0, used the (", l.eff - 1, ") moment instead"))
  }
  return(c.hat)
}



D_estimate <- function(A,clique.legend,fixed.effects){
  # rows of clique legend represent distinct cliques
  K = nrow(clique.legend)
  l = ncol(clique.legend)

  idx.pairs = expand.grid(1:K,1:K)
  idx.pairs = idx.pairs[idx.pairs[,1] < idx.pairs[,2],]
  D.hat = matrix(0, nrow = K, ncol = K)
  D.var.hat = matrix(0, nrow = K, ncol = K)
  n.idx = nrow(idx.pairs)
  for(i in 1:n.idx){
    x = idx.pairs[i,1]
    y = idx.pairs[i,2]
    x.idx = clique.legend[x,]
    y.idx = clique.legend[y,]
    p.xy.hat =  mean(A[x.idx, y.idx])
    gamma.x = mean(exp(fixed.effects[x,]))
    gamma.y = mean(exp(fixed.effects[y,]))
    d.xy.hat = log(gamma.x) + log(gamma.y) - log(p.xy.hat)
    D.hat[x,y] = d.xy.hat
  }
  D.hat <- D.hat + t(D.hat)

  out.list <- list(estimates = D.hat, variances = D.var.hat)
  return(out.list)
}


# This is a simple first take at such a function. This can be made more slick in the future
# Certainly we'd want to get rid of these for loops
# D_estimate will also estimate the variance of the d matrix
# TODO: Remove Variance Estimate
D_estimate <- function(A,clique.legend,fixed.effects){
  # rows of clique legend represent distinct cliques
  K = nrow(clique.legend)
  l = ncol(clique.legend)

  idx.pairs = expand.grid(1:K,1:K)
  idx.pairs = idx.pairs[idx.pairs[,1] < idx.pairs[,2],]
  D.hat = matrix(0, nrow = K, ncol = K)
  D.var.hat = matrix(0, nrow = K, ncol = K)
  n.idx = nrow(idx.pairs)
  d.list <- lapply(1:n.idx, function(i){
    x = idx.pairs[i,1]
    y = idx.pairs[i,2]
    x.idx = clique.legend[x,]
    y.idx = clique.legend[y,]
    p.xy.hat =  mean(A[x.idx, y.idx])
    gamma.x = mean(exp(fixed.effects[x,]))
    gamma.y = mean(exp(fixed.effects[y,]))
    d.xy.hat = log(gamma.x) + log(gamma.y) - log(p.xy.hat)
    return(d.xy.hat)})

  d.vec = unlist(d.list)
  D.hat[upper.tri(D.hat)] = d.vec
  D.hat <- D.hat + t(D.hat)


  out.list <- list(estimates = D.hat, variances = D.var.hat)
  return(out.list)
}






# best frechet mean of two points of a distance matrix
frechet_mean <- function(D,x,y){
  idx <- 1:nrow(D)
  idx <- idx[-c(x,y)]
  dxm <- D[x,idx]
  dym <- D[y,idx]
  obj = (dxm**2 + dym**2)
  m.hat = idx[which.min(obj)]
  return(m.hat)
}

# midpoint objective function
midpoint_objective <- function(dym,dzm,dyz){
  out = (dym**2 + dzm**2)/(dyz**2)
  return(out)
}

# midpoint objective function
midpoint_objective2 <- function(dxm,dym,dxy){
  pxm = min(exp(-dxm),1)
  pym = min(exp(-dym),1)
  pxy = min(exp(-dxy),1)

  out = (dxm**2 + dym**2)/(dxy**2)*(log(pxm*(1 - pxm)*pym*(1 - pym)*pxm*(1 - pxm)))
  return(out)
}



# how many configurations of the top midpoints should we look for
# Also could look for midpoints within a certain bias
# optimal midpoint objective may depend whether we include a restriction
# that D.hat is in fact a metric
# optimal_midpoint_search <- function(D,top.k = 1, d.yz.min = 1.0,d.yz.max = 2.5){
#   K <- nrow(D)
#   idx.set <- expand.grid(1:K,1:K,1:K)
#   idx.set <- idx.set[idx.set[,1] < idx.set[,2] & idx.set[,2] != idx.set[,3] & idx.set[,1] != idx.set[,3],]
#   n.idx <- nrow(idx.set)
#   obj.val <- sapply(1:n.idx, function(i){
#     y = idx.set[i,1]
#     z = idx.set[i,2]
#     m = idx.set[i,3]
#     dyz <- D[y,z]
#     dym <- D[y,m]
#     dzm <- D[z,m]
#
#     #filter out infinite distances:
#     dist.filter <- ifelse(dzm*dym*dyz == Inf, NaN, 1)
#
#     obj.tmp = midpoint_objective(dzm,dym,dyz)*dist.filter
#     # ideal midpoint will have a d\istance of 1/2
#     # this will also balance the distances
#     # removing sets with large dyz, i.e. a single point mass.
#     out <- abs(obj.tmp - 1/2) + abs(dym - dzm)/dyz + 10000*(dyz > d.yz.max) +  10000*(dyz < d.yz.min)
#     return(out)
#   })
#   # objective value must be > 1/2 by def
#   #
#
#
#   sort.obj <- obj.val[order(obj.val)]
#   ranked.idx <- idx.set[order(obj.val),]
#   #ranked.idx <- ranked.idx[sort.obj, ]
#   selected.set <- c()
#   best.pairs <- matrix(NA, nrow = 0, ncol = 6)
#
#   i = 1
#   # stops the loop if we do not have a full top k which don't overlap
#   while(nrow(best.pairs) < top.k & i <= nrow(ranked.idx)){
#     if(!(any(as.vector(ranked.idx[i,]) %in% selected.set))){
#       D.tmp <- D[as.numeric(ranked.idx[i,]),as.numeric(ranked.idx[i,])]
#       best.pairs <- rbind(best.pairs, c(as.numeric(ranked.idx[i,]), D.tmp[1,3],D.tmp[2,3], D.tmp[1,2]))
#       selected.set <- c(selected.set, c(ranked.idx[i,1],ranked.idx[i,2],ranked.idx[i,3]))
#     }
#
#     i = i+1
#   }
#   colnames(best.pairs) <- c("y","z","m","dym", "dzm", "dyz")
#   rownames(best.pairs) <- NULL
#   return(best.pairs)
# }
#
#
#
# # TODO: replace the existing one with this
# optimal_midpoint_search <- function(D,top.k = 1, d.yz.min = 1.5,d.yz.max = 2.5){
#   K <- nrow(D)
#   idx.set.yz <- expand.grid(1:K,1:K)
#   idx.set.yz <- idx.set.yz[idx.set.yz[,1] < idx.set.yz[,2],]
#   dyz.set <- c()
#   for(k in seq(nrow(idx.set.yz))){
#     dyz.set[k] <- D[as.numeric(idx.set.yz[k,1]), as.numeric(idx.set.yz[k,2])]
#   }
#
#   idx.set.yz <- idx.set.yz[dyz.set < d.yz.max & d.yz.min < dyz.set, ]
#
#   idx.set <- matrix(NA, nrow = 0, ncol = 3)
#   for(k in seq(nrow(idx.set.yz))){
#     idx.block <- expand.grid(idx.set.yz[k,1],idx.set.yz[k,2], 1:K)
#     idx.set <- rbind(idx.set,idx.block)
#   }
#   idx.set <- idx.set[idx.set[,2] != idx.set[,3] & idx.set[,1] != idx.set[,3],]
#
#   n.idx <- nrow(idx.set)
#   obj.val <- sapply(1:n.idx, function(i){
#     y = idx.set[i,1]
#     z = idx.set[i,2]
#     m = idx.set[i,3]
#     dyz <- D[y,z]
#     dym <- D[y,m]
#     dzm <- D[z,m]
#
#     #filter out infinite distances:
#     dist.filter <- ifelse(dzm*dym*dyz == Inf, NaN, 1)
#
#     obj.tmp = midpoint_objective(dzm,dym,dyz)*dist.filter
#     # ideal midpoint will have a distance of 1/2
#     # this will also balance the distances
#     # removing sets with large dyz, i.e. a single point mass.
#     out <- obj.tmp #+ abs(dym - dzm)/dyz #abs(obj.tmp - 1/2) + abs(dym - dzm)/dyz + 10000*(dyz > d.yz.max) +  10000*(dyz < d.yz.min)
#     return(out)
#   })
#   # objective value must be > 1/2 by def
#   #
#
#
#   sort.obj <- obj.val[order(obj.val)]
#   ranked.idx <- idx.set[order(obj.val),]
#   #ranked.idx <- ranked.idx[sort.obj, ]
#   selected.set <- c()
#   best.pairs <- matrix(NA, nrow = 0, ncol = 6)
#
#   i = 1
#   # stops the loop if we do not have a full top k which don't overlap
#   while(nrow(best.pairs) < top.k & i <= nrow(ranked.idx)){
#     if(!(any(as.vector(ranked.idx[i,]) %in% selected.set))){
#       D.tmp <- D[as.numeric(ranked.idx[i,]),as.numeric(ranked.idx[i,])]
#       best.pairs <- rbind(best.pairs, c(as.numeric(ranked.idx[i,]), D.tmp[1,3],D.tmp[2,3], D.tmp[1,2]))
#       selected.set <- c(selected.set, c(ranked.idx[i,1],ranked.idx[i,2],ranked.idx[i,3]))
#     }
#
#     i = i+1
#   }
#   colnames(best.pairs) <- c("y","z","m","dym", "dzm", "dyz")
#   rownames(best.pairs) <- NULL
#   return(best.pairs)
# }




# how many configurations of the top midpoints should we look for
# Also could look for midpoints within a certain bias
# optimal_midpoint_search_2 <- function(D,top.k = 1){
#   K <- nrow(D)
#   idx.set <- expand.grid(1:K,1:K,1:K)
#   idx.set <- idx.set[idx.set[,1] < idx.set[,2] & idx.set[,2] < idx.set[,3],]
#   n.idx <- nrow(idx.set)
#   obj.val <- sapply(1:n.idx, function(i){
#     y = idx.set[i,1]
#     z = idx.set[i,2]
#     m = idx.set[i,3]
#     dyz <- D[y,z]
#     dym <- D[y,m]
#     dzm <- D[z,m]
#
#     #filter out infinite distances:
#     dist.filter <- ifelse(dzm*dym*dyz == Inf, NaN, 1)
#
#     out = midpoint_objective2(dzm,dym,dyz)*dist.filter
#     return(out)
#   })
#   # objective value must be > 1/2 by def
#   #
#
#
#   sort.obj <- obj.val[order(obj.val)]
#   ranked.idx <- idx.set[order(obj.val),]
#   ranked.idx <- ranked.idx[sort.obj >= 1/2, ]
#   selected.set <- c()
#   best.pairs <- matrix(NA, nrow = 0, ncol = 3)
#
#   i = 1
#   # stops the loop if we do not have a full top k which don't overlap
#   while(nrow(best.pairs) < top.k & i <= nrow(ranked.idx)){
#     if(!(any(as.vector(ranked.idx[i,]) %in% selected.set))){
#       best.pairs <- rbind(best.pairs, as.vector(ranked.idx[i,]))
#       selected.set <- c(selected.set, c(ranked.idx[i,1],ranked.idx[i,2],ranked.idx[i,3]))
#     }
#
#     i = i+1
#   }
#   colnames(best.pairs) <- c("y","z","m")
#   rownames(best.pairs) <- NULL
#   return(best.pairs)
# }



# searches for the optimal variance x
# optimal_x_search <- function(D,y,z,m){
#   K <- ncol(D)
#   idx <- 1:K
#   idx <- idx[-c(y,z,m)]
#   obj <- abs(D[idx,y]^2 - D[idx,z]^2)/(D[idx,y]^2 + D[idx,z]^2)
#   # discarding points which have no connections
#   obj.thresh <- ifelse(D[idx,y]*D[idx,z]*D[idx,m]  == Inf, Inf, 1)
#   obj <- obj*obj.thresh
#
#   x.idx <- idx[which.min(obj)]
#
#   return(x.idx)
# }


# bias compute
# estimate of the bias of our estimates
# performs a test of whether constant curvature
# given a set of indexing points, this essentially equates to an anova test with unequal variances
# Welch's Anova
# intra_network_constant_curvature_test <- function(kappa.hat,sd.kappa.hat){
#   R <- length(kappa.hat)
#   weights <- 1/sd.kappa.hat
#   kappa.pool <- sum(weights*kappa.hat)/(sum(weights))
#   chi.square.stat <- sum(((kappa.hat - kappa.pool)/sd.kappa.hat)**2)
#   p.value <- 1 - pchisq(chi.square.stat, df = R)
#   return(p.value)
# }


# input, a 3D array of multi-plex network
# given a set of indexing points, this essentially equates to an anova test with unequal variances
# Welch's Anova
# inter_network_constant_curvature <- function(kappa.hat,sd.kappa.hat){
#   R <- length(kappa.hat)
#   weights <- 1/sd.kappa.hat
#   kappa.pool <- sum(weights*kappa.hat)/(sum(weights))
#   chi.square.stat <- sum(((kappa.hat - kappa.pool)/sd.kappa.hat)**2)
#   p.value <- 1 - pchisq(chi.square.stat, df = R)
#   return(p.value)
# }

# input, a 3D array of ordered time series
# lambda is a tuning parameter
# will be a simple linear change point detection problem Harchaoui et al 2010
# sequential_network_change_point <- function(kappa.hat.seq,sd.kappa.hat.seq, lambda){
#   #casting this as a lasso problem
#   T = length(kappa.hat.seq)
#   X.des <- matrix(1,nrow = T - 1, ncol = T - 1)
#   X.des[upper.tri(X.des)] = 0
#   X.des <- cbind(rep(0,T - 1), X.des)
#
#   # Lasso with known unequal variances
#   y = kappa.hat.seq
#
# }


#
# nuissance_regression_mgf <- function(y, s.vec,discard.0 = F, thresh = 10**(-3), Tmax = 10000){
#   c.prev = 0
#   d.prev = 0
#   b1.prev = 1
#   b2.prev = -3
#   theta.prev <- c(c.prev,d.prev,b1.prev,b2.prev)
#   diff = Inf
#   iter = 0
#   while(iter <= Tmax & diff > thresh){
#     c.prev = theta.prev[1]
#     d.prev = theta.prev[2]
#     b1.prev = theta.prev[3]
#     b2.prev = theta.prev[4]
#
#     score = nuissance_grad(s.vec,y,c.prev,d.prev,b1.prev,b2.prev)
#     H = nuissance_hess(s.vec,y,c.prev,d.prev,b1.prev,b2.prev)
#     theta.next = theta.prev + solve(H)%*%score
#     iter = iter + 1
#     diff = min(abs(theta.next - theta.prev))
#     theta.prev = theta.next
#   }
#   #print(iter)
#   c = theta.prev[1]
#   d = theta.prev[2]
#   b1 = theta.prev[3]
#   b2 = theta.prev[4]
#   #c = -.1
#   #d = 2
#   #b1 = 2
#   #b2 = -3
#   pred = nuissance_model(s.vec,c,d,b1,b2)
#   #plot(s.vec,y)
#   #lines(s.vec, pred)
#   if(iter > Tmax){
#     print("max iterations")
#   }
#   out.list <- list("c.hat" = theta.prev[1], "predictions" = pred,
#                    "d" = d, "b1" = b1, "b2" = b2)
#   #print(out.list)
#   return(out.list)
# }
#
# nuissance_model <- function(s.vec,c,d,b1,b2){
#   q = c*s.vec + d + b2*log(s.vec + b1)
#   return(q)
# }
#
# nuissance_q <- function(s.vec,y,c,d,b1,b2){
#   q = y - c*s.vec - d - b2*log(s.vec + b1)
#   return(q)
# }
#
# nuissance_grad <- function(s.vec,y,c,d,b1,b2){
#   v1 = 2*s.vec
#   v2 = rep(2,length(s.vec))
#   v3 = b2/(b1 + s.vec)
#   v4 = 2*log(b1 + s.vec)
#   q = nuissance_q(s.vec,y,c,d,b1,b2)
#   grad = c(sum(v1*q), sum(v2*q),
#            sum(v3*q), sum(v4*q))
#   return(grad)
# }
#
#
# nuissance_hess <- function(s.vec,y,c,d,b1,b2){
#   h11 <- sum(2*s.vec^2)
#   h12 <- sum(2*s.vec)
#   h21 <- h12
#   h13 <- sum(2*b2*s.vec/(b1 + s.vec))
#   h31 <- h13
#   h14 <- sum(2*s.vec*log(b1 + s.vec))
#   h41 <- h14
#   h22 <- 2*length(s.vec)
#   h23 <- 2*sum(2*b2/(b1 + s.vec))
#   h32 <- h23
#   h24 <- 2*sum(log(b1 + s.vec))
#   h42 <- h24
#   h33 <- 2*sum(b2^2/(b1 + s.vec)^2)
#   h34 <- 2*sum(b2*nuissance_q(s.vec,y,c,d,b1,b2)/(b1 + s.vec)^2)
#   h43 <- h34
#   h44 <- 2*sum(log(b1 + s.vec)^2)
#   H.vec <- c(h11,h12,h13,h14,
#              h21,h22,h23,h24,
#              h31,h32,h33,h34,
#              h41,h42,h43,h44)
#   H <- matrix(H.vec, nrow = 4, ncol = 4)
#   return(H)
# }


colSDs <- function(df, ...){
  diff.df <- t(df) - colMeans(df, ...)
  diff.df <- t(diff.df)
  var.vec <- colMeans((diff.df)^2, ...)
  sd.vec <- sqrt(var.vec)
  return(sd.vec)
}




# a filtered median estimate of kappa.
# filter_median <- function(D, y,z,m, c1 = 1/2,c2 = 2,c3 = 0.75, min.K = 3){
#   K = nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(y,z,m)]
#
#   #indicator of whether to include an x
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
#       return(F)
#     } else {
#       i1 = (c1)*D[y,z] <= D[x,y]
#       i2 = (c1)*D[y,z] <= D[x,z]
#       i3 = D[x,z] <= (c2)*D[y,z]
#       i4 = D[x,y] <= (c2)*D[y,z]
#       # filtering region change
#       #i5 = abs(D[x,y]^2 - D[x,z]^2) <= c3*(D[x,y] + D[x,z] - D[y,z])
#       i5 = abs(D[x,y] - D[x,z]) <= 2*c3*(D[y,z])
#       return(ifelse(i1*i2*i3*i4*i5 == 1, T, F))
#     }
#
#   })
#
#   x.filtered <- x.set[x.id]
#   if(length(x.filtered) < min.K){
#     kappa.med <- NA
#   } else {
#     kappa.set <- sapply(x.filtered, function(x){
#       dxy = D[x,y]
#       dxz = D[x,z]
#       dyz = D[y,z]
#       dxm = D[x,m]
#       d.vec = c(dxy,dxz,dyz,dxm)
#       if(any(d.vec == Inf)){
#         kappa.hat.x <- NA
#       } else{
#         kappa.hat.x <- estimate_kappa(dxy,dxz,dyz,dxm)
#       }
#       return(kappa.hat.x)
#     })
#     kappa.med = median(kappa.set, na.rm = T)
#   }
#
#
#
#   return(kappa.med)
# }


# D_bootstrap <- function(A,clique.legend,nus.hat.mat){
#   # rows of clique legend represent distinct cliques
#   K = nrow(clique.legend)
#   ell = ncol(clique.legend)
#
#   idx.pairs = expand.grid(1:K,1:K)
#   idx.pairs = idx.pairs[idx.pairs[,1] > idx.pairs[,2],]
#   D.boot = matrix(0, nrow = K, ncol = K)
#   n.idx = nrow(idx.pairs)
#   for(i in 1:n.idx){
#     x = idx.pairs[i,1]
#     y = idx.pairs[i,2]
#     x.idx = clique.legend[x,]
#     y.idx = clique.legend[y,]
#     A.boot = sample(A[x.idx, y.idx], ell**2, replace = T)
#     p.xy.boot =  mean(A.boot)
#     gamma.x = mean(exp(nus.hat.mat[x,]))
#     gamma.y = mean(exp(nus.hat.mat[y,]))
#     # note that here we are not accounting for the randomness in
#     # gamma.x and gamma.y
#     # this is assumed to be small as these quantities are learned from
#     # the remainder of the network
#     d.xy.boot = log(gamma.x) + log(gamma.y) - log(p.xy.boot)
#     D.boot[x,y] = d.xy.boot
#   }
#
#
#   D.boot <- D.boot + t(D.boot)
#
#   return(D.boot)
# # }
#
#
# # vectorized version (TODO: Delete the old versions)
# D_bootstrap <- function(A,clique.legend,nus.hat.mat){
#   # rows of clique legend represent distinct cliques
#   K = nrow(clique.legend)
#   ell = ncol(clique.legend)
#
#   idx.pairs = expand.grid(1:K,1:K)
#   idx.pairs = idx.pairs[idx.pairs[,1] > idx.pairs[,2],]
#   D.boot = matrix(0, nrow = K, ncol = K)
#   n.idx = nrow(idx.pairs)
#   d.list <- lapply(1:n.idx, function(i){
#     x = idx.pairs[i,1]
#     y = idx.pairs[i,2]
#     x.idx = clique.legend[x,]
#     y.idx = clique.legend[y,]
#     A.boot = sample(A[x.idx, y.idx], ell**2, replace = T)
#     p.xy.boot =  mean(A.boot)
#     gamma.x = mean(exp(nus.hat.mat[x,]))
#     gamma.y = mean(exp(nus.hat.mat[y,]))
#     # note that here we are not accounting for the randomness in
#     # gamma.x and gamma.y
#     # this is assumed to be small as these quantities are learned from
#     # the remainder of the network
#     d.xy.boot = log(gamma.x) + log(gamma.y) - log(p.xy.boot)
#     return(d.xy.boot)
#   })
#   d.vec = unlist(d.list)
#   D.boot[lower.tri(D.boot)] = d.vec
#   D.boot <- D.boot + t(D.boot)
#
#   return(D.boot)
# }


#
#
# bootstrap_bias_se <- function(A,clique.legend,nus.hat.mat,
#                               B = 1000, c1 = 1/2, c2 = 2, c3 = 1.5){
#   K <- nrow(clique.legend)
#   l <- ncol(clique.legend)
#
#   D.est = D_estimate(A,clique.legend,nus.hat.mat)
#   D.hat = D.est$estimates
#
#   y.opt = 1
#   z.opt = 2
#   m.opt = 3
#
#   kappa.med <- filter_median(D.hat, y.opt,
#                              z.opt,m.opt,
#                              c1 = c1,c2 = c2,c3 = c3)
#
#   se.sum <- 0
#   bias.sum <- 0
#   kappa.med.boot.list <- lapply(1:B, function(z){
#     D.boot <- D_bootstrap(A,clique.legend,nus.hat.mat)
#     kappa.med.boot <- filter_median(D.boot, y.opt,
#                                     z.opt,m.opt,
#                                     c1 = c1,c2 = c2,c3 = c3)
#     return(kappa.med.boot)
#   })
#   kappa.med.boot.vec = unlist(kappa.med.boot.list)
#
#   se.est = sqrt(mean((kappa.med.boot.vec - kappa.med)^2, na.rm = T))
#   bias.est = (mean(kappa.med.boot.vec, na.rm = T)) - kappa.med
#
#   out.list = list("bias" = bias.est, "se" = se.est)
#   return(out.list)
# }

#
# D_full_from_cliques <- function(D,clique.legend){
#   K = nrow(D)
#   ell = ncol(clique.legend)
#   idx.pairs = expand.grid(1:K*ell,1:K*ell)
#   idx.pairs = idx.pairs[idx.pairs[,1] > idx.pairs[,2],]
#   D.full = matrix(0, nrow = K*ell, ncol = K*ell)
#   d.list <- lapply(1:K, function(i){
#     mat.col <- rep(D[i,], each = ell)
#     mat.block <- rep(mat.col, times = ell)
#     return(mat.block)
#   })
#   d.vec = unlist(d.list)
#   D.full <- matrix(d.vec, nrow = K*ell, ncol = K*ell, byrow = F)
#   return(D.full)
# }
#
# parametric_bootstrap <- function(A,clique.legend,nus.hat.mat,
#                                  B = 1000, c1 = 1/2, c2 = 2, c3 = 1.5){
#   K <- nrow(clique.legend)
#   l <- ncol(clique.legend)
#
#   D.est = D_estimate(A,clique.legend,nus.hat.mat)
#   D.hat = D.est$estimates
#
#   D.hat.full <- D_full_from_cliques(D.hat,clique.legend)
#
#   y.opt = 1
#   z.opt = 2
#   m.opt = 3
#
#   kappa.med <- filter_median(D.hat, y.opt,
#                              z.opt,m.opt,
#                              c1 = c1,c2 = c2,c3 = c3)
#   # transposing keeps the correct order
#   rand.eff <- as.vector(t(nus.hat.mat))
#   se.sum <- 0
#   bias.sum <- 0
#   # kappa.med.boot.list <- lapply(1:B, function(z){
#   #   # sample from fitted model
#   #   A.boot <- sim_ls_network(rand.eff,D.hat.full)
#   #
#   #   D.est.boot <- D_estimate(A.boot,
#   #                            clique.legend,
#   #                            nus.hat.mat)
#   #
#   #   D.boot = D.est.boot$estimates
#   #
#   #   # to filter or not to filter
#   #   # that is the question, (i.e. filter before, or
#   #   # bootstrap the whole thing)
#   #   #kappa.med.boot <- filter_median(D.boot, y.opt,z.opt,m.opt,
#   #   #                                c1 = c1,c2 = c2, c3 = c3)
#   #   kappa.med.boot <- filter_median(D.boot, y.opt,
#   #                                   z.opt,m.opt,
#   #                                   c1 = c1,c2 = c2,
#   #                                   c3 = c3)
#   #   return(kappa.med.boot)
#   # })
#   #kappa.med.boot.vec = unlist(kappa.med.boot.list)
#   kappa.med.boot.vec <- c()
#   for(b in 1:B){
#       # sample from fitted model
#       A.boot <- sim_ls_network(rand.eff,D.hat.full)
#       D.est.boot <- D_estimate(A.boot,
#                                clique.legend,
#                                nus.hat.mat)
#       D.boot = D.est.boot$estimates
#       kappa.med.boot <- filter_median(D.boot, y.opt,
#                                       z.opt,m.opt,
#                                       c1 = c1,c2 = c2,
#                                       c3 = c3)
#       #cat("Bootstrap ", b, "/", B, end = "\r")
#       kappa.med.boot.vec[b] = kappa.med.boot
#   }
#   se.est = sqrt(mean((kappa.med.boot.vec - kappa.med)^2, na.rm = T))
#   sd.est = sqrt(mean((kappa.med.boot.vec - mean(kappa.med.boot.vec, na.rm = T))^2, na.rm = T))
#   bias.est = (mean(kappa.med.boot.vec, na.rm = T)) - kappa.med
#
#   out.list = list("bias" = bias.est, "se" = se.est, "sd" = sd.est)
#   return(out.list)
# }


# filter_indices <- function(D,y,z,m, c1 = 1/2,c2 = 2,c3 = 0.5){
#   K = nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(y,z,m)]
#
#   #indicator of whether to include an x
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
#       return(F)
#     } else {
#       i1 = (c1)*D[y,z] <= D[x,y]
#       i2 = (c1)*D[y,z] <= D[x,z]
#       i3 = D[x,z] <= (c2)*D[y,z]
#       i4 = D[x,y] <= (c2)*D[y,z]
#       #i5 = abs(D[x,y]^2 - D[x,z]^2) <= c3*(D[x,y] + D[x,z] - D[y,z])
#       i5 = abs(D[x,y] - D[x,z]) <= 2*c3*(D[y,z])
#       #i5 = T
#       return(ifelse(i1*i2*i3*i4*i5 == 1, T, F))
#     }
#
#   })
#   x.filtered <- x.set[x.id]
#   return(x.filtered)
# }



# Intuition:
# triangle inequalities must hold up to a scaling factor
filter_indices_2 <- function(D,y,z,m, tri.const = sqrt(sqrt(2))){
  K = nrow(D)
  x.set <- 1:K
  x.set <- x.set[-c(y,z,m)]

  #indicator of whether to include an x
  x.id <- sapply(x.set, function(x){
    # ensuring that the distance is not infinite to each of x.y.m
    if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
      return(F)
    } else {
      i1 = D[y,z] + D[x,z] >= tri.const*D[x,y]
      i2 = D[y,z] + D[x,y] >= tri.const*D[x,z]
      i3 = D[x,y] + D[x,z] >= tri.const*D[y,z]

      return(ifelse(i1*i2*i3 == 1, T, F))
    }

  })
  x.filtered <- x.set[x.id]
  return(x.filtered)
}

filtered_curvature_median <- function(estimates, tri.const = sqrt(sqrt(2))){
  D <- estimates$D
  mid.search <- estimates$midpoints
  y.opt = mid.search[1,1]
  z.opt = mid.search[1,2]
  m.opt = mid.search[1,3]
  x.set <- filter_indices_2(D, y.opt,
                            z.opt,m.opt,
                            tri.const = tri.const)
  kappa.est.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
  kappa.est.set <- kappa.est.set[!is.na(kappa.est.set)]
  return(median(kappa.est.set))
}

constant_curvature_test <- function(estimates, num.midpoints = 3, tri.const = sqrt(sqrt(2))){
  D <- estimates$D
  mid.search <- estimates$midpoints

  kappa.sets <- list()
  index <- c()
  for(k in seq(num.midpoints)){
    y.opt = mid.search[k,1]
    z.opt = mid.search[k,2]
    m.opt = mid.search[k,3]
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set <- filter_indices_2(D, y.opt,
                              z.opt,m.opt,
                              tri.const = tri.const)
    kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
    kappa.set <- kappa.set[!is.na(kappa.set)]
    kappa.sets[[k]] = kappa.set
    index <- c(index, rep(k, length(kappa.set)))
  }

  kappa.vec <- unlist(kappa.sets)
  if(length(unique(index)) > 1 ){
    test.dat <- data.frame("loc" = index, "est" = kappa.vec)
    test <- kruskal.test(est ~ loc, data = test.dat)
    out.list <- list("p.value" =  test$p.value, "estimates" = test.dat)
  } else {
    out.list <- list("p.value" =  NULL, "estimates" = NULL)
  }

  return(out.list)
}

normalized_constant_curvature_test <- function(estimates, num.midpoints = 3, tri.const = sqrt(sqrt(2)), curve.scale = 10){
  D <- estimates$D
  mid.search <- estimates$midpoints

  trim.kappa.sets <- list()
  kappa.sets <- list()
  index <- c()
  for(k in seq(num.midpoints)){
    y.opt = mid.search[k,1]
    z.opt = mid.search[k,2]
    m.opt = mid.search[k,3]
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set <- filter_indices_2(D, y.opt,
                              z.opt,m.opt,
                              tri.const = tri.const)
    if(length(x.set) > 1){
      kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
      kappa.set <- kappa.set[!is.na(kappa.set)]
      y = scale_curvature(kappa.set, curve.scale)
      # transformation scale
      y = (y - median(y))/mad(y) + median(y)
      trim.kappa.sets[[k]] = y
      kappa.sets[[k]] = kappa.set
      index <- c(index, rep(k, length(kappa.set)))
    }
  }

  trim.kappa.vec <- unlist(trim.kappa.sets)
  kappa.vec <- unlist(kappa.sets)
  if(length(unique(index)) > 1 ){
    trim.test.dat <- data.frame("loc" = index, "est" = trim.kappa.vec)
    est.dat <- data.frame("loc" = index, "est" = kappa.vec)
    norm.test <- kruskal.test(est ~ loc, data = trim.test.dat)
    test <-  kruskal.test(est ~ loc, data = est.dat)
    out.list <- list("p.value" =  test$p.value,"norm.p.value" =norm.test$p.value, "estimates" = est.dat, "transformed_estimates"=trim.test.dat)
  } else {
    out.list <- list("p.value" =  NULL, "norm.p.value" =  NULL,"estimates" = NULL, "transformed_estimates"= NULL)
  }
  return(out.list)
}



normalized_constant_curvature_test_seq <- function(estimates, num.midpoints = 3,
                                                   tri.const.seq = sqrt(sqrt(2)),
                                                   curve.scale = 10,
                                                   heavy.tail.scale = F){
  D <- estimates$D
  mid.search <- estimates$midpoints
  K = nrow(D)
  trim.kappa.sets <- list()
  kappa.sets <- list()
  x.sets <- list()
  index.sets <- list()

  if(nrow(mid.search) < num.midpoints) {
    n.test.points <- nrow(mid.search)
  } else {
    n.test.points <- num.midpoints
  }

  for(k in seq(n.test.points)){

    y.opt = mid.search[k,1]
    z.opt = mid.search[k,2]
    m.opt = mid.search[k,3]
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set<- filter_indices_2(D, y.opt,
                             z.opt,m.opt,
                             tri.const = 0.9)
    #x.set <- 1:K
    #x.set <- x.set[!x.set %in% c(y.opt,z.opt,m.opt)]
    x.sets[[k]] <- x.set
    if(length(x.set) > 1){
      kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
      #kappa.set <- kappa.set[!is.na(kappa.set)]
      #y = scale_curvature(kappa.set, curve.scale)
      # transformation scale
      #y = (y - median(y))/mad(y) + median(y)
      #trim.kappa.sets[[k]] = y
      kappa.sets[[k]] = kappa.set
      index.sets[[k]] = rep(k, length(kappa.set))

    }
  }

  if(num.midpoints <= 1){
    index <- unlist(index.sets)
    #plot(index, unlist(kappa.sets))
    tri.length <- length(tri.const.seq)
    long.out.list <- list()
    for(j in seq(tri.length)){
      tri.const <- tri.const.seq[j]
      trim.kappa.sets.sub <- list()
      kappa.sets.sub <- list()
      index.sets.sub <- list()

      y.opt = mid.search[1,1]
      z.opt = mid.search[1,2]
      m.opt = mid.search[1,3]
      # x.set <- filter_indices(D.hat, y.opt,
      #                         z.opt,m.opt,
      #                         c1 = c1,c2 = c2,c3 = c3)
      x.set.reduced <- filter_indices_2(D, y.opt,
                                        z.opt,m.opt,
                                        tri.const = tri.const)
      x.set <- x.sets[[1]]
      idx.sub <- which(x.set %in% x.set.reduced)
      if(length(idx.sub) == 0){
        next
      }
      x.sub <- x.set[idx.sub]
      kappa.set.sub <- kappa.sets[[1]][idx.sub]
      kappa.set.sub <- kappa.set.sub[!is.na(kappa.set.sub)]

      y = scale_curvature(kappa.set.sub, curve.scale)
      # transformation scale
      if(length(y) > 0 ){

        # q.y <- quantile(y, c(.1,.9))
        # scale.y <- as.numeric(q.y[2] - q.y[1])
        # if(scale.y > 0){
        #   y = (y - median(y))/scale.y + median(y)
        # } else {
        #   y = rep(median(y), length(y))
        # }

        if(mad(y) > 0){
          y = (y - median(y))/mad(y) + median(y)
        } else {
          y = rep(median(y), length(y))
        }

        trim.kappa.sets.sub[[1]] = y
        kappa.sets.sub[[1]] = kappa.set.sub
        index.sets.sub[[1]] = rep(1, length(kappa.set.sub))
      }

      trim.kappa.vec <- unlist(trim.kappa.sets.sub)
      kappa.vec <- unlist(kappa.sets.sub)
      index.vec <- unlist(index.sets.sub)
      trim.test.dat <- data.frame("loc" = index.vec, "trim.est" = trim.kappa.vec,"est" = kappa.vec)
      out.list <- list("p.value" =  1,"norm.p.value" = 1,
                       "estimates" = trim.test.dat, "tri.const" = tri.const)

      long.out.list[[j]] <- out.list
    }
  } else {
    index <- unlist(index.sets)
    #plot(index, unlist(kappa.sets))
    tri.length <- length(tri.const.seq)
    long.out.list <- list()
    for(j in seq(tri.length)){
      tri.const <- tri.const.seq[j]
      trim.kappa.sets.sub <- list()
      kappa.sets.sub <- list()
      index.sets.sub <- list()


      for(k in seq(n.test.points)){

        y.opt = mid.search[k,1]
        z.opt = mid.search[k,2]
        m.opt = mid.search[k,3]
        # x.set <- filter_indices(D.hat, y.opt,
        #                         z.opt,m.opt,
        #                         c1 = c1,c2 = c2,c3 = c3)
        x.set.reduced <- filter_indices_2(D, y.opt,
                                          z.opt,m.opt,
                                          tri.const = tri.const)
        x.set <- x.sets[[k]]
        idx.sub <- which(x.set %in% x.set.reduced)
        if(length(idx.sub) == 0){
          next
        }
        x.sub <- x.set[idx.sub]
        kappa.set.sub <- kappa.sets[[k]][idx.sub]
        kappa.set.sub <- kappa.set.sub[!is.na(kappa.set.sub)]

        y = scale_curvature(kappa.set.sub, curve.scale)
        # transformation scale
        if(length(y) > 0 ){

          # q.y <- quantile(y, c(.1,.9))
          # scale.y <- as.numeric(q.y[2] - q.y[1])
          # if(scale.y > 0){
          #   y = (y - median(y))/scale.y + median(y)
          # } else {
          #   y = rep(median(y), length(y))
          # }

          if(mad(y) > 0){
            y = (y - median(y))/mad(y) + median(y)
          } else {
            y = rep(median(y), length(y))
          }

          trim.kappa.sets.sub[[k]] = y
          kappa.sets.sub[[k]] = kappa.set.sub
          index.sets.sub[[k]] = rep(k, length(kappa.set.sub))
        }


      }

      trim.kappa.vec <- unlist(trim.kappa.sets.sub)
      kappa.vec <- unlist(kappa.sets.sub)
      index.vec <- unlist(index.sets.sub)
      if(length(unique(index.vec)) > 1){
        at.least.two.groups <- length(as.numeric(table(index.vec))) > 1
        if(at.least.two.groups){
          trim.test.dat <- data.frame("loc" = index.vec, "trim.est" = trim.kappa.vec,"est" = kappa.vec)
          # normalized test
          norm.test <- tryCatch({
            kruskal.test(trim.est ~ loc, data = trim.test.dat)
          }, error = function(error_condition) {
            return(NULL)
          })
          test <- tryCatch({
            kruskal.test(est ~ loc, data = trim.test.dat)
          }, error = function(error_condition) {
            return(NULL)
          })
          if(is.null(norm.test) | is.null(test)){
            out.list <- list("p.value" =  NULL,"norm.p.value" =NULL, "estimates" = NULL, "tri.const" = tri.const)
          } else {
            out.list <- list("p.value" =  test$p.value,"norm.p.value" =norm.test$p.value,
                             "estimates" = trim.test.dat, "tri.const" = tri.const)
          }

        }
      } else {
        out.list <- list("p.value" =  NULL,"norm.p.value" =NULL, "estimates" = NULL, "tri.const" = tri.const)
      }
      long.out.list[[j]] <- out.list
    }
  }



  return(long.out.list)
}



test_multi_network_constant_curve <- function(estimates, num.midpoints = 3,
                                              tri.const.seq = sqrt(2),
                                              curve.scale = 10){


  D <- estimates$D
  mid.search <- estimates$midpoints
  K = nrow(D)
  trim.kappa.sets <- list()
  kappa.sets <- list()
  x.sets <- list()
  index.sets <- list()
  for(k in seq(num.midpoints)){
    y.opt = mid.search[k,1]
    z.opt = mid.search[k,2]
    m.opt = mid.search[k,3]
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set<- filter_indices_2(D, y.opt,
                             z.opt,m.opt,
                             tri.const = 0.9)
    #x.set <- 1:K
    #x.set <- x.set[!x.set %in% c(y.opt,z.opt,m.opt)]
    x.sets[[k]] <- x.set
    if(length(x.set) > 1){
      kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
      #kappa.set <- kappa.set[!is.na(kappa.set)]
      #y = scale_curvature(kappa.set, curve.scale)
      # transformation scale
      #y = (y - median(y))/mad(y) + median(y)
      #trim.kappa.sets[[k]] = y
      kappa.sets[[k]] = kappa.set
      index.sets[[k]] = rep(k, length(kappa.set))

    }
  }

  index <- unlist(index.sets)
  tri.length <- length(tri.const.seq)
  long.out.list <- list()
  for(j in seq(tri.length)){
    tri.const <- tri.const.seq[j]
    trim.kappa.sets.sub <- list()
    kappa.sets.sub <- list()
    index.sets.sub <- list()
    for(k in seq(num.midpoints)){

      y.opt = mid.search[k,1]
      z.opt = mid.search[k,2]
      m.opt = mid.search[k,3]
      # x.set <- filter_indices(D.hat, y.opt,
      #                         z.opt,m.opt,
      #                         c1 = c1,c2 = c2,c3 = c3)
      x.set.reduced <- filter_indices_2(D, y.opt,
                                        z.opt,m.opt,
                                        tri.const = tri.const)
      x.set <- x.sets[[k]]
      idx.sub <- which(x.set %in% x.set.reduced)
      x.sub <- x.set[idx.sub]
      kappa.set.sub <- kappa.sets[[k]][idx.sub]
      kappa.set.sub <- kappa.set.sub[!is.na(kappa.set.sub)]
      y = scale_curvature(kappa.set.sub, curve.scale)
      # transformation scale
      y = (y - median(y))/mad(y) + median(y)
      trim.kappa.sets[[k]] = y
      kappa.sets.sub[[k]] = kappa.set.sub
      index.sets.sub[[k]] = rep(k, length(kappa.set.sub))


    }

    trim.kappa.vec <- unlist(trim.kappa.sets)
    kappa.vec <- unlist(kappa.sets.sub)
    index.vec <- unlist(index.sets.sub)
    if(length(unique(index.vec)) > 1){
      at.least.two.groups <- min(as.numeric(table(index.vec))) > 1
      if(at.least.two.groups){
        trim.test.dat <- data.frame("loc" = index.vec, "trim.est" = trim.kappa.vec,"est" = kappa.vec)
        # normalized test
        norm.test <- tryCatch({
          kruskal.test(trim.est ~ loc, data = trim.test.dat)
        }, error = function(error_condition) {
          return(NULL)
        })
        test <- tryCatch({
          kruskal.test(est ~ loc, data = trim.test.dat)
        }, error = function(error_condition) {
          return(NULL)
        })
        if(is.null(norm.test) | is.null(test)){
          out.list <- list("p.value" =  NULL,"norm.p.value" =NULL, "estimates" = NULL, "tri.const" = tri.const)
        } else {
          out.list <- list("p.value" =  test$p.value,"norm.p.value" =norm.test$p.value,
                           "estimates" = trim.test.dat, "tri.const" = tri.const)
        }

      }
    } else {
      out.list <- list("p.value" =  NULL,"norm.p.value" =NULL, "estimates" = NULL, "tri.const" = tri.const)
    }
    long.out.list[[j]] <- out.list
  }

  return(long.out.list)
}
# bootstrap_constant_curvature_test <- function(estimates, num.midpoints = 3, tri.const = sqrt(sqrt(2))){
#   D <- estimates$D
#   mid.search <- estimates$midpoints
#
#   kappa.sets <- list()
#   index <- c()
#   for(k in seq(num.midpoints)){
#     y.opt = mid.search[k,1]
#     z.opt = mid.search[k,2]
#     m.opt = mid.search[k,3]
#     # x.set <- filter_indices(D.hat, y.opt,
#     #                         z.opt,m.opt,
#     #                         c1 = c1,c2 = c2,c3 = c3)
#     x.set <- filter_indices_2(D, y.opt,
#                               z.opt,m.opt,
#                               tri.const = tri.const)
#     kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
#     kappa.set <- kappa.set[!is.na(kappa.set)]
#     kappa.sets[[k]] = kappa.set
#     index <- c(index, rep(k, length(kappa.set)))
#   }
#
#   kappa.vec <- unlist(kappa.sets)
#   overall.med <- median(kappa.vec)
#   if(length(unique(index)) > 1 ){
#
#     boot.med <- c()
#     for(k in seq(num.midpoints)){
#       kappa.block <- kappa.sets[[k]]
#       block.med[k] <- median(kappa.block)
#     }
#     dev.block <- abs(outer(block.med,block.med, "-" ))
#     test.stat <- sum(dev.block[lower.tri(dev.block)])
#     offset.deviation <- block.med - overall.med
#
#     dev.boot <- c()
#     for(b in seq(B)){
#       boot.med <- c()
#       for(k in seq(num.midpoints)){
#         # includes an offset term to adjust for the
#         # fact that we resample.
#         kappa.block <- kappa.sets[[k]] - offset.deviation[[k]]
#         block.size <- length(kappa.block)
#         boot.idx <- sample(seq(block.size), block.size, replace = T)
#         kappa.boot <- kappa.block[boot.idx]
#         boot.med[k] <- median(kappa.boot)
#
#       }
#       dev.block <- abs(outer(boot.med,boot.med, "-" ))
#       dev.boot[b] <- sum(dev.block[lower.tri(dev.block)])
#     }
#     test.dat <- data.frame("loc" = index, "est" = kappa.vec)
#     p.val <- mean(test.stat > dev.boot)
#     out.list <- list("p.value" =  p.val, "estimates" = test.dat)
#   } else {
#     out.list <- list("p.value" =  NULL, "estimates" = NULL)
#   }
#
#   return(out.list)
# }


mv_constant_curvature_test <- function(estimate.set, tri.const = sqrt(sqrt(2))){

  K <- length(estimate.set)
  kappa.sets <- list()
  index <- c()
  for(k in seq(K)){
    estimates <- estimate.set[[k]]
    D <- estimates$D
    mid.search <- estimates$midpoints
    y.opt = mid.search[1,1]
    z.opt = mid.search[1,2]
    m.opt = mid.search[1,3]
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set <- filter_indices_2(D, y.opt,
                              z.opt,m.opt,
                              tri.const = tri.const)
    kappa.set <- estimate_kappa_set(D,y.opt,z.opt,m.opt,x.set)
    kappa.set <- kappa.set[!is.na(kappa.set)]
    kappa.sets[[k]] = kappa.set
    index <- c(index, rep(k, length(kappa.set)))
  }

  kappa.vec <- unlist(kappa.sets)

  test.dat <- data.frame("loc" = index, "est" = kappa.vec)
  test <- kruskal.test(est ~ loc, data = test.dat)
  out.list <- list("p.value" =  test$p.value, "estimates" = test.dat)
  return(out.list)
}




# estimate_kappa_median <- function(D,y,z,m,x.set){
#   # ensure x's exist on the right side
#   # we should have the distance matrix always ensure
#   # y = 1
#   # z = 2
#   # m = 3
#   d1 = nrow(D)
#   d2 = ncol(D)
#   if(d1 > d2){
#     D = t(D)
#   }
#
#   kappa.set <- sapply(x.set, function(x){
#     dxy = D[y,x]
#     dxz = D[z,x]
#     dyz = D[y,z]
#     dxm = D[m,x]
#     d.vec = c(dxy,dxz,dyz,dxm)
#     if(any(d.vec == Inf)){
#       kappa.hat.x <- NA
#     } else{
#       kappa.hat.x <- estimate_kappa(dxy,dxz,dyz,dxm)
#     }
#     return(kappa.hat.x)
#   })
#   kappa.med = median(kappa.set, na.rm = T)
#   return(kappa.med)
# }


estimate_kappa_set <- function(D,y,z,m,x.set){
  # ensure x's exist on the right side
  # we should have the distance matrix always ensure
  # y = 1
  # z = 2
  # m = 3
  d1 = nrow(D)
  d2 = ncol(D)
  if(d1 > d2){
    D = t(D)
  }

  kappa.set <- sapply(x.set, function(x){
    dxy = D[y,x]
    dxz = D[z,x]
    dyz = D[y,z]
    dxm = D[m,x]
    d.vec = c(dxy,dxz,dyz,dxm)
    if(any(d.vec == Inf)){
      kappa.hat.x <- NA
    } else{
      kappa.hat.x <- estimate_kappa(dxy,dxz,dyz,dxm)
    }
    return(kappa.hat.x)
  })

  return(kappa.set)
}


#
# estimate_kappa_mean <- function(D,y,z,m,x.set){
#   # ensure x's exist on the right side
#   # we should have the distance matrix always ensure
#   # y = 1
#   # z = 2
#   # m = 3
#   d1 = nrow(D)
#   d2 = ncol(D)
#   if(d1 > d2){
#     D = t(D)
#   }
#
#   kappa.set <- sapply(x.set, function(x){
#     dxy = D[y,x]
#     dxz = D[z,x]
#     dyz = D[y,z]
#     dxm = D[m,x]
#     d.vec = c(dxy,dxz,dyz,dxm)
#     if(any(d.vec == Inf)){
#       kappa.hat.x <- NA
#     } else{
#       kappa.hat.x <- estimate_kappa(dxy,dxz,dyz,dxm)
#     }
#     return(kappa.hat.x)
#   })
#   kappa.mean = mean(kappa.set, na.rm = T)
#   return(kappa.mean)
# }

# simulate adjacency matrix in the smallest possible manner
# I.e. eliminate all other information
# outputs a bootstrapped distance matrix
# lean_ls_sim <- function(rand.eff.mat, D){
#   # check that dimensions of rand.eff and D match, and check the positivity conditions
#   max.dist = 10**3 # cap on distances
#   # we can use a for loop to initialize the tensor
#
#   ell = ncol(rand.eff.mat)
#   K = nrow(rand.eff.mat)
#   rand.eff.tensor = array(NA, dim = c(3,K,ell^2))
#   D.tensor = array(NA, dim = c(3,K,ell^2))
#
#   # these for loops are only over l arrays
#   # TODO: can we improve the efficiency here
#   # in general, ell^2 won't be too big so this shouldn't be too
#   # worrysome
#   for(i in 1:ell){
#     for(j in 1:ell){
#       slice.idx = ell*(i - 1) + j
#       slice.mat = outer(rand.eff.mat[1:3,i],rand.eff.mat[,j], "+")
#       rand.eff.tensor[,,slice.idx] = slice.mat
#       D.tensor[,,slice.idx] = D[1:3,]
#     }
#   }
#
#
#   P.tensor = exp(rand.eff.tensor - D.tensor)
#
#   X.sim <- apply(P.tensor, 1:3, function(p){
#     out <- rbinom(1,1,p)
#     return(out)
#   })
#   # sum along indices
#   P.sim = apply(X.sim, c(1,2), mean)
#
#
#
#   avg.rand.eff = apply(exp(rand.eff.tensor), c(1,2), mean)
#
#   D.sim = log(avg.rand.eff) - log(P.sim)
#   #
#   D.sim[!is.finite(D.sim)] <- max.dist
#   # adding a maximum distance value for the some occurances which give a negative distance
#   # matrix.
#   # returns the minimal element for the simulated distance matrix
#   D.sim[2,1] = NA
#   D.sim[3,2] = NA
#   D.sim[3,1] = NA
#   D.sim[1,1] = NA
#   D.sim[2,2] = NA
#   D.sim[3,3] = NA # replacing the redundant connections
#   return(D.sim)
# }


# TODO: matrix must be reordered accounting for y,z,m indices.
# TODO: Naming consistency
# lean_network_bootstrap <- function(D,nus.hat.mat,
#                                    B = 1000,y,z,m,x.set){
#
#   K <- length(x.set) + 3
#   l <- ncol(nus.hat.mat)
#
#
#   # if true distances this will just be the true value
#   kappa.med <- estimate_kappa_median(D,y,
#                                      z,m,x.set)
#   kappa.med.boot.vec <- c()
#   for(b in 1:B){
#     D.boot <- lean_ls_sim(nus.hat.mat,D)
#
#
#     # to filter or not to filter
#     # that is the question, (i.e. filter before, or
#     # bootstrap the whole thing)
#     #kappa.med.boot <- filter_median(D.boot, y.opt,z.opt,m.opt,
#     #                                c1 = c1,c2 = c2, c3 = c3)
#     kappa.med.boot <- estimate_kappa_median(D.boot,y,z,m,x.set)
#     kappa.med.boot.vec[b] = kappa.med.boot
#     if(b %% 100 == 0){
#       cat("Bootstrap ", b, "/", B, end = "\r")
#     }
#
#   }
#   # Is this going to be faster?
#
#   #kappa.med.boot.list <- lapply(1:B, function(z){
#     # sample from fitted model
#     #D.boot <- lean_ls_sim(nus.hat.mat,D)
#
#
#     # to filter or not to filter
#     # that is the question, (i.e. filter before, or
#     # bootstrap the whole thing)
#     #kappa.med.boot <- filter_median(D.boot, y.opt,z.opt,m.opt,
#     #                                c1 = c1,c2 = c2, c3 = c3)
#     #kappa.med.boot <- estimate_kappa_median(D.boot,y,z,m,x.set)
#     #return(kappa.med.boot)
#   #})
#   #kappa.med.boot.vec = unlist(kappa.med.boot.list)
#
#   se.est = sqrt(mean((kappa.med.boot.vec - kappa.med)^2, na.rm = T))
#   sd.est = sqrt(mean((kappa.med.boot.vec - mean(kappa.med.boot.vec, na.rm = T))^2, na.rm = T))
#   bias.est = (mean(kappa.med.boot.vec, na.rm = T)) - kappa.med
#
#   out.list = list("bias" = bias.est, "se" = se.est, "sd" = sd.est)
#   return(out.list)
# }

#
# filter_median <- function(D, y,z,m, c1 = 1/2,c2 = 2,c3 = 4, min.K = 3){
#   K = nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(y,z,m)]
#
#   #indicator of whether to include an x
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
#       return(F)
#     } else {
#       i1 = (c1)*D[y,z] <= D[x,y]
#       i2 = (c1)*D[y,z] <= D[x,z]
#       i3 = D[x,z] <= (c2)*D[y,z]
#       i4 = D[x,y] <= (c2)*D[y,z]
#       i5 = abs(D[x,y]^2 - D[x,z]^2) <= c3*(D[x,y] + D[x,z] - D[y,z])
#       return(ifelse(i1*i2*i3*i4*i5 == 1, T, F))
#     }
#
#   })
#
#   x.filtered <- x.set[x.id]
#   if(length(x.filtered) < min.K){
#     kappa.med <- NA
#   } else {
#     kappa.set <- sapply(x.filtered, function(x){
#       dxy = D.hat[x,y.opt]
#       dxz = D.hat[x,z.opt]
#       dyz = D.hat[y.opt,z.opt]
#       dxm = D.hat[x,m.opt]
#       d.vec = c(dxy,dxz,dyz,dxm)
#       if(any(d.vec == Inf)){
#         kappa.hat.x <- NA
#       } else{
#         kappa.hat.x <- estimate_kappa(dxy,dxz,dyz,dxm)
#       }
#       return(kappa.hat.x)
#     })
#     kappa.med = median(kappa.set, na.rm = T)
#   }
#   return(kappa.med)
# }


cm_f <- function(d,kappa){
  if(kappa == 0){
    return(d**2)
  } else if(kappa > 0){
    return(2*(1/kappa)*(1 - cos(d*sqrt(kappa))))
  } else if(kappa < 0){
    return(2*(1/-kappa)*(cosh(d*sqrt(-kappa)) - 1))
  }
}

# gives the Cayley menger determinant (up to a multiplicative factor)
cm_determinant <- function(D, kappa){
  K = nrow(D)
  if(kappa == 0){
    M = D**2 # elementwise square
    M = rbind(M, rep(1,K))
    M = cbind(M, c(rep(1,K),0))
  } else if(kappa > 0){
    M = cm_f(D,kappa)
    M = rbind(M, rep(1,K))
    M = cbind(M, c(rep(1,K),kappa/2))
  } else if(kappa < 0){
    M = cm_f(D,kappa)
    M = rbind(M, rep(1,K))
    M = cbind(M, c(rep(1,K),-kappa/2))
  }
  n = nrow(D)
  const <- (((-1)**(n))/(((factorial(n - 1)**2)**2)*(2**(n - 1))))
  out <- det(sqrt(sqrt(sqrt(const)))*M)
  out <- det(const*M)
  return(out)
}




#
# oracle_filter <- function(Z, kappa.true, theta.max, d.max, d.min){
#   D = pos_to_dist(Z,kappa.true)
#   # only true when midpoint is the third index
#   dvec <- D[3,]
#   K <- nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(1,2,3)]
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     i1 = d.max >= dvec[x]
#     i2 = d.min <= dvec[x]
#     i3 = abs(atan(Z[x,2]/Z[x,3])) <= theta.max
#     return(ifelse(i1*i2*i3 == 1, T, F))
#   })
#
#   x.filtered <- x.set[x.id]
#   print(paste("Fraction filtered: ", round(mean(x.id),3)))
#   return(x.filtered)
# }

#
# angle_filter <- function(Z, kappa.true, theta.max, d.max, d.min){
#   D = pos_to_dist(Z,kappa.true)
#   # only true when midpoint is the third index
#   dvec <- D[3,]
#   K <- nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(1,2,3)]
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     i1 = d.max >= dvec[x]
#     i2 = d.min <= dvec[x]
#     i3 = abs(atan(Z[x,2]/Z[x,3])) <= theta.max
#     return(ifelse(i1*i2*i3 == 1, T, F))
#   })
#
#   x.filtered <- x.set[x.id]
#   return(x.filtered)
# }

# jacknife for a single set
#
# jacknife_estimate <- function(A,nus.hat.mat,clique.legend,
#                           B = 1000,y,z,m,x){
#   A <- sim_ls_network(rand.eff,D)
#   D.est = D_estimate(A,clique.legend,fixed.effects = nus.hat.mat)
#   D.hat = D.est$estimates
#   triang <- c(y,z,m,x)
#   D.hat.small <- D.hat[c(y,z,m,x),c(y,z,m,x)]
#   nus.hat.mat.small <- nus.hat.mat[c(y,z,m,x),]
#   clique.subset <- clique.legend[c(y,z,m,x),]
#
#
#   kappa.hat <- estimate_kappa(D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#
#   K <- nrow(clique.legend)
#   ell <- ncol(clique.legend)
#
#   # permute the subset randomly
#   permutations <- t(sapply(1:4, function(j){
#     return(sample(1:ell,ell))
#   }))
#
#   permuted.subset <- t(sapply(1:4, function(j){
#     return(clique.subset[j,permutations[j,]])
#   }))
#
#   permuted.nus <- t(sapply(1:4, function(j){
#     return(nus.hat.mat.small[j,permutations[j,]])
#   }))
#
#   # rearranging A to the proper subset:
#
#   # TODO: check that this is in fact using the correct rows
#   # TODO: Unit tests
#   A.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   nu.sum.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   pair.set <- matrix(NA,nrow = 4, ncol = 2)
#   pair.set[1,] <- c(x,y)
#   pair.set[2,] <- c(x,z)
#   pair.set[3,] <- c(y,z)
#   pair.set[4,] <- c(x,m)
#   for(j in 1:4){
#     # expands the grid and allows the
#     # number of points to be selected
#     grid.idx <- expand.grid(permuted.subset[pair.set[j,1],],permuted.subset[pair.set[j,2],])
#
#     a.sub.vec <- apply(grid.idx,1, function(q){
#       A[q[1],q[2]]})
#     #a.sub.vec <- as.vector(A[clique.subset[pair.set[j,1],],clique.subset[pair.set[j,2],]])
#     A.subset[j,] <- a.sub.vec
#
#     nu.sum.sub.vec <- rowSums(expand.grid(permuted.nus[pair.set[j,1],],permuted.nus[pair.set[j,2],]))
#     nu.sum.subset[j,] <- nu.sum.sub.vec
#   }
#
#   jack.kappas <- sapply(1:ell**2, function(j){
#     idx <- 1:(ell**2)
#     idx <- idx[-j]
#     p.hat.jack <- rowMeans(A.subset[,idx])
#     nu.exp.mean.jack <- rowMeans(exp(nu.sum.subset[,idx]))
#     d.hat.jack <- log(nu.exp.mean.jack) - log(p.hat.jack)
#     print(c(d.hat.jack[1],d.hat.jack[2],d.hat.jack[3],d.hat.jack[4]))
#     kappa.jack <- tryCatch({
#       estimate_kappa(d.hat.jack[1],d.hat.jack[2],d.hat.jack[3],d.hat.jack[4])
#     }, error = function(e) {return(NA)})
#
#     return(kappa.jack)
#   })
#
#
#   # jack.kappas <- sapply(1:ell, function(j){
#   #   idx <- 1:ell
#   #   idx <- idx[-j]
#   #   permuted.subset.jack <- permuted.subset[,idx]
#   #   nus.hat.mat.jack <- permuted.nus[,idx]
#   #   D.jack <- D_estimate(A,permuted.subset.jack,fixed.effects = nus.hat.mat.jack)
#   #   D.jack.hat <- D.jack$estimates
#   #   print(j)
#   #   kappa.jack <- tryCatch({
#   #     estimate_kappa(D.jack.hat[x,y],D.jack.hat[x,z],D.jack.hat[y,z],D.jack.hat[x,m])
#   #   }, error = function(e) {return(NA)})
#   #
#   #   return(kappa.jack)
#   # })
#   #print(jack.kappas)
#   #print(mean(jack.kappas,na.rm = T))
#   bias.est <- ((ell**2 - 1))*(mean(jack.kappas,na.rm = T) - kappa.hat)
#   #TODO: Remove
#   print(kappa.hat)
#   print(mean(jack.kappas,na.rm = T))
#   print(bias.est)
#   #mean(jack.kappas,na.rm = T)
#   return(mean(jack.kappas,na.rm = T))
# }




#
# bootstrap_bias_estimate <- function(A,nus.hat.mat,clique.legend,
#                                     B.boot = 200,y,z,m,x, parallel = F){
#   #A <- sim_ls_network(rand.eff,D)
#   #B.boot <- 200
#   D.est = D_estimate(A,clique.legend,fixed.effects = nus.hat.mat)
#   D.hat = D.est$estimates
#   triang <- c(y,z,m,x)
#   D.hat.small <- D.hat[c(y,z,m,x),c(y,z,m,x)]
#   nus.hat.mat.small <- nus.hat.mat[c(y,z,m,x),]
#   clique.subset <- clique.legend[c(y,z,m,x),]
#
#
#   kappa.hat <- estimate_kappa(D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#
#   K <- nrow(clique.legend)
#   ell <- ncol(clique.legend)
#
#   # permute the subset randomly
#   permutations <- t(sapply(1:4, function(j){
#     return(sample(1:ell,ell))
#   }))
#
#   permuted.subset <- t(sapply(1:4, function(j){
#     return(clique.subset[j,permutations[j,]])
#   }))
#
#   permuted.nus <- t(sapply(1:4, function(j){
#     return(nus.hat.mat.small[j,permutations[j,]])
#   }))
#
#   # rearranging A to the proper subset:
#
#   # TODO: check that this is in fact using the correct rows
#   # TODO: Unit tests
#   A.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   nu.sum.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   pair.set <- matrix(NA,nrow = 4, ncol = 2)
#
#   pair.set[1,] <- c(x,y)
#   pair.set[2,] <- c(x,z)
#   pair.set[3,] <- c(y,z)
#   pair.set[4,] <- c(x,m)
#
#   pair.set[1,] <- c(4,1)
#   pair.set[2,] <- c(4,2)
#   pair.set[3,] <- c(1,2)
#   pair.set[4,] <- c(4,3)
#
#
#   # should we do a parametric bootstrap
#   if(parallel){
#     boot.kappas <- mclapply(1:B.boot, function(b){
#       for(j in 1:4){
#         # expands the grid and allows the
#         # number of points to be selected
#         grid.idx <- expand.grid(permuted.subset[pair.set[j,1],],permuted.subset[pair.set[j,2],])
#         a.sub.vec <- apply(grid.idx,1, function(q){
#           A[q[1],q[2]]})
#         a.sub.vec <- as.vector(A[clique.subset[pair.set[j,1],],clique.subset[pair.set[j,2],]])
#         A.subset[j,] <- a.sub.vec
#
#         nu.sum.sub.vec <- rowSums(expand.grid(nus.hat.mat.small[pair.set[j,1],],nus.hat.mat.small[pair.set[j,2],]))
#         nu.sum.subset[j,] <- nu.sum.sub.vec
#       }
#       idx <- sample(1:ell**2,ell**2, replace = T)
#       p.hat.boot <- rowMeans(A.subset[,idx])
#       #print(p.hat.boot)
#       nu.exp.mean.boot <- rowMeans(exp(nu.sum.subset[,idx]))
#       d.hat.boot <- log(nu.exp.mean.boot) - log(p.hat.boot)
#       #print(c(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4]))
#       kappa.boot <- tryCatch({
#         estimate_kappa(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4])
#       }, error = function(e) {return(NA)})
#       kappa.boot
#     })
#     boot.kappas <- unlist(boot.kappas)
#   } else{
#     boot.kappas <- rep(NA, B.boot)
#     for(b in 1:B.boot){
#       if(b %% 100 == 0){
#         #cat(paste("Boot", b, "/", B.boot), end = "\r")
#       }
#       for(j in 1:4){
#         # expands the grid and allows the
#         # number of points to be selected
#         grid.idx <- expand.grid(permuted.subset[pair.set[j,1],],permuted.subset[pair.set[j,2],])
#
#         a.sub.vec <- apply(grid.idx,1, function(q){
#           A[q[1],q[2]]})
#         a.sub.vec <- as.vector(A[clique.subset[pair.set[j,1],],clique.subset[pair.set[j,2],]])
#         A.subset[j,] <- a.sub.vec
#
#         nu.sum.sub.vec <- rowSums(expand.grid(nus.hat.mat.small[pair.set[j,1],],nus.hat.mat.small[pair.set[j,2],]))
#         nu.sum.subset[j,] <- nu.sum.sub.vec
#       }
#       idx <- sample(1:ell**2,ell**2, replace = T)
#       p.hat.boot <- rowMeans(A.subset[,idx])
#       #print(p.hat.boot)
#       nu.exp.mean.boot <- rowMeans(exp(nu.sum.subset[,idx]))
#       d.hat.boot <- log(nu.exp.mean.boot) - log(p.hat.boot)
#       #print(c(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4]))
#       kappa.boot <- tryCatch({
#         estimate_kappa(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4])
#       }, error = function(e) {return(NA)})
#       boot.kappas[b] <- kappa.boot
#     }
#   }
#
#   kappa.boot.mean <- mean(boot.kappas,na.rm = T)
#
#
#   # jack.kappas <- sapply(1:ell, function(j){
#   #   idx <- 1:ell
#   #   idx <- idx[-j]
#   #   permuted.subset.jack <- permuted.subset[,idx]
#   #   nus.hat.mat.jack <- permuted.nus[,idx]
#   #   D.jack <- D_estimate(A,permuted.subset.jack,fixed.effects = nus.hat.mat.jack)
#   #   D.jack.hat <- D.jack$estimates
#   #   print(j)
#   #   kappa.jack <- tryCatch({
#   #     estimate_kappa(D.jack.hat[x,y],D.jack.hat[x,z],D.jack.hat[y,z],D.jack.hat[x,m])
#   #   }, error = function(e) {return(NA)})
#   #
#   #   return(kappa.jack)
#   # })
#   #print(jack.kappas)
#   #print(mean(jack.kappas,na.rm = T))
#   bias.est <- (kappa.boot.mean - kappa.hat)
#   #TODO: Remove
#   # print(kappa.hat)
#   #print(kappa.boot.mean)
#   # print(bias.est)
#   #mean(jack.kappas,na.rm = T)
#   bias.corrected.estimate <- kappa.hat - bias.est
#   return(bias.est)
# }

#
# parametric_bootstrap_dev_estimate <- function(A,nus.hat.mat,clique.legend,
#                                               B.boot = 200,y,z,m,x, parallel = F, verbose = F){
#   #A <- sim_ls_network(rand.eff,D)
#   #B.boot <- 200
#   D.est = D_estimate(A,clique.legend,fixed.effects = nus.hat.mat)
#   D.hat = D.est$estimates
#   triang <- c(y,z,m,x)
#   D.hat.small <- D.hat[c(y,z,m,x),c(y,z,m,x)]
#   nus.hat.mat.small <- nus.hat.mat[c(y,z,m,x),]
#   clique.subset <- clique.legend[c(y,z,m,x),]
#
#
#   kappa.hat <- estimate_kappa(D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#   #print(kappa.hat)
#   K <- nrow(clique.legend)
#   ell <- ncol(clique.legend)
#
#   # permute the subset randomly
#   permutations <- t(sapply(1:4, function(j){
#     return(sample(1:ell,ell))
#   }))
#
#   permuted.subset <- t(sapply(1:4, function(j){
#     return(clique.subset[j,permutations[j,]])
#   }))
#
#   permuted.nus <- t(sapply(1:4, function(j){
#     return(nus.hat.mat.small[j,permutations[j,]])
#   }))
#
#   # rearranging A to the proper subset:
#
#   # TODO: check that this is in fact using the correct rows
#   # TODO: Unit tests
#   A.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   nu.sum.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   pair.set <- matrix(NA,nrow = 4, ncol = 2)
#
#   pair.set[1,] <- c(x,y)
#   pair.set[2,] <- c(x,z)
#   pair.set[3,] <- c(y,z)
#   pair.set[4,] <- c(x,m)
#
#   pair.set[1,] <- c(4,1)
#   pair.set[2,] <- c(4,2)
#   pair.set[3,] <- c(1,2)
#   pair.set[4,] <- c(4,3)
#
#
#   for(j in 1:4){
#     # expands the grid and allows the
#     # number of points to be selected
#     grid.idx <- expand.grid(permuted.subset[pair.set[j,1],],permuted.subset[pair.set[j,2],])
#     a.sub.vec <- apply(grid.idx,1, function(q){
#       A[q[1],q[2]]})
#     a.sub.vec <- as.vector(A[clique.subset[pair.set[j,1],],clique.subset[pair.set[j,2],]])
#     A.subset[j,] <- a.sub.vec
#
#     nu.sum.sub.vec <- rowSums(expand.grid(nus.hat.mat.small[pair.set[j,1],],nus.hat.mat.small[pair.set[j,2],]))
#     nu.sum.subset[j,] <- nu.sum.sub.vec
#   }
#   d.hat.vec <- c( D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#   p.hat.subset <- exp(-d.hat.vec + nu.sum.subset)
#   vec.p.hat.subset <- as.vector(p.hat.subset)
#   # should we do a parametric bootstrap
#
#   boot.kappas <- rep(NA, B.boot)
#   for(b in 1:B.boot){
#     if(verbose & b %% 100 == 0){
#       cat(paste("Boot", b, "/", B.boot), end = "\r")
#     }
#
#     U.boot <- runif(n = length(vec.p.hat.subset))
#     A.boot <- matrix(1*(U.boot <= vec.p.hat.subset), nrow = 4)
#     p.hat.boot <- rowMeans(A.boot)
#     #print(p.hat.boot)
#     nu.exp.mean.boot <- rowMeans(exp(nu.sum.subset))
#     d.hat.boot <- log(nu.exp.mean.boot) - log(p.hat.boot)
#     #print(c(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4]))
#     kappa.boot <- tryCatch({
#       estimate_kappa(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4])
#     }, error = function(e) {return(NA)})
#     boot.kappas[b] <- kappa.boot
#   }
#
#
#   kappa.boot.med <- median(boot.kappas,na.rm = T)
#   dev.est <- (kappa.boot.med - kappa.hat)
#   #TODO: Remove
#   # print(kappa.hat)
#   #print(kappa.boot.mean)
#   # print(bias.est)
#   #mean(jack.kappas,na.rm = T)
#   dev.corrected.estimate <- kappa.hat - dev.est
#   return(dev.est)
# }




#
# parametric_double_bootstrap_dev_cor_estimate <- function(A,nus.hat.mat,clique.legend,
#                                                          B.boot1 = 200,B.boot2 = 200,
#                                                          y,z,m,x, parallel = F, verbose = T){
#
#   D.est = D_estimate(A,clique.legend,fixed.effects = nus.hat.mat)
#   D.hat = D.est$estimates
#   triang <- c(y,z,m,x)
#   D.hat.small <- D.hat[c(y,z,m,x),c(y,z,m,x)]
#   nus.hat.mat.small <- nus.hat.mat[c(y,z,m,x),]
#   clique.subset <- clique.legend[c(y,z,m,x),]
#
#
#   kappa.hat <- estimate_kappa(D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#   #print(kappa.hat)
#   K <- nrow(clique.legend)
#   ell <- ncol(clique.legend)
#
#   # permute the subset randomly
#   permutations <- t(sapply(1:4, function(j){
#     return(sample(1:ell,ell))
#   }))
#
#   permuted.subset <- t(sapply(1:4, function(j){
#     return(clique.subset[j,permutations[j,]])
#   }))
#
#   permuted.nus <- t(sapply(1:4, function(j){
#     return(nus.hat.mat.small[j,permutations[j,]])
#   }))
#
#   # rearranging A to the proper subset:
#
#   # TODO: check that this is in fact using the correct rows
#   # TODO: Unit tests
#   A.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   nu.sum.subset <- matrix(NA, nrow = 4, ncol = ell**2)
#   pair.set <- matrix(NA,nrow = 4, ncol = 2)
#
#   pair.set[1,] <- c(x,y)
#   pair.set[2,] <- c(x,z)
#   pair.set[3,] <- c(y,z)
#   pair.set[4,] <- c(x,m)
#
#   pair.set[1,] <- c(4,1)
#   pair.set[2,] <- c(4,2)
#   pair.set[3,] <- c(1,2)
#   pair.set[4,] <- c(4,3)
#
#
#   for(j in 1:4){
#     # expands the grid and allows the
#     # number of points to be selected
#     grid.idx <- expand.grid(permuted.subset[pair.set[j,1],],permuted.subset[pair.set[j,2],])
#     a.sub.vec <- apply(grid.idx,1, function(q){
#       A[q[1],q[2]]})
#     a.sub.vec <- as.vector(A[clique.subset[pair.set[j,1],],clique.subset[pair.set[j,2],]])
#     A.subset[j,] <- a.sub.vec
#
#     nu.sum.sub.vec <- rowSums(expand.grid(nus.hat.mat.small[pair.set[j,1],],nus.hat.mat.small[pair.set[j,2],]))
#     nu.sum.subset[j,] <- nu.sum.sub.vec
#   }
#   d.hat.vec <- c( D.hat[x,y],D.hat[x,z],D.hat[y,z],D.hat[x,m])
#   p.hat.subset <- exp(-d.hat.vec + nu.sum.subset)
#   vec.p.hat.subset <- as.vector(p.hat.subset)
#   # should we do a parametric bootstrap
#
#   boot.kappas1 <- rep(NA, B.boot1)
#   boot.kappas2 <- matrix(NA, nrow = B.boot1, ncol = B.boot2)
#   # double bootstrap for a higher order bias correction.
#   for(b1 in 1:B.boot1){
#     if(verbose & b1 %% 100 == 0){
#       cat(paste("Outer Boot", b1, "/", B.boot1), end = "\n")
#     }
#
#     U.boot1 <- runif(n = length(vec.p.hat.subset))
#     A.boot1 <- matrix(1*(U.boot1 <= vec.p.hat.subset), nrow = 4)
#     p.hat.boot1 <- rowMeans(A.boot1)
#     #print(p.hat.boot)
#     nu.exp.mean.boot1 <- rowMeans(exp(nu.sum.subset))
#     d.hat.boot1 <- log(nu.exp.mean.boot1) - log(p.hat.boot1)
#     #print(c(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4]))
#     # creating a new array using the well estimated random effects and the
#     # computed distance matrix
#     #
#     p.hat.subset.boot1 <- exp(-d.hat.boot1 + nu.sum.subset)
#     vec.p.hat.subset.boot1 <- as.vector(p.hat.subset.boot1)
#     kappa.boot1 <- tryCatch({
#       estimate_kappa(d.hat.boot1[1],d.hat.boot1[2],d.hat.boot1[3],d.hat.boot1[4])
#     }, error = function(e) {return(NA)})
#
#     boot.kappas1[b1] <- kappa.boot1
#     for(b2 in 1:B.boot2){
#       if(verbose & b2 %% 100 == 0){
#         cat(paste("Inner Boot", b2, "/", B.boot2), end = "\r")
#       }
#
#       U.boot2 <- runif(n = length(vec.p.hat.subset.boot1))
#       A.boot2 <- matrix(1*(U.boot2 <= vec.p.hat.subset.boot1), nrow = 4)
#       p.hat.boot2 <- rowMeans(A.boot2)
#       #print(p.hat.boot)
#       nu.exp.mean.boot2 <- rowMeans(exp(nu.sum.subset))
#       d.hat.boot2 <- log(nu.exp.mean.boot2) - log(p.hat.boot2)
#       #print(c(d.hat.boot[1],d.hat.boot[2],d.hat.boot[3],d.hat.boot[4]))
#       kappa.boot2 <- tryCatch({
#         estimate_kappa(d.hat.boot2[1],d.hat.boot2[2],d.hat.boot2[3],d.hat.boot2[4])
#       }, error = function(e) {return(NA)})
#       #print(kappa.boot2)
#       boot.kappas2[b1,b2] <- kappa.boot2
#     }
#   }
#
#   row.dev <- apply(boot.kappas2, 1, median, na.rm = T) - boot.kappas1
#
#   kappa.row.corrected.kappa.1 <- boot.kappas1 - row.dev
#
#
#   kappa.boot.med <- median(kappa.row.corrected.kappa.1,na.rm = T)
#   dev.est <- (kappa.boot.med - kappa.hat)
#   #TODO: Remove
#   # print(kappa.hat)
#   #print(kappa.boot.mean)
#   # print(bias.est)
#   #mean(jack.kappas,na.rm = T)
#   #dev.corrected.estimate <- 3*kappa.hat - 3*mean(boot.kappas1,na.rm = T) + mean(boot.kappas2,na.rm = T)
#   dev.corrected.estimate <- kappa.hat - dev.est
#   return(dev.corrected.estimate)
# }





# a filtered median estimate of kappa.
# filter_median_bootstrap_debias <- function(A, y,z,m, c1 = 1/2,c2 = 2,c3 = 0.75, min.K = 3){
#   K = nrow(D)
#   x.set <- 1:K
#   x.set <- x.set[-c(y,z,m)]
#
#   #indicator of whether to include an x
#   x.id <- sapply(x.set, function(x){
#     # ensuring that the distance is not infinite to each of x.y.m
#     if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
#       return(F)
#     } else {
#       i1 = (c1)*D[y,z] <= D[x,y]
#       i2 = (c1)*D[y,z] <= D[x,z]
#       i3 = D[x,z] <= (c2)*D[y,z]
#       i4 = D[x,y] <= (c2)*D[y,z]
#       i5 = abs(D[x,y] - D[x,z]) <= 2*c3*(D[y,z])
#       return(ifelse(i1*i2*i3*i4*i5 == 1, T, F))
#     }
#
#   })
#
#   x.filtered <- x.set[x.id]
#
#   kappa.set <- sapply(x.filtered, function(x){
#
#
#
#     dxy = D[x,y]
#     dxz = D[x,z]
#     dyz = D[y,z]
#     dxm = D[x,m]
#     d.vec = c(dxy,dxz,dyz,dxm)
#     if(any(d.vec == Inf)){
#       kappa.hat.x <- NA
#     } else{
#       kappa.hat.x <- parametric_bootstrap_dev_estimate(A,nus.hat.mat,clique.legend,
#                                                        B.boot = 200,y,z,m,x,
#                                                        parallel = F, verbose = F)
#     }
#     return(kappa.hat.x)})
#   kappa.med = median(kappa.set, na.rm = T)
#
#   return(kappa.med)
# }


#
# outer_objective <- function(x){
#   y = 1/x + 1/(1 - x)
#   y = ifelse(abs(x - 1/2) <= 1/2, y,Inf)
#   return(y)
# }
#


scale_curvature <- function(x,c = 1){
  out <- c*(exp(2*x/c) - 1)/(exp(2*x/c) + 1)
  out[is.infinite(x) & x > 0] = c
  return(out)
}





# TODO: Remove time lines
# G <- G.sub
# cliques <- subset.ids
# fixed.effect.vec <- fixed.effect.vec
estimate_D_restricted_fast <- function(G, cliques, fixed.effect.vec, D0,
                                       thresh = 10**(-3), max.iter = 50, solver = "MOSEK",
                                       verbose = F){
  # numerical smoothing for some gradient terms.
  eps = 10**(-8) # precision for terms

  if(solver %in% c("MOSEK","GUROBI")){
    if(! solver %in% installed_solvers()){
      cvx_solver = "OSQP"
    } else {
      cvx_solver = solver
    }
  } else {
    cvx_solver = solver
  }

  K = max(cliques)
  if(missing(D0)){
    D0 <- init_D0(G, cliques, fixed.effect.vec) # provide an initial value
  }


  D <- Variable(K,K, name = "Distance Matrix")

  #D.big <- Variable(n.subg,n.subg, name = "Distance Matrix")
  nu.big <- outer(fixed.effect.vec,fixed.effect.vec, "+")
  #d.vec <- Variable(K^2, name = "Distance Vector")

  d.vec = vec(D) #define a vectorization of the distance matrix
  constraints <- list(
    diag(D) == 0,
    D == t(D),
    D >= 0
  )

  # the triangle inequalities

  index.block <- expand.grid(seq(K),seq(K),seq(K))
  ib.1 <- index.block[,1] < index.block[,2]
  ib.2 <- index.block[,2] < index.block[,3]
  index.block <- index.block[ib.1 & ib.2, ]


  ref.block <- index.block
  ref.block[,1] <- (index.block[,1] - 1)*K + index.block[,2]
  ref.block[,2] <- (index.block[,2] - 1)*K + index.block[,3]
  ref.block[,3] <- (index.block[,1] - 1)*K + index.block[,3]

  E.block.1 <- matrix(0,nrow(index.block),K^2)
  E.block.1 <- as(E.block.1, "sparseMatrix")
  n.ref <- nrow(ref.block)
  E.block.1 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                            j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                            x = c(rep(1,n.ref),rep(1,n.ref),rep(-1,n.ref)),
                            dims = c(n.ref,K^2))
  E.block.2 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                            j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                            x = c(rep(1,n.ref),rep(-1,n.ref),rep(1,n.ref)),
                            dims = c(n.ref,K^2))
  E.block.3 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                            j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                            x = c(rep(-1,n.ref),rep(1,n.ref),rep(1,n.ref)),
                            dims = c(n.ref,K^2))

  E.mat <- rbind(E.block.1,
                 E.block.2,
                 E.block.3)

  rm(E.block.1)
  rm(E.block.2)
  rm(E.block.3)

  # turn the many restrictions into a single one.
  constraints <- append(constraints, E.mat %*% d.vec >= 0)



  ######
  # for fast computation of M and B matrices
  clique.sizes = table(cliques)
  l.max = max(clique.sizes)
  Id.tens = array(0,c(K,K,l.max**2))
  G.tens = array(0,c(K,K,l.max**2)) #tensor version of the clique adjacency
  nu.tens = array(0,c(K,K,l.max**2))

  # each pair has a length
  for(k1 in seq(K)){
    for(k2 in seq(K)){
      if(k1 != k2){
        n.potential.connections = clique.sizes[k1]*clique.sizes[k2]
        idx1 = which(cliques == k1)
        idx2 = which(cliques == k2)
        nu.block <- nu.big[idx1,idx2]

        Id.tens[k1,k2,1:n.potential.connections] = 1
        nu.tens[k1,k2,1:n.potential.connections] = as.vector(nu.block)
        G.tens[k1,k2,1:n.potential.connections] = as.vector(G[idx1,idx2])
        #D.prev.tens[k1,k2,] = D.prev[k1,k2]
      }
    }
  }

  ######
  # successive second order approximation
  # allows for a fast implementation of QP solutions
  D.prev <- D0
  D.prev.prev <- D0
  diff = Inf
  iter = 1
  lik.prev = -Inf
  while(diff > thresh & iter < max.iter){
    time.1 <- Sys.time()
    D.prev.tens = array(0,c(K,K,l.max**2))
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        if(k1 != k2){
          D.prev.tens[k1,k2,] = D.prev[k1,k2]
        }
      }
    }

    M.tens <- ((1 - G.tens)*(-1/2)*(exp(nu.tens - D.prev.tens - eps))/((1 - exp(nu.tens - D.prev.tens - eps))**2))*Id.tens

    B.tens <- (-G.tens + (1 - G.tens)*(exp(nu.tens - D.prev.tens - eps))/(1 - exp(nu.tens - D.prev.tens - eps)) - 2*M.tens*D.prev.tens)*Id.tens


    M = rowSums(M.tens, dims = 2)
    B = rowSums(B.tens, dims = 2)

    diag(M) = 0
    diag(B) = 0

    time.2 <- Sys.time()
    # print(time.2 - time.1)


    b = as.vector(B)
    w = sqrt(as.vector(-M))
    obj.arg = -sum_squares(w*d.vec) + sum(b * d.vec)

    obj <- Maximize(obj.arg)
    prob <- Problem(obj, constraints)

    time.1 <- Sys.time()
    # can get weirdly stuck in numerical error

    result <- tryCatch({
      solve(prob, solver = cvx_solver)
    }, error = function(e) {

      return(NULL)})

    time.2 <- Sys.time()
    #print(time.2 - time.1)
    # unstuck the problem sometimes

    if(!is.null(result)){
      if(result$status == "solver_error"){
        D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
        lik.prev = -Inf
      } else {
        D.next <- result$getValue(D)
      }

    } else {
      D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
      lik.prev <- -Inf
    }

    D.next[D.next < eps] = 0
    # large fraction at zero, we should reset
    if(mean(D.next == 0) > 0.5){
      D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
      lik.prev <- -Inf
    }

    D.next[D.next < eps] = 0
    lik.next = likelihood_value(G, cliques, fixed.effect.vec, D.next)
    diff = lik.next - lik.prev
    if(is.nan(diff) | lik.prev == -Inf){
      diff = Inf
    }
    lik.prev = lik.next
    D.prev.prev = D.prev
    D.prev = D.next

    if(verbose){
      cat(paste("Num Steps:", iter, "Diff Likelihood:", round(diff,4)), end = "\r")
    }
    iter = iter + 1
  }

  D.hat = D.next

  return(D.hat)
}


# Explain this is due to numeric precision
# trimmed version seems to be a good arguement for an initialize
init_D0_old <- function(G, cliques,fixed.effect.vec){
  thresh = 10**(-7)
  K = max(cliques)
  D0 = matrix(0, K,K)
  for(k1 in seq(K)){
    for(k2 in seq(K)){
      idx1 = which(cliques == k1)
      idx2 = which(cliques == k2)
      p.xy = mean(G[idx1,idx2])
      D0[k1,k2] = -log(p.xy) + log(mean(exp(fixed.effect.vec[idx1]))) + log(mean(exp(fixed.effect.vec[idx2])))
    }
  }
  D0 = D0 - min(D0)
  diag(D0) = 0
  trim = 100
  D0[D0 >trim] = trim
  while(max_triangle_deviation(D0) > thresh){
    trim = max(D0) - max_triangle_deviation(D0)
    D0[D0 >trim] = trim
  }
  return(D0)
}


# essentially a variant of the Floyd-Warshal Algorithm with an extra pre-processing step.
# G = G.sub
# cliques = subset.ids
# fixed.effect.vec
init_D0 <- function(G, cliques,fixed.effect.vec, trim.0 = T, global.randeff = T, add.offset = 0.001){
  if(global.randeff){
    K = max(cliques)
    D0 = matrix(0, K,K)
    glb.rndf <- log(mean(exp(fixed.effect.vec)))
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        idx1 = which(cliques == k1)
        idx2 = which(cliques == k2)
        p.xy = mean(G[idx1,idx2])
        D0[k1,k2] = -log(p.xy) + 2*glb.rndf
      }
    }
  } else {
    K = max(cliques)
    D0 = matrix(0, K,K)
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        idx1 = which(cliques == k1)
        idx2 = which(cliques == k2)
        p.xy = mean(G[idx1,idx2])
        D0[k1,k2] = -log(p.xy) + log(mean(exp(fixed.effect.vec[idx1]))) + log(mean(exp(fixed.effect.vec[idx2])))
      }
    }
  }

  D0[D0 < 0 ] = 0
  D0 = fwa(D0)

  D0 = D0 + add.offset
  diag(D0) = 0

  return(D0)
}


fwa <- function(D0.init){
  # first pass D0 may not be a metric
  min.d0 <- min(D0.init)
  D0 = D0.init - min.d0
  diag(D0) = 0
  K <- nrow(D0)
  for(k in seq(K)){
    for(j in seq(K)){
      for(i in seq(j)){
        if(D0[i,j] > D0[i,k] + D0[k,j]){
          D0[i,j] = D0[i,k] + D0[k,j]
          D0[j,i] = D0[i,j]
        }
      }
    }
  }
  max.d = max(D0[!is.infinite(D0)])
  # trimming unconnected cliques
  D0[D0 > max.d] = max.d

  return(D0)
}

estimate_nus <- function(A,clique.set){
  K <- length(clique.set)
  nu.set <- list()
  for(k in 1:K){
    p.hats <- colMeans(A[,as.vector(clique.set[[k]])])
    nu.hats <- log(p.hats) - log(max(p.hats))
    nu.set[[k]] <- nu.hats
  }
  return(nu.set)
}

likelihood_value <- function(G, cliques, fixed.effect.vec, D0){

  K = max(cliques)
  nu.big <- outer(fixed.effect.vec,fixed.effect.vec, "+")
  obj.val <- 0
  K.large <- dim(G)[1]
  D.big <- matrix(0,K.large,K.large)

  #full.idx <- seq(length(cliques))
  for(k1 in seq(K)){
    for(k2 in seq(K)){
      idx1 = which(cliques == k1)
      idx2 = which(cliques == k2)
      D.big[idx1,idx2] = D0[k1,k2]
      # nu.block <- nu.big[idx1,idx2]
      # A.block <- G[idx1,idx2]
    }
  }


  lik.mat = G*(nu.big - D.big) + (1 - G)*log(1 - exp(nu.big - D.big))
  lik.mat[is.infinite(lik.mat)] = 0 # not counting self-clique connections
  lik.mat[is.na(lik.mat)] = 0 # not counting self-clique connections
  obj.val <- sum(lik.mat, na.rm = T) # removing the erroneous terms
  return(obj.val)
}

max_triangle_deviation <- function(D){
  K <- dim(D)[1]

  grid = expand.grid(seq(K),seq(K),seq(K))
  idx1 = grid[,1] < grid[,2]
  idx2 = grid[,2] < grid[,3]
  grid <- grid[idx1 & idx2, ]
  dev.vec <- sapply(seq(nrow(grid)),function(i){
    j = grid[i,1]
    k = grid[i,2]
    l = grid[i,3]
    dev1 <- -(D[j,k] + D[k,l] - D[j,l])
    dev2 <- -(D[j,k] - D[k,l] + D[j,l])
    dev3 <- -(-D[j,k] + D[k,l] + D[j,l])
    dev.max <- max(c(dev1,dev2,dev3))
    return(dev.max)
  })
  return(max(dev.vec))
}

clique_partition <- function(clique.set, randomize = F){
  K <- length(clique.set)
  idx <- 1:K
  if(randomize){
    idx <- sample(idx,K)
  }

  included.cliques <- c()
  included.idx <- c()
  for(k in idx){
    clq.tmp <- clique.set[[k]]
    if(!(any(clq.tmp %in% included.cliques))){
      included.cliques <- c(included.cliques,clq.tmp)
      included.idx <- c(included.idx, k)
    }
  }
  reduced.cliques <- list()
  k.sub <- 1
  for(k in included.idx){
    reduced.cliques[[k.sub]] <- clique.set[[k]]
    k.sub <- k.sub + 1
  }

  return(reduced.cliques)
}




rowSDs <- function(df){
  diff.df <- df - rowMeans(df, na.rm = T)
  var.vec <- rowMeans((diff.df)^2)
  sd.vec <- sqrt(var.vec)
  return(sd.vec)
}


rowMedians <- function(df){
  medians <- apply(df, 1, function(z){
    return(median(z, na.rm = T))
  })
  return(medians)
}

perm_median_test <- function(y, x, n.perm = 1000){
  fit = fit_medians(y, x)
  theta = fit$theta
  loss = fit$loss

  loss.dist = rep(0,n.perm)
  n = length(x)
  for(i in seq(n.perm)){
    perm <- sample(seq(n), replace = F)
    x.perm <- x[perm]
    fit.perm = fit_medians(y, x.perm)
    loss.dist[i] = fit.perm$loss
  }
  p.val <- mean(loss.dist <= loss)
  out.list <- list("theta" = theta, "loss" = loss, "p-value" = p.val)
  return(out.list)
}

fit_medians <- function(y, x){
  J = max(x)
  theta = rep(NA,J)
  ell <- 0
  for(j in 1:J){
    idx = which(x == j)
    theta[j] = median(y[idx])
    ell <- ell + sum(abs(y[idx] - theta[j]))
  }
  out.list <- list("theta" = theta, "loss" = ell)
  return(out.list)
}


# wrapped hyperbolic normal
rwhn <- function(n,mu,sigma = diag(1,length(mu) - 1)){
  p <- length(mu) - 1

  mu0 <- c(1,rep(0,p))
  nu.set <- Rfast::rmvnorm(n, rep(0,p), sigma)
  nu.set <- cbind(rep(0,n),nu.set)

  z <- apply(nu.set,1, function(nu){
    u.tmp <- pt(nu,mu0,mu)
    z.tmp <- exp_map(u.tmp,mu)

    return(z.tmp)
  })

  z <- t(z)
  return(z)
}



pt <- function(v,mu0,mu){
  alpha <- -inner_l(mu0,mu)
  out <- v + (1/(alpha + 1))*inner_l(mu - alpha*mu0,v)*(mu0 + mu)
  return(out)
}

exp_map <- function(u,mu){
  norm.u <- norm_l(u)
  out <- cosh(norm.u)*mu + sinh(norm.u)*(u/norm.u)
  return(out)
}

inner_l <- function(z1,z2){
  p = length(z1)
  # if(z1[1] <= 0 | z2[1] <= 0 ){
  #   stop("first index must be positive")
  # }
  I = diag(1,p)
  I[1,1] = -1
  ip = as.numeric(t(z1) %*% I %*% z2)
  return(ip)
}

norm_l <- function(z){
  return(sqrt(inner_l(z,z)))
}


# warped manifold versions of the normal distributions.
# All are isotropic



rmanifn <- function(n,mu,var.scale,kappa){
  if(kappa == 0){
    Z <- Rfast::rmvnorm(n,mu,diag(var.scale, length(mu)))
  } else if(kappa > 0){
    Z <- Rfast::rvmf(n,mu,1/max(var.scale,10**(-7)))
  } else if(kappa < 0){
    Z <- rwhn(n,mu,diag(var.scale, length(mu) - 1))
  }
  return(Z)
}



# kappa = 0
# Z <- rmanifn(1000,c(0,0), var.scale,kappa)
# D <- pos_to_dist(Z, kappa)
# mean(D)
#
# kappa = 1
# Z <- rmanifn(1000,c(1,0,0), var.scale,kappa)
# D <- pos_to_dist(Z, kappa)
# mean(D)
#
# kappa = -1
# Z <- rmanifn(1000,c(1,0,0), var.scale,kappa)
# D <- pos_to_dist(Z, kappa)
# mean(D)

latent_position_cluster_model <- function(n,n.centers, p, centers.radius, kappa,
                                          variance.scales = rep(0.05,n.centers),
                                          PI = rep(1/n.centers,n.centers), flatness = 1){

  cluster.sizes <- as.numeric(rmultinom(n = 1, size = n, prob = PI))
  # TODO: Replace this back to the old version if it doesnt work
  if(kappa > 3){
    centers <- sim_latent_uniform_ball(n.centers,p,kappa,centers.radius, flatness)
  } else {
    # n.inner <- round(n.centers/2)
    # n.outer <- n.centers - n.inner
    # centers1 <- sim_projected_uniform_ball(n.outer,p,kappa,centers.radius, flatness)
    # # inner radius
    # inner.radius <- centers.radius/(2.5)
    # #centers2 <- sim_projected_uniform_ball(n.inner,p,kappa,inner.radius)
    # centers2 <- sim_projected_conic_distribution(n.centers,p,kappa,inner.radius)
    # centers <- rbind(centers1,centers2)

    ## Alternative:
    centers <- sim_projected_conic_distribution(n.centers,p,kappa,centers.radius, flatness)

    }



  clust.labels <- c()
  Z <- NULL
  for(clust.idx in seq(n.centers)){
    clust.labels <- c(clust.labels, rep(clust.idx,cluster.sizes[clust.idx]))
    cent <- centers[clust.idx,]
    scale <- variance.scales[clust.idx]
    Z.block <- rmanifn(cluster.sizes[clust.idx],cent,scale,kappa)

    #print(paste0(cluster.sizes[clust.idx]," ::: ", nrow(Z.block)))
    if(is.null(Z)){
      Z <- Z.block
    } else {

      if(dim(Z.block)[2] == dim(Z)[2] ){
        Z <- rbind(Z,Z.block)
      }

    }

  }
  out.list <- list("Z" = Z, "cluster_labels" = clust.labels)
  return(out.list)
}


latent_position_cluster_model_2 <- function(n,n.centers, p, kappa,
                                            centers.variance = 0.5**2,
                                            cluster.variance = 0.25**2,
                                            PI = rep(1/n.centers,n.centers)){

  cluster.variance.vec = rgamma(n.centers, shape = cluster.variance)
  cluster.sizes <- as.numeric(rmultinom(n = 1, size = n, prob = PI))


  if(kappa != 0 ){
    ref.center <- c(1,rep(0,p))
    centers <- rmanifn(n.centers,ref.center,centers.variance,kappa)
  } else {
    ref.center <- rep(0,p)
    centers <- rmanifn(n.centers,ref.center,centers.variance,kappa)
  }


  clust.labels <- c()
  Z <- NULL
  for(clust.idx in seq(n.centers)){
    clust.labels <- c(clust.labels, rep(clust.idx,cluster.sizes[clust.idx]))
    cent <- centers[clust.idx,]
    scale <- cluster.variance.vec[clust.idx]
    Z.block <- rmanifn(cluster.sizes[clust.idx],cent,scale,kappa)

    #print(paste0(cluster.sizes[clust.idx]," ::: ", nrow(Z.block)))
    if(is.null(Z)){
      Z <- Z.block
    } else {

      if(dim(Z.block)[2] == dim(Z)[2] ){
        Z <- rbind(Z,Z.block)
      }

    }

  }
  out.list <- list("Z" = Z, "cluster_labels" = clust.labels)
  return(out.list)
}






prod_space_dist <- function(coord.array,curvatures){
  D <- matrix(0,nrow = dim(coord.array)[1], ncol = dim(coord.array)[1])
  for(k in seq(length(curvatures))){
    kappa <- curvatures[k]
    Z <- coord.array[,,k]
    na.cols <- is.na(Z[1,])
    Z <- as.matrix(Z[,!na.cols])
    D.mat <- pos_to_dist(Z,kappa)
    D <- D +  D.mat^2
  }
  D <- sqrt(D)
  diag(D) = 0
  return(D)
}


normal_point_mass_test <- function(x,sigma, B = 1000){
  mu = mean(x)
  n = length(x)
  loss <- mean( (x - mu)^2)
  loss.vec <- rep(0,B)
  for(b in seq(B)){
    x.sim <- rnorm(n,mu,sigma)
    mu.sum <- mean(x.sim)
    loss.sim <- mean( (x.sim - mu.sum)^2)
    loss.vec[b] <- loss.sim
  }
  p.value <- mean(loss.vec > loss)
  return(p.value)
}



normal_point_mass_test <- function(x,sigma, B = 1000){
  mu = mean(x)
  n = length(x)
  loss <- mean( (x - mu)^2)
  loss.vec <- rep(0,B)
  for(b in seq(B)){
    x.sim <- rnorm(n,mu,sigma)
    mu.sum <- mean(x.sim)
    loss.sim <- mean( (x.sim - mu.sum)^2)
    loss.vec[b] <- loss.sim
  }
  p.value <- mean(loss.vec > loss)
  return(p.value)
}

rlaplace <- function(n,mu = 0,rate = 1){
  U <- runif(n)
  Y <- rexp(n,rate)
  X <- (U < 1/2)*Y + (U >= 1/2)*(-Y) + mu
  return(X)
}

laplace_point_mass_test <- function(x,rate, B = 1000){
  mu = median(x)
  n = length(x)
  loss <- mean( abs(x - mu))
  loss.vec <- rep(0,B)
  for(b in seq(B)){
    x.sim <- rlaplace(n,mu,rate)
    mu.sum <- median(x.sim)
    loss.sim <- mean( abs(x.sim - mu.sum))
    loss.vec[b] <- loss.sim
  }
  p.value <- mean(loss.vec > loss)
  return(p.value)
}

cdf <- function(x){
  out <- rep(0,length(x))
  for(j in x){
    out[j] <- mean(x <= x[j])
  }
  return(out)
}

surv <- function(x){
  out <- rep(0,length(x))
  for(j in x){
    out[j] <- mean(x >= x[j])
  }
  return(out)
}


connected_spheres_uniform_sim <- function(n1, n2, p1,p2,kappa1 = 1,kappa2 = 1.5){
  R1 <- pi/sqrt(kappa1)
  R2 <- pi/sqrt(kappa2)
  Z1 <- sim_latent_uniform_ball(n1,p1,kappa1,R1)
  Z2 <- sim_latent_uniform_ball(n2,p2,kappa2,R2)
  D1 <- pos_to_dist(Z1,kappa1)
  D2 <- pos_to_dist(Z2,kappa2)
  #O1 <- sim_latent_uniform_ball(n1,p1,kappa1,0)
  # O1 <- t(Z1)
  # O1[,] <-
  # O2 <- t(Z2)
  # O2[,] <- c(1,rep(0,p2))
  # O1 <- t(O1)
  # O2 <- t(O2)
  O1 <- matrix(c(1,rep(0,p1)), nrow = 1)
  O2 <- matrix(c(1,rep(0,p2)), nrow = 1)

  do1 <- as.numeric(pos_to_dist_pair(Z1, O1, kappa1))
  do2 <- as.numeric(pos_to_dist_pair(Z2, O2, kappa2))

  D12 <- outer(do1, do2, "+")

  D <- matrix(0, nrow = n1+n2, ncol = n1+n2)
  D[1:n1,1:n1] <- D1
  D[(n1 + (1:n2)),(n1 + (1:n2))] <- D2
  D[1:n1,(n1 + (1:n2))] <- D12
  diag(D) = 0
  shuffle <- sample(1:(n1 + n2))
  D <- D[shuffle,shuffle]
  return(D)
}


connected_spheres_lpcm <- function(n, n.centers1, n.centers2, p1,p2,PI1, PI2, rand.eff = rep(0,n), q = 0.5, kappa1 = 1,kappa2 = 1.5, approximate.variance = 0.25**2, max.rad = 2.5, shuffle.out = F, sim.A = T, compute.D = F){
  n1 <- rbinom(1,n,prob = q)
  n2 <- n - n1

  R1 <- min(pi/sqrt(kappa1) - 10**(-6), max.rad)
  R2 <- min(pi/sqrt(kappa2) - 10**(-6), max.rad)


  variance.scales1 <- rgamma(n.centers1, shape = approximate.variance)
  variance.scales2 <- rgamma(n.centers2, shape = approximate.variance)
  lpcm1 <- latent_position_cluster_model(n1,n.centers1, p1, R1, kappa1, variance.scales = variance.scales1, PI = PI1)
  lpcm2 <- latent_position_cluster_model(n2,n.centers2, p2, R2, kappa2, variance.scales = variance.scales2, PI = PI2)


  Z1 <- lpcm1$Z
  Z2 <- lpcm2$Z

  rand.eff1 <- rand.eff[1:n1]
  rand.eff2 <- rand.eff[(n1 + 1):(n)]

  if(sim.A){
    A1 <- sim_ls_network_fast_2(rand.eff1, Z1, kappa1)
    A2 <- sim_ls_network_fast_2(rand.eff2, Z2, kappa2)

    O1 <- matrix(c(-1,rep(0,p1)), nrow = 1)
    O2 <- matrix(c(-1,rep(0,p2)), nrow = 1)

    do1 <- as.numeric(pos_to_dist_pair(Z1, O1, kappa1))
    do2 <- as.numeric(pos_to_dist_pair(Z2, O2, kappa2))

    D12 <- outer(do1, do2, "+")
    P12 <- exp(outer(rand.eff1,rand.eff2,"+") - D12)

    U <- runif(length(P12))
    input <- 1*(U <= P12)
    A12 <- matrix(input, n1,n2)
    A <- matrix(0,n,n)
    A[1:n1,1:n1] <- as.vector(A1)
    A[(n1 + 1):(n),(n1 + 1):(n)] = as.vector(A2)
    A[1:n1,(n1 + 1):(n)] = as.vector(A12)
    A[(n1 + 1):(n),1:n1] = as.vector(t(A12))
  } else{
    A <- NULL
  }
  if(compute.D){
    D1 <- pos_to_dist_pair(Z1, Z1, kappa1)
    D2 <- pos_to_dist_pair(Z2, Z2, kappa2)
    #O1 <- sim_latent_uniform_ball(n1,p1,kappa1,0)
    # O1 <- t(Z1)
    # O1[,] <-
    # O2 <- t(Z2)
    # O2[,] <- c(1,rep(0,p2))
    # O1 <- t(O1)
    # O2 <- t(O2)

    # minus one for less connection through the manifold
    O1 <- matrix(c(-1,rep(0,p1)), nrow = 1)
    O2 <- matrix(c(-1,rep(0,p2)), nrow = 1)

    do1 <- as.numeric(pos_to_dist_pair(Z1, O1, kappa1))
    do2 <- as.numeric(pos_to_dist_pair(Z2, O2, kappa2))

    D12 <- outer(do1, do2, "+")

    D <- matrix(0, nrow = n1+n2, ncol = n1+n2)
    D[1:n1,1:n1] <- D1
    D[(n1 + (1:n2)),(n1 + (1:n2))] <- D2
    D[1:n1,(n1 + (1:n2))] <- D12
    D[(n1 + (1:n2)),1:n1] <- t(D12)
  } else {
    D <- NULL
  }


  if(shuffle.out){
    shuffle <- sample(1:(n1 + n2))
    D <- D[shuffle,shuffle]
    clust.labels <- c(lpcm1$cluster_labels, max(lpcm1$cluster_labels) + lpcm2$cluster_labels)
    clust.labels <- clust.labels[shuffle]
    A <- A[shuffle,shuffle]
    out.list <- list("D" = D, "cluster_labels" = clust.labels, "A" = A)
  } else{
    clust.labels <- c(lpcm1$cluster_labels, max(lpcm1$cluster_labels) + lpcm2$cluster_labels)
    out.list <- list("D" = D, "cluster_labels" = clust.labels, "A" = A)
  }

  return(out.list)
}



clique_split <- function(clique.set, min_clique_size){
  K <- length(clique.set)
  new.clique.set <- list()
  for(k in seq(K)){
    ellk <- length(clique.set[[k]])
    if(ellk >= 2*min_clique_size){
      list1 <- list(clique.set[[k]][1:min_clique_size])
      list2 <- list(clique.set[[k]][min_clique_size:ellk])

      new.clique.set <- append(new.clique.set, list1)
      new.clique.set <- append(new.clique.set, list2)
    } else {
      list1 <- list(clique.set[[k]])
      new.clique.set <- append(new.clique.set, list1)
    }
  }
  return(new.clique.set)
}


# cliques are known in this case
# clique searching is too much work
# or use the spectral clustering heuristic
lpcm_spherical_rand_walk_dist <- function(n,n.centers, p, kappa = 1,
                                          approximate.variance = 0.05,
                                          PI = rep(1/n.centers,n.centers),
                                          centers.radius = 2.5, Time.steps = 50, rho = 1/2){
  cluster.model.variance = rgamma(n.centers, shape = approximate.variance)

  #centers.radius = pi/sqrt(kappa) - 10**(-6)
  centers <- sim_latent_uniform_ball(n.centers,p,kappa,centers.radius)
  cluster.sizes <- as.numeric(rmultinom(n = 1, size = n, prob = PI))
  #mean(cluster.model.variance)

  Z.set <- list()
  Z.prev = NULL
  for(time in seq(Time.steps)){
    clust.labels <- c()
    Z <- NULL
    for(clust.idx in seq(n.centers)){
      clust.labels <- c(clust.labels, rep(clust.idx,cluster.sizes[clust.idx]))
      cent <- centers[clust.idx,]
      var.scale <- cluster.model.variance[clust.idx]
      var.scale = max(10**(-8), var.scale)
      Z.block <- rmanifn(cluster.sizes[clust.idx],cent,var.scale,kappa)
      #print(paste0(cluster.sizes[clust.idx]," ::: ", nrow(Z.block)))
      if(is.null(Z)){
        Z <- Z.block
      } else {
        Z <- rbind(Z,Z.block)
      }

    }

    if(is.null(Z.prev)){
      Z.prev = Z
    } else {
      Z = (1 - rho)*Z.prev + rho*(Z)
      Z = Z/sqrt(rowSums(Z^2))
      Z.prev = Z
    }
    Z.set[[time]] <- Z


    cat(paste("Time:", time, "/", Time.steps), end = "\r")
  }

  out.list <- list("Z.set" = Z.set, "cluster_labels" = clust.labels)
  rm(Z.set)
  return(out.list)
}




# search for cliques withing a clustered set.
# this can be used with any clustering,
# in our case, for the time series case, we just assume these are known.
guided_clique_set <- function(G,labels, min_clique_size = 8){
  K <- max(labels)
  clique.set = list()
  for(k in seq(K)){
    idx.k <- which(labels == k)
    Gk <- G[idx.k, idx.k]
    gk <- igraph::graph_from_adjacency_matrix(Gk, mode = "undirected")
    clique <- igraph::largest.cliques(gk)
    if(length(clique) > 0 ){
      if(length(clique[[1]]) > min_clique_size){
        new_clique <- list(idx.k[clique[[1]]])
        clique.set <- append(clique.set, new_clique)
      }
    }
  }
  return(clique.set)
}


#TODO: Estimate the curvature here,

# average_noise_estimate <- function(D, nu.vec){
#
#   # use the base D as the initial estimate
#   # clique locations are preserved
#
#
# }
#
#G = A.sim

estimate_curvature <- function(G, clique.set,
                               c1 = 0.5,
                               c2 = 2,
                               c3 = 0.25,
                               d.yz.min = 1.5,
                               d.yz.max,
                               verbose = T,
                               tri.const = 1.4,
                               rand.eff.0 = F,
                               max.iter = 50,
                               no.refit = F,
                               D0){
  nu.hats <- estimate_nus(G, clique.set)
  K <- length(clique.set)
  ell <- Inf
  clique.idx <- c()
  subset.ids <- c()
  fixed.effect.vec <- c()
  for(k in seq(K)){
    clique.idx <- c(clique.idx,clique.set[[k]])
    subset.ids <- c(subset.ids, rep(k,length(clique.set[[k]])))
    fixed.effect.vec <- c(fixed.effect.vec,nu.hats[[k]])
    if(ell > length(clique.set[[k]])){
      ell = length(clique.set[[k]])
    }
  }
  if(rand.eff.0){
    fixed.effect.vec[] = 0
  }
  G.sub <- G[clique.idx,clique.idx]
  diag(G.sub) = 0
  if(missing(d.yz.max)){
    d.yz.max = max(log(ell), log(ell^2/10))
  }

  if(missing(D0)){
    D0 = init_D0(G.sub,subset.ids,fixed.effect.vec)
  }
  # this seems to be the fastest LCQP that I can use.
  #D0 <- 0.5*D0
  if(no.refit){
    D.hat <- D0
  } else {
    D.hat <- estimate_D_restricted_fast(G.sub,subset.ids,
                                        fixed.effect.vec,
                                        D0,thresh = 10**(-3),
                                        max.iter = max.iter,
                                        solver = "MOSEK",
                                        verbose = T)
  }


  l1 <- likelihood_value(G.sub,subset.ids,fixed.effect.vec,D0)

  l2 <- likelihood_value(G.sub,subset.ids,fixed.effect.vec,D.hat)
  if(l1 >= l2){
    warning("Likelihood did not increase, returning initial estimate")
    D.hat <- D0
  }
  mid.search <- optimal_midpoint_search(D.hat,top.k = 10,
                                        d.yz.min = min(d.yz.min, max(D.hat)/2),
                                        d.yz.max = d.yz.max)
  if(verbose){
    print("Midpoints: ")
    max.row <- min(3,nrow(mid.search))
    print(mid.search[1:max.row,])
  }


  y.opt = mid.search[1,1]
  z.opt = mid.search[1,2]
  m.opt = mid.search[1,3]



  opt.vec <- c(y.opt,z.opt,m.opt)
  if(!any(is.na(opt.vec))){
    x.set <- filter_indices_2(D.hat,
                              y.opt,
                              z.opt,
                              m.opt,
                              tri.const = tri.const)

    kappa.set <- estimate_kappa_set(D.hat,
                                    y.opt,
                                    z.opt,
                                    m.opt,
                                    x.set)

    out.set <- list("kappas" = kappa.set,
                    "kappa.med" = median(kappa.set, na.rm = T),
                    "D" = D.hat,
                    "midpoints" = mid.search)

  } else {
    out.set <- list("kappas" = NULL,
                    "kappa.med" = NULL,
                    "D" = D.hat,
                    "midpoints" = mid.search)
  }

  return(out.set)
}


plot_cliques <- function(G,clique.set){

  labels <- c()
  clique.ids <- c()
  for(lab in 1:length(clique.set)){
    labels <- c(labels, rep(lab,length(clique.set[[lab]])))
    clique.ids <- c(clique.ids,clique.set[[lab]])
  }

  A.sub <- G[clique.ids,clique.ids]


  #D.hat <- estimate_D(A,clique.set,nu.hats)

  g <- graph_from_adjacency_matrix(A.sub,mode = "undirected")


  # doesn't seem to work well
  # c28

  n <- length(clique.set)
  if(n <= 28){
    col_vector <- c("dodgerblue2", "#E31A1C", # red
                    "green4",
                    "#6A3D9A", # purple
                    "#FF7F00", # orange
                    "brown", "gold1",
                    "skyblue2", "#FB9A99", # lt pink
                    "palegreen2",
                    "#CAB2D6", # lt purple
                    "#FDBF6F", # lt orange
                    "gray70", "khaki2",
                    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                    "darkturquoise", "green1", "yellow4", "yellow3",
                    "darkorange4", "red", "white", "darkgrey"
                  )
  } else {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  }


  V(g)$color <- col_vector[labels]
  plt <- plot(g,vertex.size= 6,vertex.label=NA)
  return(plt)
}



plot_cliques_curvature <- function(G,D,mid.set,clique.set, num.midpoints = 3, tri.const = 1.4){

  #D <- estimates$D
  #mid.set <- estimates$midpoints
  K = nrow(D)


  m.set <- c()
  y.set <- c()
  z.set <- c()
  kappas <- c()
  for(k in seq(num.midpoints)){
    y.opt = mid.set[k,1]
    z.opt = mid.set[k,2]
    m.opt = mid.set[k,3]
    m.set[k] <- m.opt
    y.set[k] <- y.opt
    z.set[k] <- z.opt
    # x.set <- filter_indices(D.hat, y.opt,
    #                         z.opt,m.opt,
    #                         c1 = c1,c2 = c2,c3 = c3)
    x.set<- lolaR::filter_indices(D, y.opt,
                                  z.opt, m.opt,
                                  tri.const = tri.const)
      #filter_indices_2(D, y.opt,
            #                 z.opt,m.opt,
            #                 tri.const = tri.const)

    kappa.set <- lolaR::EstimateKappaSet(D,y.opt,z.opt,m.opt,x.set)
    kappas[k] <- median(kappa.set, na.rm = T)

  }



  kappas <- round(kappas,4)

  # maximum for plotting
  max.curve = 10
  kappas[kappas < -max.curve] = -max.curve
  kappas[kappas > max.curve] = max.curve

  labels <- c()
  clique.ids <- c()
  for(lab in 1:length(clique.set)){
    labels <- c(labels, rep(lab,length(clique.set[[lab]])))
    clique.ids <- c(clique.ids,clique.set[[lab]])
  }

  closest.curves <- c()

  for(j in seq(K)){
    d.m.mini <- D[j,m.set]
    d.y.mini <- D[j,y.set]
    d.z.mini <- D[j,z.set]
    d.mini <- c()
    for(k in seq(num.midpoints)){
      d.mini[k] <- min(d.m.mini[k],d.y.mini[k],d.z.mini[k])
    }
    j.hat <- which.min(d.mini)
    closest.curves[j] = kappas[j.hat]
  }



  R <- sum(abs(unique(kappas)) > 0)

  col.pal <- diverge_hsv(2*R + 1) #
  kappas.mirror <- unique(c(kappas,-kappas,0))
  mirror.kappas <- sort(kappas.mirror)

  curve.color.idx <- c()

  for(j in seq(K)){
    curve.color.idx[j] = which(mirror.kappas == closest.curves[j])
  }

  label.set <- rep(NA, K)
  for(k in seq(num.midpoints)){
    label.set[y.set[k]] = k
    label.set[z.set[k]] = k
    label.set[m.set[k]] = k
  }



  A.sub <- G[clique.ids,clique.ids]
  col.short <- col.pal[curve.color.idx] # cliques indicator
  label.short <- label.set
  col.long <- c()
  label.long <- c()
  for(lab in 1:length(clique.set)){
    col.long <- c(col.long, rep(col.short[lab],length(clique.set[[lab]])))
    #label.long <- c(label.long,rep(label.short[lab],length(clique.set[[lab]])))
    label.long <- c(label.long,label.short[lab],rep(NA,length(clique.set[[lab]]) - 1))
  }

  #plot(x = 1:length(col.short), y = col.short)
  #plot(x=1:length(col.pal), y=rep(0,length(col.pal)), pch=19, cex=15, col=col.pal)
  #D.hat <- estimate_D(A,clique.set,nu.hats)
  #show_palette(col.pal)
  #display.brewer.pal(col.pal)
  g <- graph_from_adjacency_matrix(A.sub,mode = "undirected")

  #label.long
  V(g)$color <- col.long

  shape.long <- ifelse(is.na(label.long), "circle", "square")

  #V(g)$label <- label.long
  V(g)$label <- NA
  V(g)$shape <- shape.long

  V(g)$label.cex <- 2
  V(g)$vertex.label.color <- "black"
  #plt <- plot(g,vertex.size= 6,vertex.label=NA)
  plt <- plot(g,vertex.size= 6, vertex.label.color = "gold1")

  return(plt)
}




# heavy_tail_transformation <- function(x){
#   out <- x*(x >= 0) + -log(1 + abs(x))*(x < 0 & !is.infinite(x))
#   return(out)
# }


# Lubold_Estimate_Curvature_Sphere <- function(D,a = 0.1, b = 10){
#   b.max <- min(b,max((pi/max(D))^2))
#   kappaVec = seq(a, b.max, length.out =10000)
#   K <- nrow(D)
#   index <- K
#   objFun = rep(0,(length(kappaVec)))
#   for(i in seq(length(kappaVec))){
#     C = cos(D * sqrt(kappaVec[i]))
#     eig <- eigen(C, symmetric = T)
#     objFun[i] = abs(eig$values[index])
#   }
#   index_min = which.min(objFun)
#   #plot(objFun)
#   return(kappaVec[index_min])
# }


# Lubold_Estimate_Curvature_Hyperbolic <- function(D,a = -10, b = -0.1){
#   kappaVec = seq(a, b, length.out =10000)
#   K <- nrow(D)
#   index <- 2
#
#   objFun = rep(0,(length(kappaVec)))
#   for(i in seq(length(kappaVec))){
#     C = cosh(D * sqrt(-kappaVec[i]))
#     eig <- eigen(C, symmetric = T)
#     objFun[i] = abs(eig$values[index])
#   }
#
#   index_min = which.min(objFun)
#   #plot(objFun)
#   return(kappaVec[index_min])
# }


show_palette <- function(colors) {
  n = length(colors)
  image(1:n, 1, as.matrix(1:n), col = colors,
        xlab = "", ylab = "", xaxt = "n",
        yaxt = "n", bty = "n")
}





