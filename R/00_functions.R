


library(igraph)
library(dplyr)

library(ggplot2)
library(CVXR)
#
library(truncnorm)
library(Matrix)

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



# sampling from a uniform ball in various geometries
sim_latent_uniform_ball <- function(n,p,kappa,radius){
  # equal angle sampler
  X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = diag(rep(1,p)))
  X <- cbind(rep(0,n),X)
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
  A <- A + t(A)
  return(A)
}


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



colSDs <- function(df, ...){
  diff.df <- t(df) - colMeans(df, ...)
  diff.df <- t(diff.df)
  var.vec <- colMeans((diff.df)^2, ...)
  sd.vec <- sqrt(var.vec)
  return(sd.vec)
}



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




scale_curvature <- function(x,c = 1){
  out <- c*(exp(2*x/c) - 1)/(exp(2*x/c) + 1)
  out[is.infinite(x) & x > 0] = c
  return(out)
}




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


latent_position_cluster_model <- function(n,n.centers, p, centers.radius, kappa,
                                          variance.scales = rep(0.05,n.centers),
                                          PI = rep(1/n.centers,n.centers), flatness = 1,
                                          radius.ratio = 1/centers.radius, sample.ratio = 0.5, force.midpoints = F, midpoint.ratio  = .2){

  cluster.sizes <- as.numeric(rmultinom(n = 1, size = n, prob = PI))

  if(force.midpoints){

  } else {
    if(kappa > 0){
      centers <- sim_latent_uniform_ball(n.centers,p,kappa,centers.radius, flatness)
    } else {
      n.inner <- round(n.centers*sample.ratio)
      n.outer <- n.centers - n.inner
      centers1 <- sim_latent_uniform_ball(n.outer,p,kappa,centers.radius, flatness)
      # inner radius
      inner.radius <- centers.radius*radius.ratio
      centers2 <- sim_latent_uniform_ball(n.inner,p,kappa,inner.radius, flatness)
      #centers2 <- sim_projected_conic_distribution(n.centers,p,kappa,inner.radius)
      centers <- rbind(centers1,centers2)

      ## Alternative:
      #centers <- sim_projected_conic_distribution(n.centers,p,kappa,centers.radius, flatness)

    }
  }



  clust.labels <- c()
  Z <- NULL
  for(clust.idx in seq(n.centers)){
    clust.labels <- c(clust.labels, rep(clust.idx,cluster.sizes[clust.idx]))
    cent <- centers[clust.idx,]
    scale <- variance.scales[clust.idx]
    Z.block <- rmanifn(cluster.sizes[clust.idx],cent,scale,kappa)

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


show_palette <- function(colors) {
  n = length(colors)
  image(1:n, 1, as.matrix(1:n), col = colors,
        xlab = "", ylab = "", xaxt = "n",
        yaxt = "n", bty = "n")
}





