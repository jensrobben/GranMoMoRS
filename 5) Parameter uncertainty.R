##### 0) Read global environment ----
# Clear environment
rm(list = ls())
gc()

# Packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble',
              'Polychrome','tidyr','ISOweek', 'ggplot2','RColorBrewer',
              'data.table', 'viridis', 'imputeTS', 'imputeTS',
              'parallel','doParallel', 'spdep', 'tidyr','abind',
              'mvtnorm','randtoolbox','MASS','snow','doSNOW',
              'progress')
invisible(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load
load("Results/EMresult_tau10.RDATA")


##### 1) Parameter estimates RS model -----

# Object
obj[[length(obj)]] <- NULL
iteropt <- length(obj)

# Iteration results
alpha1.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.1))
alpha2.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.2))
beta0.iter  <- do.call('rbind', lapply(obj, function(x) x$beta.0))
beta11.iter <- do.call('rbind', lapply(obj, function(x) x$beta.11))
beta22.iter <- do.call('rbind', lapply(obj, function(x) x$beta.22))
gamma.iter  <- do.call('rbind', lapply(obj, function(x) x$gamma))
U.iter      <- do.call('rbind', lapply(obj, function(x) x$U)) 

theta.iter <- cbind(alpha1.iter, alpha2.iter, beta0.iter, 
                    beta11.iter, beta22.iter, gamma.iter)

# Length vector
d <- ncol(theta.iter) - 2

# Optimal parameter vector
alpha1.opt <- alpha1.iter[iteropt,]
alpha2.opt <- alpha2.iter[iteropt,]
beta0.opt  <- beta0.iter[iteropt,]
beta01.opt <- beta0.opt[1:ncol(xx01)]
beta02.opt <- beta0.opt[-c(1:ncol(xx01))]
beta11.opt <- beta11.iter[iteropt,]
beta22.opt <- beta22.iter[iteropt,]
gamma.opt  <- gamma.iter[iteropt,1] 
U.opt      <- U.iter[iteropt,]

alpha.opt  <- c(alpha1.opt, alpha2.opt)
beta.opt   <- c(beta0.opt, beta11.opt, beta22.opt)
theta.opt  <- unname(c(alpha1.opt, alpha2.opt, beta0.opt, beta11.opt,
                       beta22.opt))

##### 2) Gradient functions of the log-likelihood ----

# Gradient functions
S.alpha   <- function(alpha, b01, b02, b1, b2, dtsim){
  # Length vectors
  ln   <- sapply(list(alpha1.opt, alpha2.opt), function(v) length(v))
  cln  <- c(1,cumsum(ln))
  
  # Parameter values
  a1  <- alpha[c(cln[1]):cln[2]]
  a2  <- alpha[c(cln[2]+1):cln[3]]
  b0  <- c(b01, b02)
  b11 <- b1
  b22 <- b2
  
  # Conditional Poisson densities
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  fdts1 <- matrix(dpois(dtsim, lambda = bt), nrow = nT*nA, ncol = R)
  fdts2 <- matrix(dpois(dtsim, lambda = bt * exp(x1a1)), nrow = nT*nA, ncol = R)
  fdts3 <- matrix(dpois(dtsim, lambda = bt * exp(x2a2)), nrow = nT*nA, ncol = R)
  fdts  <- cbind(fdts1, fdts2, fdts3)
  
  # Change 0 values by very small positive value to avoid numerical issues
  fdts[which(fdts == 0)] <- 10^(-30)
  
  x1b01 <- pen01 + sapply(1:R, function(r) xx01[,,r] %*% b01 + U.opt[r])
  x2b02 <- pen02 + sapply(1:R, function(r) xx02[,,r] %*% b02 + U.opt[r]) 
  x1b11 <- pen11 + sapply(1:R, function(r) xx11[,,r] %*% b11 + U.opt[r])
  x2b22 <- pen22 + sapply(1:R, function(r) xx22[,,r] %*% b22 + U.opt[r]) 
  
  # All info
  exp.x1b01 <- exp(x1b01) 
  exp.x2b02 <- exp(x2b02)
  exp.x1b11 <- exp(x1b11) 
  exp.x2b22 <- exp(x2b22) 
  
  p00 <- 1/(1 + exp.x1b01 + exp.x2b02)
  p01 <- exp.x1b01/(1 + exp.x1b01 + exp.x2b02)
  p02 <- exp.x2b02/(1 + exp.x1b01 + exp.x2b02)
  
  p10 <- 1/(1 + exp.x1b11)
  p11 <- exp.x1b11/(1 + exp.x1b11)
  
  p20 <- 1/(1 + exp.x2b22)
  p22 <- exp.x2b22/(1 + exp.x2b22)
  
  # avoid numerical issues
  p00[which(p00 == 0)] <- 10^(-300)
  p01[which(p01 == 0)] <- 10^(-300)
  p02[which(p02 == 0)] <- 10^(-300)
  p10[which(p10 == 0)] <- 10^(-300)
  p11[which(p11 == 0)] <- 10^(-300)
  p20[which(p20 == 0)] <- 10^(-300)
  p22[which(p22 == 0)] <- 10^(-300)
  
  P <- unname(cbind(p00, p01, p02, p10, p11, p11*0, p20, p22*0, p22))
  P[1,] <- NA
  
  # Calculate joint state probabilities
  lst <- rep(list(list('a' = rep(NA, R*9),
                       'b' = rep(NA, R),
                       'c' = rep(NA, R*9))),nT)
  
  for(t in 2:nT){
    if(t == 2) {
      gR <- rep(gamma.iter[iteropt,], each = R)
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        c(rep(gR[1:R], times = 3), rep(gR[(R+1):(2*R)], times = 3),
          rep(gR[(2*R+1):(3*R)], times = 3))
    } else {
      lst_slice <- lst[[t-1]][['c']]
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        (lst_slice[indl0] + lst_slice[indl1] + lst_slice[indl2])
    }
    
    b <- rowSums(matrix(a, nrow = R))
    b[b == 0] <- 10^(-30)
    c <- a/rep(b, times = 9)
    
    if(any(is.na(c))) stop('na')
    
    lst[[t]][['a']][1:(9*R)] <- a
    lst[[t]][['b']][1:R]     <- b
    lst[[t]][['c']][1:(9*R)] <- c
  }
  
  js.prob <- do.call('rbind', lapply(1:length(lst), function(x) lst[[x]][['c']]))
  ms.prob <- rbind(c(sapply(1:R, function(j) sum(js.prob[2,c(j,R+j,2*R+j)])),
                     sapply(1:R, function(j) sum(js.prob[2,c(3*R+j,4*R+j,5*R+j)])),
                     sapply(1:R, function(j) sum(js.prob[2,c(6*R+j,7*R+j,8*R+j)]))),
                   cbind(sapply(1:R, function(j) rowSums(js.prob[-1,c(j,3*R+j,6*R+j)])),
                         sapply(1:R, function(j) rowSums(js.prob[-1,c(R+j,4*R+j,7*R+j)])),
                         sapply(1:R, function(j) rowSums(js.prob[-1,c(2*R+j,5*R+j,8*R+j)]))))
  
  ###  Derivatives
  id   <- 1 + nT*(0:(nA-1))
  
  # Alpha 1
  msprob1 <- do.call('rbind', replicate(nA, ms.prob[-1,(R+1):(2*R)], 
                                        simplify = FALSE))
  vec1    <- (- bt[-id,] * exp(x1a1[-id,]) + dtsim[-id,]) * msprob1
  res1    <- apply(x1[-id,,], 2, function(x) sum(x * vec1)) %>% unname()
  
  # Alpha 2
  msprob2 <- do.call('rbind', replicate(nA, ms.prob[-1,(2*R+1):(3*R)], 
                                        simplify = FALSE))
  vec2    <- (- bt[-id,] * exp(x2a2[-id,]) + dtsim[-id,]) * msprob2
  res2    <- apply(x2[-id,,], 2, function(x) sum(x * vec2)) %>% unname()
  
  # Derivatives
  c(res1, res2)
}
S.beta   <- function(beta, a1, a2, g, dtsim){
  ln <- sapply(list(beta01.opt, beta02.opt, beta11.opt, beta22.opt), 
               function(v) length(v))
  cln  <- c(1,cumsum(ln))
  
  b01 <- beta[c(cln[1]):cln[2]]
  b02 <- beta[c(cln[2]+1):cln[3]]
  b0  <- c(b01, b02)
  b11 <- beta[c(cln[3]+1):cln[4]]
  b22 <- beta[c(cln[4]+1):cln[5]]
  
  # Conditional Poisson densities
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  fdts1 <- matrix(dpois(dtsim, lambda = bt), nrow = nT*nA, ncol = R)
  fdts2 <- matrix(dpois(dtsim, lambda = bt * exp(x1a1)), nrow = nT*nA, ncol = R)
  fdts3 <- matrix(dpois(dtsim, lambda = bt * exp(x2a2)), nrow = nT*nA, ncol = R)
  fdts  <- cbind(fdts1, fdts2, fdts3)
  
  # Change 0 values by very small positive value to avoid numerical issues
  fdts[which(fdts == 0)] <- 10^(-30)
  
  x1b01 <- pen01 + sapply(1:R, function(r) xx01[,,r] %*% b01 + U.opt[r])
  x2b02 <- pen02 + sapply(1:R, function(r) xx02[,,r] %*% b02 + U.opt[r]) 
  x1b11 <- pen11 + sapply(1:R, function(r) xx11[,,r] %*% b11 + U.opt[r])
  x2b22 <- pen22 + sapply(1:R, function(r) xx22[,,r] %*% b22 + U.opt[r]) 
  
  # All info
  exp.x1b01 <- exp(x1b01) 
  exp.x2b02 <- exp(x2b02)
  exp.x1b11 <- exp(x1b11) 
  exp.x2b22 <- exp(x2b22) 
  
  p00 <- 1/(1 + exp.x1b01 + exp.x2b02)
  p01 <- exp.x1b01/(1 + exp.x1b01 + exp.x2b02)
  p02 <- exp.x2b02/(1 + exp.x1b01 + exp.x2b02)
  
  p10 <- 1/(1 + exp.x1b11)
  p11 <- exp.x1b11/(1 + exp.x1b11)
  
  p20 <- 1/(1 + exp.x2b22)
  p22 <- exp.x2b22/(1 + exp.x2b22)
  
  # avoid numerical issues
  p00[which(p00 == 0)] <- 10^(-300)
  p01[which(p01 == 0)] <- 10^(-300)
  p02[which(p02 == 0)] <- 10^(-300)
  p10[which(p10 == 0)] <- 10^(-300)
  p11[which(p11 == 0)] <- 10^(-300)
  p20[which(p20 == 0)] <- 10^(-300)
  p22[which(p22 == 0)] <- 10^(-300)
  
  P <- unname(cbind(p00, p01, p02, p10, p11, p11*0, p20, p22*0, p22))
  P[1,] <- NA
  
  # Calculate joint state probabilities
  lst <- rep(list(list('a' = rep(NA, R*9),
                       'b' = rep(NA, R),
                       'c' = rep(NA, R*9))),nT)
  
  for(t in 2:nT){
    if(t == 2) {
      gR <- rep(c(g,0, 1-sum(g)-10^(-30)), each = R)
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        c(rep(gR[1:R], times = 3), rep(gR[(R+1):(2*R)], times = 3),
          rep(gR[(2*R+1):(3*R)], times = 3))
    } else {
      lst_slice <- lst[[t-1]][['c']]
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        (lst_slice[indl0] + lst_slice[indl1] + lst_slice[indl2])
    }
    
    b <- rowSums(matrix(a, nrow = R))
    b[b == 0] <- 10^(-30)
    c <- a/rep(b, times = 9)
    
    if(any(is.na(c))) stop('na')
    
    lst[[t]][['a']][1:(9*R)] <- a
    lst[[t]][['b']][1:R]     <- b
    lst[[t]][['c']][1:(9*R)] <- c
  }
  
  js.prob <- do.call('rbind', lapply(1:length(lst), function(x) lst[[x]][['c']]))
  ms.prob <- rbind(c(sapply(1:R, function(j) sum(js.prob[2,c(j,R+j,2*R+j)])),
                     sapply(1:R, function(j) sum(js.prob[2,c(3*R+j,4*R+j,5*R+j)])),
                     sapply(1:R, function(j) sum(js.prob[2,c(6*R+j,7*R+j,8*R+j)]))),
                   cbind(sapply(1:R, function(j) rowSums(js.prob[-1,c(j,3*R+j,6*R+j)])),
                         sapply(1:R, function(j) rowSums(js.prob[-1,c(R+j,4*R+j,7*R+j)])),
                         sapply(1:R, function(j) rowSums(js.prob[-1,c(2*R+j,5*R+j,8*R+j)]))))
  
  ### Derivatives
  
  ### Beta's
  
  # Hessian
  H1 <- diag(colSums((js.prob[-1,1:R] + js.prob[-1,(R+1):(2*R)] +
                        js.prob[-1,(2*R+1):(3*R)]) * p00[-1,] * (1-p00[-1,]) + 
                       (js.prob[-1,(3*R+1):(4*R)] + js.prob[-1,(4*R+1):(5*R)]) * 
                       p10[-1,]*(1-p10[-1,]) + 
                       (js.prob[-1,(6*R+1):(7*R)] + js.prob[-1,(8*R+1):(9*R)]) * 
                       p20[-1,]*(1-p20[-1,])))
  H2 <- Q
  H <- H1 + H2
  
  sgm.b01.log.pt00 <- sapply(1:dim(xx01)[2], function(j) 
    - sum(xx01[-1,j,] * p01[-1,] * js.prob[-1,1:R]))
  sgm.b01.log.pt01 <- sapply(1:dim(xx01)[2], function(j) 
    sum(xx01[-1,j,] * (1 - p01[-1,]) * js.prob[-1,(R+1):(2*R)]))
  sgm.b01.log.pt02 <- sapply(1:dim(xx01)[2], function(j) 
    - sum(xx01[-1,j,] * p01[-1,] * js.prob[-1,(2*R+1):(3*R)]))
  
  sgm.b02.log.pt00 <- sapply(1:dim(xx02)[2], function(j) 
    - sum(xx02[-1,j,] * p02[-1,] * js.prob[-1,1:R]))
  sgm.b02.log.pt01 <- sapply(1:dim(xx02)[2], function(j) 
    - sum(xx02[-1,j,] * p02[-1,] * js.prob[-1,(R+1):(2*R)]))
  sgm.b02.log.pt02 <- sapply(1:dim(xx02)[2], function(j) 
    sum(xx02[-1,j,] * (1 - p02[-1,]) * js.prob[-1,(2*R+1):(3*R)]))
  
  sgm.b11.log.pt10 <- sapply(1:dim(xx11)[2], function(j) 
    - sum(xx11[-1,j,] * p11[-1,] * js.prob[-1,(3*R+1):(4*R)]))
  sgm.b11.log.pt11 <- sapply(1:dim(xx11)[2], function(j) 
    sum(xx11[-1,j,] * (1 - p11[-1,]) * js.prob[-1,(4*R+1):(5*R)]))
  
  sgm.b22.log.pt20 <- sapply(1:dim(xx22)[2], function(j) 
    - sum(xx22[-1,j,] * p22[-1,] * js.prob[-1,(6*R+1):(7*R)]))
  sgm.b22.log.pt22 <- sapply(1:dim(xx22)[2], function(j) 
    sum(xx22[-1,j,] * (1 - p22[-1,]) * js.prob[-1,(8*R+1):(9*R)]))
  
  sgm.Hes.b01 <- sapply(1:dim(xx01)[2], function(j) 
    colSums(exp.x1b01[-1,] * (1 - exp.x1b01[-1,] - exp.x2b02[-1,])/
              ((1 + exp.x1b01[-1,] + exp.x2b02[-1,])^3) * 
              xx01[-1,j,] * (js.prob[-1,1:R] + js.prob[-1,(R+1):(2*R)] + js.prob[-1,(2*R+1):(3*R)])))
  sgm.Hes.b02 <- sapply(1:dim(xx02)[2], function(j) 
    colSums(exp.x2b02[-1,] * (1 - exp.x1b01[-1,] - exp.x2b02[-1,])/
              ((1 + exp.x1b01[-1,] + exp.x2b02[-1,])^3) * 
              xx02[-1,j,] * (js.prob[-1,1:R] + js.prob[-1,(R+1):(2*R)] + js.prob[-1,(2*R+1):(3*R)])))
  
  sgm.Hes.b11 <- sapply(1:dim(xx11)[2], function(j) 
    colSums(exp.x1b11[-1,] * (1 - exp.x1b11[-1,])/((1 + exp.x1b11[-1,])^3) * 
              xx11[-1,j,] * (js.prob[-1,(3*R+1):(4*R)] + js.prob[-1,(4*R+1):(5*R)])))
  
  sgm.Hes.b22 <- sapply(1:dim(xx22)[2], function(j) 
    colSums(exp.x2b22[-1,] * (1 - exp.x2b22[-1,])/((1 + exp.x2b22[-1,])^3) * 
              xx22[-1,j,] * (js.prob[-1,(6*R+1):(7*R)] + js.prob[-1,(8*R+1):(9*R)])))
  
  Hinv <- ginv(H)
  hb01 <- -0.5 * colSums(sweep(sgm.Hes.b01, 1, diag(Hinv), '*'))
  hb02 <- -0.5 * colSums(sweep(sgm.Hes.b02, 1, diag(Hinv), '*'))
  hb11 <- -0.5 * colSums(sweep(sgm.Hes.b11, 1, diag(Hinv), '*'))
  hb22 <- -0.5 * colSums(sweep(sgm.Hes.b22, 1, diag(Hinv), '*'))
  
  
  res3 <- c(sgm.b01.log.pt00 + sgm.b01.log.pt01 + sgm.b01.log.pt02 + hb01, 
            sgm.b02.log.pt00 + sgm.b02.log.pt01 + sgm.b02.log.pt02 + hb02, 
            sgm.b11.log.pt10 + sgm.b11.log.pt11 + hb11,
            sgm.b22.log.pt20 + sgm.b22.log.pt22 + hb22)
  
  # Derivatives
  res3
}

# Check if gradient is approx zero around optimal parameter values
S.alpha(alpha = alpha.opt, b01 = beta01.opt, b02 = beta02.opt,
        b1 = beta11.opt, b2 = beta22.opt, dtsim = dt) %>% unname()
S.beta(beta = beta.opt, a1 = alpha1.opt, a2 = alpha2.opt, 
       g = gamma.opt, dtsim = dt)

##### 3) Fisher information matrix ----

# Sample death counts from sampled RS trajectory (Function)
gen_st_dt <- function(r) {
  # start in baseline state 
  St <- matrix(NA, nrow = nT, ncol = 1)
  St[1,] <- 0
  
  # Sample paths regime states
  for(t in 2:nT){
    j0 <- which(St[t-1,] == 0)
    j1 <- which(St[t-1,] == 1)
    j2 <- which(St[t-1,] == 2)
    
    nj0 <- length(j0)
    nj1 <- length(j1)
    nj2 <- length(j2)
    
    s0 <- sample(0:2, size = nj0, replace = TRUE, prob = c(pt.00[t,r], pt.01[t,r], pt.02[t,r]))
    s1 <- sample(0:2, size = nj1, replace = TRUE, prob = c(pt.10[t,r], pt.11[t,r], 0))
    s2 <- sample(0:2, size = nj2, replace = TRUE, prob = c(pt.20[t,r], 0, pt.22[t,r]))
    
    St[t,j0] <- s0
    St[t,j1] <- s1
    St[t,j2] <- s2
  }
  
  # Death counts
  k0  <- do.call('rbind', replicate(nA, St == 0, simplify = FALSE))
  k1  <- do.call('rbind', replicate(nA, St == 1, simplify = FALSE))
  k2  <- do.call('rbind', replicate(nA, St == 2, simplify = FALSE))
  
  x1a1  <- x1[,,r] %*% alpha1.opt
  x2a2  <- x2[,,r] %*% alpha2.opt
  
  dt0 <- rpois(sum(k0), lambda = bt[row(k0)[k0],r])
  dt1 <- rpois(sum(k1), lambda = bt[row(k1)[k1],r] * exp(x1a1[row(k1)[k1],]))
  dt2 <- rpois(sum(k2), lambda = bt[row(k2)[k2],r] * exp(x2a2[row(k2)[k2],]))
  
  Dt  <- matrix(NA, nrow = nT*nA, ncol = 1)
  Dt[k0] <- dt0
  Dt[k1] <- dt1
  Dt[k2] <- dt2
  
  list('St' = St, 'Dt' = Dt)
}

# Number of bootstrap samples
B <- 25000

# Simultaneous perturbation stochastic approximation to estimate Hessian
cl <- snow::makeCluster(detectCores() - 1)
registerDoSNOW(cl)
snow::clusterEvalQ(cl, c(library(dplyr), library(MASS)))
snow::clusterExport(cl, list = c('nT','nA','B', 'R', 'd', 'x1','x2', 'alpha1.opt',
                                 'alpha2.opt', 'dt', 'bt', 'theta.opt', 'beta.opt',
                                 'gen_st_dt', 'pen01', 'pen02','pen11','pen22'))

pb <- progress_bar$new(format = "Hessian = :letter [:bar] :elapsed | eta: :eta",
                       total = B, width = 60)
progress_letter <- 1:B
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
}
opts <- list(progress = progress)

obj <- foreach(k = 1:B, .options.snow = opts) %dopar% {
  # Calculate samples for death counts
  Dt    <- lapply(1:R, gen_st_dt)
  Dt    <- do.call('cbind', lapply(Dt, function(x) x$Dt[,1]))
  
  # Perturbations
  eps1 <- 10^(-6)
  eps2 <- 10^(-8)
  
  # Perturbations
  Dlt1 <- (rbinom(n = length(alpha.opt), size = 1, p = 0.5)*2 - 1)*eps1
  Dlt2 <- (rbinom(n = length(beta.opt), size = 1, p = 0.5)*2 - 1)*eps2
  
  # S-vector
  Sk1pD <- S.alpha(alpha = alpha.opt + Dlt1, dtsim = Dt, b01 = beta01.opt,
                   b02 = beta02.opt, b1 = beta11.opt, b2 = beta22.opt)
  Sk1mD <- S.alpha(alpha = alpha.opt - Dlt1, dtsim = Dt, b01 = beta01.opt,
                   b02 = beta02.opt, b1 = beta11.opt, b2 = beta22.opt)
  Sk2pD <- S.beta(beta = beta.opt + Dlt2, dtsim = Dt, a1 = alpha1.opt,
                  a2 = alpha2.opt, g = gamma.opt)
  Sk2mD <- S.beta(beta = beta.opt - Dlt2, dtsim = Dt, a1 = alpha1.opt,
                  a2 = alpha2.opt, g = gamma.opt)
  
  # Differences
  dSk1 <- Sk1pD - Sk1mD
  dSk2 <- Sk2pD - Sk2mD
  
  # Hessian
  Hk1   <- - 0.5 * ((dSk1/2) %*% t(1/(Dlt1)) + t((dSk1/2) %*% t(1/Dlt1)))
  Hk2   <- - 0.5 * ((dSk2/2) %*% t(1/(Dlt2)) + t((dSk2/2) %*% t(1/Dlt2)))
  
  list('alpha' = Hk1, 'beta' = Hk2)
}
stopCluster(cl)

# Save
saveRDS(obj, file = 'Results/FI.rds')