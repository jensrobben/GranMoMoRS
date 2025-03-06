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


##### 1) Load baseline fit ----
fitb   <- readRDS('Results/GAM.rds')
dgamb  <- readRDS('Results/df_GAM.rds') %>% 
  dplyr::arrange(Region, Age, ISOYear, ISOWeek) %>%
  dplyr::filter(Date %in% unique(df$Date))

Mu     <- fitb$coefficients %>% unname()
Sigma  <- vcov(fitb) %>% unname()
Expo   <- dgamb$Expo
Deaths <- dgamb$Deaths
Xd     <- dgamb %>% dplyr::select(-c('Date', 'ISOYear', 'ISOWeek', 'Region', 'Age', 
                                     'Deaths', 'Expo', 'fsin52', 'fcos52', 'fsin26', 
                                     'fcos26', 'Time')) %>% SparseM::as.matrix()
rm(fitb)
gc()
##### 2) Parameter estimates RS model -----

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

##### 3) Parameter uncertainty using Fisher information ----

# Results
FI <- readRDS(file = 'Results/FI.rds')

# Covariance matrices
B <- 25000
FI.alpha <- Reduce('+',lapply(FI, function(x) x$alpha))/B
FI.beta  <- Reduce('+',lapply(FI, function(x) x$beta))/B

COV.alpha <- ginv(FI.alpha)
COV.beta  <- ginv(FI.beta)

SD.alpha <- sqrt(diag(COV.alpha))
SD.beta  <- sqrt(diag(COV.beta))

##### 4) State uncertainty only ----

# Number of bootstrap samples
B <- 25000

cl <- snow::makeCluster(detectCores() - 1) 
registerDoSNOW(cl)
snow::clusterEvalQ(cl, c(library(dplyr), library(MASS)))
snow::clusterExport(cl, list = c('Mu','Sigma','Expo','Deaths','Xd'))   

pb <- progress_bar$new(format = "Hessian = :letter [:bar] :elapsed | eta: :eta",
                       total = B, width = 60)
progress_letter <- 1:B 
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 
opts <- list(progress = progress)

obj <- foreach(k = 1:B, .options.snow = opts, .inorder = FALSE) %dopar% {
  # Simulate new baseline 
  btsim <- bt
  
  # Calculate transition probabilities, marginal state probs
  alpha.sim <- alpha.opt
  beta.sim  <- beta.opt
  
  a1  <- alpha.sim[1:length(alpha1.opt)]
  a2  <- alpha.sim[-c(1:length(alpha1.opt))]
  b0  <- beta.sim[1:length(beta0.opt)]
  b01 <- b0[1:length(beta01.opt)]
  b02 <- b0[-c(1:length(beta01.opt))]
  b11 <- beta.sim[(length(beta0.opt)+1):(length(beta0.opt)+length(beta11.opt))]
  b22 <- beta.sim[-c(1:length(c(beta0.opt,beta11.opt)))]
  
  ### (1.1) Conditional densities (Poisson)
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  
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
  
  # Simulate RS trajectory
  
  St <- matrix(NA, nrow = nT, ncol = R)
  St[1,] <- 0
  
  # Sample paths regime states
  for(t in 2:nT){
    j0 <- which(St[t-1,] == 0)
    j1 <- which(St[t-1,] == 1)
    j2 <- which(St[t-1,] == 2)
    
    nj0 <- length(j0)
    nj1 <- length(j1)
    nj2 <- length(j2)
    
    s0 <- sapply(seq_along(j0), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p00[t,j0[k]],
                                                            p01[t,j0[k]],
                                                            p02[t,j0[k]])))
    s1 <- sapply(seq_along(j1), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p10[t,j1[k]], 
                                                            p11[t,j1[k]], 0)))
    s2 <- sapply(seq_along(j2), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p20[t,j2[k]], 0, 
                                                            p22[t,j2[k]])))
    
    St[t,j0] <- unlist(s0)
    St[t,j1] <- unlist(s1)
    St[t,j2] <- unlist(s2)
  }
  
  # Death counts
  k0  <- do.call('rbind', replicate(nA, St == 0, simplify = FALSE))
  k1  <- do.call('rbind', replicate(nA, St == 1, simplify = FALSE))
  k2  <- do.call('rbind', replicate(nA, St == 2, simplify = FALSE))
  
  dt0 <- btsim[k0]
  dt1 <- btsim[k1] * exp(x1a1[k1])
  dt2 <- btsim[k2] * exp(x2a2[k2])
  
  Dt  <- matrix(NA, nrow = nT*nA, ncol = R)
  Dt[k0] <- dt0
  Dt[k1] <- dt1
  Dt[k2] <- dt2
  
  as.vector(Dt)
}

stopCluster(cl)

res  <- as.data.table(obj) %>% as.matrix()
rm(obj)
gc()
qmat <- matrixStats::rowQuantiles(res, probs = c(0.025,0.5,0.975), na.rm = TRUE)
rm(res)
gc()

df.UC1 <- cbind(dgamb, qmat) %>% 
  dplyr::mutate('Method' = 'Method1') %>%
  dplyr::select(c('Date','ISOYear','ISOWeek','Region','Age','Deaths','Expo',
                  '2.5%','50%','97.5%','Method'))

saveRDS(df.UC1,'Results/UncQ1.rds')

##### 5) Parameter + spatial + state uncertainty ----

# Number of bootstrap samples
B <- 25000

cl <- snow::makeCluster(detectCores() - 1) 
registerDoSNOW(cl)
snow::clusterEvalQ(cl, c(library(dplyr), library(MASS)))
snow::clusterExport(cl, list = c('Mu','Sigma','Expo','Deaths','Xd'))   
pb <- progress_bar$new(format = "Hessian = :letter [:bar] :elapsed | eta: :eta",
                       total = B, width = 60)
progress_letter <- 1:B 
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 
opts <- list(progress = progress)
obj <- foreach(k = 1:B, .options.snow = opts, .inorder = FALSE) %dopar% {
  # Simulate new baseline 
  Mub   <- mvnfast::rmvn(1, Mu, Sigma)
  Dxtb  <- exp(Xd %*% t(Mub))*Expo
  btsim <- matrix(Dxtb[,1], ncol = R)
  Usim  <- mvtnorm::rmvnorm(1, rep(0,R), cov)
  
  # Calculate transition probabilities, marginal state probs
  alpha.sim <- mvnfast::rmvn(n = 1, mu = c(alpha1.opt, alpha2.opt),
                             sigma = COV.alpha)
  beta.sim  <- mvnfast::rmvn(n = 1, mu = c(beta0.opt, beta11.opt, beta22.opt),
                             sigma = COV.beta)
  
  a1  <- alpha.sim[1:length(alpha1.opt)]
  a2  <- alpha.sim[-c(1:length(alpha1.opt))]
  b0  <- beta.sim[1:length(beta0.opt)]
  b01 <- b0[1:length(beta01.opt)]
  b02 <- b0[-c(1:length(beta01.opt))]
  b11 <- beta.sim[(length(beta0.opt)+1):(length(beta0.opt)+length(beta11.opt))]
  b22 <- beta.sim[-c(1:length(c(beta0.opt,beta11.opt)))]
  
  ### Matrices
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  x1b01 <- pen01 + sapply(1:R, function(r) xx01[,,r] %*% b01 + Usim[r])
  x2b02 <- pen02 + sapply(1:R, function(r) xx02[,,r] %*% b02 + Usim[r]) 
  x1b11 <- pen11 + sapply(1:R, function(r) xx11[,,r] %*% b11 + Usim[r])
  x2b22 <- pen22 + sapply(1:R, function(r) xx22[,,r] %*% b22 + Usim[r]) 
  
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
  
  # Simulate RS trajectory
  
  St <- matrix(NA, nrow = nT, ncol = R)
  St[1,] <- 0
  
  # Sample paths regime states
  for(t in 2:nT){
    j0 <- which(St[t-1,] == 0)
    j1 <- which(St[t-1,] == 1)
    j2 <- which(St[t-1,] == 2)
    
    nj0 <- length(j0)
    nj1 <- length(j1)
    nj2 <- length(j2)
    
    s0 <- sapply(seq_along(j0), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p00[t,j0[k]],
                                                            p01[t,j0[k]],
                                                            p02[t,j0[k]])))
    s1 <- sapply(seq_along(j1), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p10[t,j1[k]], 
                                                            p11[t,j1[k]], 0)))
    s2 <- sapply(seq_along(j2), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p20[t,j2[k]], 0, 
                                                            p22[t,j2[k]])))
    
    St[t,j0] <- unlist(s0)
    St[t,j1] <- unlist(s1)
    St[t,j2] <- unlist(s2)
  }
  
  # Death counts
  k0  <- do.call('rbind', replicate(nA, St == 0, simplify = FALSE))
  k1  <- do.call('rbind', replicate(nA, St == 1, simplify = FALSE))
  k2  <- do.call('rbind', replicate(nA, St == 2, simplify = FALSE))
  
  dt0 <- btsim[k0]
  dt1 <- btsim[k1] * exp(x1a1[k1])
  dt2 <- btsim[k2] * exp(x2a2[k2])
  
  Dt  <- matrix(NA, nrow = nT*nA, ncol = R)
  Dt[k0] <- dt0
  Dt[k1] <- dt1
  Dt[k2] <- dt2
  
  as.vector(Dt)
}
stopCluster(cl)

res  <- as.data.table(obj) %>% as.matrix()
rm(obj)
gc()
qmat <- matrixStats::rowQuantiles(res, probs = c(0.025,0.5,0.975), na.rm = TRUE)
rm(res)
gc()

df.UC2 <- cbind(dgamb, qmat) %>% 
  dplyr::mutate('Method' = 'Method2') %>%
  dplyr::select(c('Date','ISOYear','ISOWeek','Region','Age','Deaths','Expo',
                  '2.5%','50%','97.5%','Method'))

# Save
saveRDS(df.UC2, 'Results/UncQ2.rds')

##### 6) Parameter + spatial + state + Poisson uncertainty ----

# Number of bootstrap samples
B <- 25000

cl <- snow::makeCluster(detectCores() - 1) 
registerDoSNOW(cl)
snow::clusterEvalQ(cl, c(library(dplyr), library(MASS)))
snow::clusterExport(cl, list = c('Mu','Sigma','Expo','Deaths','Xd'))   

pb <- progress_bar$new(format = "Hessian = :letter [:bar] :elapsed | eta: :eta",
                       total = B, width = 60)
progress_letter <- 1:B 
progress <- function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 
opts <- list(progress = progress)

obj <- foreach(k = 1:B, .options.snow = opts, .inorder = FALSE) %dopar% {
  # Simulate new baseline 
  Mub   <- mvnfast::rmvn(1, Mu, Sigma)
  Dxtb  <- exp(Xd %*% t(Mub))*Expo # No poisson uncertainty in baseline - otherwise double
  btsim <- matrix(Dxtb, ncol = R)
  Usim  <- mvtnorm::rmvnorm(1, rep(0,R), cov)
  
  # Calculate transition probabilities, marginal state probs
  alpha.sim <- mvnfast::rmvn(n = 1, mu = c(alpha1.opt, alpha2.opt),
                             sigma = COV.alpha)
  beta.sim  <- mvnfast::rmvn(n = 1, mu = c(beta0.opt, beta11.opt, beta22.opt),
                             sigma = COV.beta)
  
  a1  <- alpha.sim[1:length(alpha1.opt)]
  a2  <- alpha.sim[-c(1:length(alpha1.opt))]
  b0  <- beta.sim[1:length(beta0.opt)]
  b01 <- b0[1:length(beta01.opt)]
  b02 <- b0[-c(1:length(beta01.opt))]
  b11 <- beta.sim[(length(beta0.opt)+1):(length(beta0.opt)+length(beta11.opt))]
  b22 <- beta.sim[-c(1:length(c(beta0.opt,beta11.opt)))]
  
  ### Matrices
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  x1b01 <- pen01 + sapply(1:R, function(r) xx01[,,r] %*% b01 + Usim[r])
  x2b02 <- pen02 + sapply(1:R, function(r) xx02[,,r] %*% b02 + Usim[r]) 
  x1b11 <- pen11 + sapply(1:R, function(r) xx11[,,r] %*% b11 + Usim[r])
  x2b22 <- pen22 + sapply(1:R, function(r) xx22[,,r] %*% b22 + Usim[r]) 
  
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
  
  # Simulate RS trajectory
  
  St <- matrix(NA, nrow = nT, ncol = R)
  St[1,] <- 0
  
  # Sample paths regime states
  for(t in 2:nT){
    j0 <- which(St[t-1,] == 0)
    j1 <- which(St[t-1,] == 1)
    j2 <- which(St[t-1,] == 2)
    
    nj0 <- length(j0)
    nj1 <- length(j1)
    nj2 <- length(j2)
    
    s0 <- sapply(seq_along(j0), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p00[t,j0[k]],
                                                            p01[t,j0[k]],
                                                            p02[t,j0[k]])))
    s1 <- sapply(seq_along(j1), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p10[t,j1[k]], 
                                                            p11[t,j1[k]], 0)))
    s2 <- sapply(seq_along(j2), function(k) sample(0:2, size = 1, replace = TRUE,
                                                   prob = c(p20[t,j2[k]], 0, 
                                                            p22[t,j2[k]])))
    
    St[t,j0] <- unlist(s0)
    St[t,j1] <- unlist(s1)
    St[t,j2] <- unlist(s2)
  }
  
  # Death counts
  k0  <- do.call('rbind', replicate(nA, St == 0, simplify = FALSE))
  k1  <- do.call('rbind', replicate(nA, St == 1, simplify = FALSE))
  k2  <- do.call('rbind', replicate(nA, St == 2, simplify = FALSE))
  
  dt0 <- rpois(sum(k0), lambda = btsim[k0])
  dt1 <- rpois(sum(k1), lambda = btsim[k1] * exp(x1a1[k1]))
  dt2 <- rpois(sum(k2), lambda = btsim[k2] * exp(x2a2[k2]))
  
  Dt  <- matrix(NA, nrow = nT*nA, ncol = R)
  Dt[k0] <- dt0
  Dt[k1] <- dt1
  Dt[k2] <- dt2
  
  as.vector(Dt)
}

stopCluster(cl)

res  <- as.data.table(obj) %>% as.matrix()
rm(obj)
gc()
qmat <- matrixStats::rowQuantiles(res, probs = c(0.025,0.5,0.975), na.rm = TRUE)
rm(res)
gc()

df.UC3 <- cbind(dgamb, qmat) %>% 
  dplyr::mutate('Method' = 'Method3') %>%
  dplyr::select(c('Date','ISOYear','ISOWeek','Region','Age','Deaths','Expo',
                  '2.5%','50%','97.5%','Method'))
# Save
saveRDS(df.UC3, 'Results/UncQ3.rds')

