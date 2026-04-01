##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'stats',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr','httr',
              'utils', 'sp', 'spdep', 'RColorBrewer')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Colour scales
ablue  <- '#56B4E9'
agreen <- '#A3D596'
ared   <- '#F8766D'

##### 1) Load shape file NUTS 2 European regions ----

# Read shape file - downloaded from Eurostat
shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp')

# Extract shapefile for NUTS 2 regions in the countries of interest 
nuts.spec <- 2
ctry.spec <- 'FR'
shapef    <- shapefile[shapefile$CNTR_CODE %in% ctry.spec & 
                         shapefile$LEVL_CODE == nuts.spec,]

# Remove overseas areas
overseas  <- c(paste0('FRY',1:5), 'FRM0')
ind.rm    <- unlist(sapply(overseas, function(x) 
  which(grepl(x, shapef$NUTS_ID, fixed = TRUE))))
shapef    <- shapef[-c(ind.rm), ] %>%
  dplyr::arrange(NUTS_ID)

# Regions
regions <- sort(unique(shapef$NUTS_ID))


##### 2) Load data sets and preprocess -----

# Merge mortality data with feature data set
df.feat    <- readRDS('Results/df_NUTS2_FR.rds')
df.bDeaths <- readRDS('Results/df_bDeaths.rds')

df <- df.bDeaths %>% 
  left_join(df.feat, 
            by = c('Date', 'ISOYear', 'ISOWeek', 'Region')) %>%
  na.omit()

# Quantile excess hospital admission
df$Nhospq  <- pmax(df$Nhosp - quantile(df$Nhosp[df$Date > as.Date('2020-03-16')], 0.75))

# Filter dataset and select threshold quantile
vars <- c('ia100q', 'Nhospq', 'w_avg_ehi2', 'w_avg_eci2')
df   <- df %>%
  dplyr::select(all_of(c('Date', 'ISOYear', 'ISOWeek', 'Region', 'geo_name', 
                         'Age', 'Deaths', 'Expo', 'bDeaths', 'avg_ia100', vars)))
colnames(df)[grepl('w_avg_', colnames(df), fixed = TRUE)] <- c('w_avg_ehi', 'w_avg_eci')

# Create lagged features
vars <- c('ia100q', 'Nhospq', 'w_avg_ehi', 'w_avg_eci')
for(v in vars){
  df <- df %>% 
    group_by(Region, Age) %>%
    mutate(!!paste0(v,'_l1')   := lag(!!sym(v), 1),
           !!paste0(v,'_l2.3') := (lag(!!sym(v), 2) + lag(!!sym(v), 3))/2) %>% 
    ungroup()
}

# State restrictions (from state 0 -> 1)
df$penalty01 <- ifelse(df$w_avg_ehi > 0, 0, -1000)
df$penalty02 <- ifelse(df$ia100q > 0 | df$Nhospq > 0 | df$w_avg_eci > 0 | 
                         df$ia100q_l1 | df$Nhospq_l1 > 0 | df$w_avg_eci_l1 > 0,
                       0, -1000)

# Remove missing observations (weeks 1-3 of year 2019)
df <- df %>% na.omit()

# Arrange data frame by (Region, Age, Time)
df <- df %>% 
  dplyr::arrange(Region, Age, ISOYear, ISOWeek) %>%
  mutate(Region = as.factor(Region),
         Time = ISOYear - min(ISOYear))

##### 3) ICAR process for spatial effect ----

# Adjacency matrix
adjmat             <- nb2mat(poly2nb(shapef$geometry), zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat)   <- list(shapef$NUTS_ID, shapef$NUTS_ID)
diag(adjmat)       <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))

# ICAR set-up
D   <- diag(diag(adjmat))
W   <- - (adjmat - D)
tau <- 1
Q   <- tau * (D - W)
cov <- MASS::ginv(Q)

# Function to calculate Moore-Penroze generalized determinant
det.mp <- function(x) {
  sigma <- zapsmall(svd(x)$d)  
  prod(sigma[sigma != 0])       
}

##### 4) Design matrices -----

# Define parameters
R   <- length(unique(df$Region))
nA  <- length(unique(df$Age))
nT  <- nrow(df)/(nA*R)

# State-specific Poisson model specifications
# State 1
form1.0 <- as.formula(Deaths ~ offset(log(bDeaths)) + 
                        w_avg_ehi:Age + 
                        w_avg_ehi_l1:Age)
fit1.0  <- glm(form1.0, 
               data = df, family = poisson(link = 'log'))
m10     <- model.matrix(fit1.0)
x1      <- simplify2array(lapply(1:R, function(j) 
  matrix(m10[(nA*nT*(j-1) + 1):(nA*nT*j),], ncol = ncol(m10))))
alpha.1 <- unname(fit1.0$coefficients)

# State 2
form2.0 <- as.formula(Deaths ~ offset(log(bDeaths)) + 
                        w_avg_eci:Age + 
                        ia100q:Age + 
                        Nhospq:Age + 
                        w_avg_eci_l1:Age +
                        ia100q_l1:Age + 
                        Nhospq_l1:Age + 
                        w_avg_eci_l2.3:Age + 
                        ia100q_l2.3:Age +
                        w_avg_eci:Nhospq:Age + 
                        w_avg_eci_l1:Nhospq_l1:Age)
fit2.0  <- glm(form2.0, data = df, family = poisson(link = 'log'))
m20     <- model.matrix(fit2.0)
x2      <- simplify2array(lapply(1:R, function(j) 
  matrix(m20[(nA*nT*(j-1) + 1):(nA*nT*j),], ncol = ncol(m20))))
alpha.2 <- unname(fit2.0$coefficients)

# Transition model matrices
m01 <- model.matrix(Deaths ~ offset(log(bDeaths)) + w_avg_ehi,
                    data = df %>% dplyr::filter(Age == 'Y_GE90'))
xx01 <- simplify2array(lapply(1:R, function(j) 
  matrix(m01[(nT*(j-1) + 1):(nT*j),], ncol = ncol(m01))))

m11 <- model.matrix(Deaths ~ offset(log(bDeaths)) + w_avg_ehi,
                    data = df %>% dplyr::filter(Age == 'Y_GE90'))
xx11 <- simplify2array(lapply(1:R, function(j) 
  matrix(m11[(nT*(j-1) + 1):(nT*j),], ncol = ncol(m11))))

m02 <- model.matrix(Deaths ~ offset(log(bDeaths)) + w_avg_eci,
                    data = df %>% dplyr::filter(Age == 'Y_GE90'))
xx02 <- simplify2array(lapply(1:R, function(j) 
  matrix(m02[(nT*(j-1) + 1):(nT*j),], ncol = ncol(m02))))

m22 <- model.matrix(Deaths ~ offset(log(bDeaths)) + w_avg_eci_l2.3,
                    data = df %>% dplyr::filter(Age == 'Y_GE90'))
xx22 <- simplify2array(lapply(1:R, function(j) 
  matrix(m22[(nT*(j-1) + 1):(nT*j),], ncol = ncol(m22))))

# Transition probs
beta.01 <- rep(0, dim(xx01)[2])
beta.02 <- rep(0, dim(xx02)[2])
beta.0  <- c(beta.01, beta.02)
beta.11 <- rep(0, dim(xx11)[2])
beta.22 <- rep(0, dim(xx22)[2])

# Initial probs
gamma.0 <- rep(0.3,R)
gamma.1 <- rep(0,R)
gamma.2 <- rep(0.7,R)
gamma   <- c(gamma.0, gamma.1, gamma.2)

##### 5) Objective functions for M-step ----

# Put information in matrix form (nT x R)
dt    <- matrix(df$Deaths, nrow = nT*nA, ncol = R, 
                dimnames = list(NULL, regions))
bt    <- matrix(df$bDeaths, nrow = nT*nA, ncol = R, 
                dimnames = list(NULL, regions))
pen01 <- matrix(df$penalty01, nrow = nT*nA, ncol = R, 
                dimnames = list(NULL, regions))
pen02 <- matrix(df$penalty02, nrow = nT*nA, ncol = R, 
                dimnames = list(NULL, regions))

# Objective functions to be minimized in EM-algorithm
obj.alpha.1 <- function(alpha, msprobs, xmat1 = x1[,,],
                        dtr = dt, btr = bt, pen = NULL){
  Pt1  <- do.call('rbind', (replicate(6, msprobs, simplify = FALSE)))
  x1a1 <- sapply(1:dim(xmat1)[3], function(j) t(t(xmat1[,,j])) %*% alpha)
  Et1  <- exp(x1a1)
  
  - sum(Pt1 * (-btr*Et1 + dtr*x1a1 + dtr*log(btr) - lgamma(dtr+1)))
}
obj.alpha.2 <- function(alpha, msprobs, xmat2 = x2[,,],
                        dtr = dt, btr = bt, pen = NULL){
  Pt2  <- do.call('rbind', (replicate(6, msprobs, simplify = FALSE)))
  x2a2 <- sapply(1:dim(xmat2)[3], function(j) t(t(xmat2[,,j])) %*% alpha)
  Et2  <- exp(x2a2)
  
  - sum(Pt2 * (-btr*Et2 + dtr*x2a2 + dtr*log(btr) - lgamma(dtr+1)))
}
obj.beta   <- function(beta, jsprobs, U, xmat01 = xx01, xmat02 = xx02,
                       xmat11 = xx11, xmat22 = xx22, dtr = dt, btr = bt, 
                       pena01 = pen01, pena02 = pen02, pena11 = pen11, 
                       pena22 = pen22){
  i01  <- dim(xmat01)[2]
  i11  <- dim(xmat11)[2]
  i02  <- dim(xmat02)[2]
  i22  <- dim(xmat22)[2]
  b01 <- beta[1:i01]
  b02 <- beta[(i01+1):(i01+i02)]
  b11 <- beta[(i01+i02+1):(i01+i02+i11)]
  b22 <- beta[(i01+i02+i11+1):(i01+i02+i11+i22)]
  x1b01 <- pena01 + sapply(1:R, function(r) t(t(xmat01[,,r])) %*% b01 + U[r])
  x2b02 <- pena02 + sapply(1:R, function(r) t(t(xmat02[,,r])) %*% b02 + U[r]) 
  x1b11 <- pena11 + sapply(1:R, function(r) t(t(xmat11[,,r])) %*% b11 + U[r])
  x2b22 <- pena22 + sapply(1:R, function(r) t(t(xmat22[,,r])) %*% b22 + U[r]) 
  
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
  
  # Hessian
  H1 <- diag(colSums((jsprobs[-1,1:R] + jsprobs[-1,(R+1):(2*R)] +
                        jsprobs[-1,(2*R+1):(3*R)]) * p00[-1,] * (1-p00[-1,]) +
                       (jsprobs[-1,(3*R+1):(4*R)] + jsprobs[-1,(4*R+1):(5*R)]) *
                       p10[-1,]*(1-p10[-1,]) + (jsprobs[-1,(6*R+1):(7*R)] +
                                                  jsprobs[-1,(8*R+1):(9*R)]) * p20[-1,]*(1-p20[-1,])))
  H2 <- Q
  
  
  
  - sum(jsprobs[-1,1:R] * log(p00[-1,]) + 
          jsprobs[-1,(R+1):(2*R)] * log(p01[-1,]) +
          jsprobs[-1,(2*R+1):(3*R)] * log(p02[-1,]) + 
          jsprobs[-1,(3*R+1):(4*R)] * log(p10[-1,]) + 
          jsprobs[-1,(4*R+1):(5*R)] * log(p11[-1,]) +
          jsprobs[-1,(6*R+1):(7*R)] * log(p20[-1,]) + 
          jsprobs[-1,(8*R+1):(9*R)] * log(p22[-1,])) + 
    0.5*log(det.mp(H1 + H2))
}
obj.gamma <- function(gamma, msprobs, dens1, dens2, dens3){
  g1   <- gamma
  g2   <- 0
  g3   <- 1 - g1 - g2
  gvec <- rep(c(g1, g2, g3), each = R)
  
  term.g0 <- - sum(msprobs * log(gvec + 10^(-300))) 
  term.g0
}

# Laplace-approximated log-likelihood (function of spatial effect vector)
Uopt <- function(U){
  # Specs
  nA <- length(unique(df$Age))
  nT <- length(unique(df$Date))
  R  <- length(U) + 1
  
  # Make sure random effect sums to zero
  U <- c(U, -sum(U))
  
  # Gamma
  # Indices corresponding to first observation
  id  <- 1 + nT*(0:(nA-1))
  
  g1 <- gamma[1]
  g2 <- 0 # do not start in state 1 (only in summer)
  g3 <- 1 - g1 - g2
  gamma <- rep(c(g1, g2, g3), each = R)
  term.g0 <- - sum(ms.prob[1,] * log(gamma + 10^(-300)))
  
  # Alpha 0
  Pt0     <- do.call(rbind, replicate(nA, ms.prob[,1:R], simplify = FALSE))
  term.a0 <- - sum(Pt0 * (-bt + dt * log(bt) - lgamma(dt + 1)))
  
  # Alpha 1
  term.a1 <- obj.alpha.1(alpha.1, ms.prob[,(R+1):(2*R)], xmat1 = x1,
                         dtr = dt, btr = bt)
  
  # Alpha 2
  term.a2 <- obj.alpha.2(alpha.2, ms.prob[,(2*R+1):(3*R)], xmat2 = x2,
                         dtr = dt, btr = bt)
  
  # Beta
  term.b <- obj.beta(beta = c(beta.0, beta.11, beta.22), jsprobs = js.prob, U = U,
                     xmat01 = xx01, xmat02 = xx02, xmat11 = xx11, xmat22 = xx22, 
                     dtr = dt, btr = bt, pena01 = pen01, pena02 = pen02,
                     pena11 = pen11, pena22 = pen22)
  
  # Total
  Reduce('+', list(term.g0, term.a0, term.a1, term.a2, term.b)) +
    0.5 * as.numeric(t(U) %*% Q %*% t(t(U))) + 
    0.5 * log(det.mp(cov))
}

##### 6) Run the EM algorithm! ----

# Object to save results
plst <- rep(list(NULL),7)
names(plst) <- c('alpha.1', 'alpha.2', 'gamma', 'beta.0', 'beta.11', 'beta.22', 'U')
obj <- list(plst)

# Create index list for for-loop in EM algorithm 
create_ind_lst <- function(i) {
  unlist(rep(list((3*R*i+1):(3*R*i+R), (3*R*i+R+1):(3*R*i+2*R), 
                  (3*R*i+2*R+1):(3*R*i+3*R)),each = 3))
}

indl0 <- create_ind_lst(0)
indl1 <- create_ind_lst(1)
indl2 <- create_ind_lst(2)

# Initial values (random)
ll0   <- 20*10^8
ll1   <- 10*10^8
iter  <- 0
U.new <- rep(0,R)

# EM loop
while(abs(ll0 - ll1) > 10^(-12)){
  # Update iteration and log-likelihood
  ll0   <- ll1
  iter  <- iter + 1
  
  # Conditional densities (Poisson)
  x1a1  <- sapply(1:dim(x1)[3], function(j) t(t(x1[,,j])) %*% alpha.1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) t(t(x2[,,j])) %*% alpha.2)
  fdts1 <- dpois(dt, lambda = bt)
  fdts2 <- dpois(dt, lambda = bt * exp(x1a1))
  fdts3 <- dpois(dt, lambda = bt * exp(x2a2))
  fdts  <- cbind(fdts1, fdts2, fdts3)
  
  # Change 0 values by very small positive value to avoid numerical issues
  fdts1[which(fdts1 == 0)] <- 10^(-300)
  fdts2[which(fdts2 == 0)] <- 10^(-300)
  fdts3[which(fdts3 == 0)] <- 10^(-300)
  fdts[which(fdts == 0)]   <- 10^(-300)
  
  # Transition model matrices
  x1b01 <- sapply(1:dim(xx01)[3], function(j) t(t(xx01[,,j])) %*% beta.01 + U.new[j])
  x2b02 <- sapply(1:dim(xx02)[3], function(j) t(t(xx02[,,j])) %*% beta.02 + U.new[j])
  x1b11 <- sapply(1:dim(xx11)[3], function(j) t(t(xx11[,,j])) %*% beta.11 + U.new[j])
  x2b22 <- sapply(1:dim(xx22)[3], function(j) t(t(xx22[,,j])) %*% beta.22 + U.new[j])
  
  # Calculate the transition probabilities
  pen01  <- pen01[1:nT,]
  pen02  <- pen02[1:nT,]
  pen11  <- pen01*0
  pen22  <- pen02*0
  pt.01 <- exp(pen01 + x1b01)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
  pt.02 <- exp(pen02 + x2b02)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
  pt.00 <- 1 - pt.01 - pt.02
  
  pt.11 <- exp(pen11 + x1b11)/(1 + exp(pen11 + x1b11))
  pt.12 <- matrix(0, nrow = nT, ncol = R)
  pt.10 <- 1 - pt.11
  
  pt.21 <- matrix(0, nrow = nT, ncol = R)
  pt.22 <- exp(pen22 + x2b22)/(1 + exp(pen22 + x2b22))
  pt.20 <- 1 - pt.22
  
  P <- unname(cbind(pt.00, pt.01, pt.02, pt.10, pt.11, pt.12, pt.20, pt.21, pt.22))
  P[1,] <- NA
  
  # Calculate joint state probabilities (check paper for iterative procedure)
  lst <- rep(list(list('a' = rep(NA, R*9),
                       'b' = rep(NA, R),
                       'c' = rep(NA, R*9))), nT)
  
  for(t in 2:nT){
    if(t == 2) {
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        c(rep(gamma[1:R], times = 3), rep(gamma[(R+1):(2*R)], times = 3),
          rep(gamma[(2*R+1):(3*R)], times = 3))
    } else {
      lst_slice <- lst[[t-1]][['c']]
      a <- rep(Rfast::colprods(fdts[t + nT*(1:nA-1),]), times = 3) * P[t,] * 
        (lst_slice[indl0] + lst_slice[indl1] + lst_slice[indl2])
    }
    
    b <- rowSums(matrix(a, nrow = R))
    b[b == 0] <- 10^(-300)
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
  
  ###### Maximization step 
  
  # Indices corresponding to first observation
  id  <- 1 + nT*(0:(nA-1))
  
  # Alpha 0
  Pt0     <- do.call(rbind, replicate(nA, ms.prob[,1:R], simplify = FALSE))
  term.a0 <- - sum(Pt0 * log(fdts1))
  
  # Gamma
  fit0 <- optim(par = gamma[1], fn = obj.gamma, msprobs = ms.prob[1,],
                dens1 = fdts1, dens2 = fdts2, dens3 = fdts3,
                control = list(factr = 10^(-14)),
                method = 'L-BFGS-B', lower = 0.001, upper = 0.999)
  g1 <- fit0$par
  g2 <- 0 # do not start in state 1 (only in summer)
  g3 <- 1 - g1 - g2
  gamma <- rep(c(g1, g2, g3), each = R)
  term.g0 <- fit0$value
  
  # Alpha 1
  fit1 <- glm(form1.0,
              family = poisson(link = 'log'), 
              data = df,
              start = alpha.1, 
              control = list(epsilon = 10^(-14)),
              weights = as.vector(do.call('rbind', replicate(nA, ms.prob[,(R+1):(2*R)], simplify = FALSE))))
  alpha.1 <- fit1$coefficients
  
  # Alpha 2
  fit2 <- glm(form2.0,
              family = poisson(link = 'log'), data = df,
              start = alpha.2,
              control = list(epsilon = 10^(-14)),
              weights = as.vector(do.call('rbind', replicate(nA, ms.prob[,(2*R+1):(3*R)], simplify = FALSE))))
  alpha.2 <- fit2$coefficients
  
  # Transition parameters
  fit3 <- optim(par = c(beta.0, beta.11, beta.22), fn = obj.beta, U = U.new,
                jsprobs = js.prob, xmat01 = xx01, xmat02 = xx02, xmat11 = xx11, 
                xmat22 = xx22, dtr = dt, btr = bt, pena01 = pen01,
                pena02 = pen02, pena11 = pen11, pena22 = pen22,
                method = 'BFGS', control = list(reltol = 10^(-14)))
  
  beta.0  <- fit3$par[1:length(beta.0)]
  beta.01 <- fit3$par[1:length(beta.01)]
  beta.02 <- fit3$par[(length(beta.01)+1):(length(beta.01) + length(beta.02))]
  beta.11 <- fit3$par[(length(beta.01)+length(beta.02)+1):(length(beta.01) + length(beta.02) + length(beta.11))]
  beta.22 <- fit3$par[(length(beta.01)+length(beta.02)+length(beta.11)+1):(length(beta.01) + length(beta.02) + length(beta.11) + length(beta.22))]
  
  if(abs(Uopt(U.new[-R]) - ll0)/ll0 > 10^(-8)){
    test <- nlm(f = Uopt, p = U.new[-R])
    U.new <- c(test$estimate, -sum(test$estimate))
    ll1   <- test$minimum
  } else{
   ll1 <- Uopt(U.new[-R])
  }
  
  # Store results
  obj[[iter]][['gamma']]   <- unique(gamma)
  obj[[iter]][['alpha.1']] <- alpha.1
  obj[[iter]][['alpha.2']] <- alpha.2
  obj[[iter]][['beta.0']]  <- beta.0
  obj[[iter]][['beta.11']] <- beta.11
  obj[[iter]][['beta.22']] <- beta.22
  obj[[iter]][['U']]       <- U.new
  obj <- append(obj, list(plst)) 
  
  print(c(abs(ll0-ll1), round(ll1, 10)))
}

# AIC, BIC
npar <- ncol(x1) + ncol(x2) + ncol(xx01) + ncol(xx02) + 
  ncol(xx11) + ncol(xx22)
bic  <- as.numeric(2*ll1 + log(nrow(df)) * (npar))
c(npar, bic)

# Save output + global environment
saveRDS(object = obj, file = paste0('Results/EMoutput_tau', tau, '.rds'))
save(list = ls(.GlobalEnv), file = paste0("Results/EMresult_tau", tau, ".Rdata"))
