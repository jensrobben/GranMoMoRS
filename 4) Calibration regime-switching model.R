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

# All times and ISO week/year information
t_min <- as.Date('2013-01-01', format = '%Y-%m-%d') # Start of ISO week 1 in 2013
t_max <- as.Date('2024-06-30', format = '%Y-%m-%d') # End of ISO week 26 in 2024

time_frame <- seq(t_min, t_max, by="days")
df_itime   <- data.frame('Date' = time_frame, 'ISOWeek' = isoweek(time_frame),
                         'ISOYear' = isoyear(time_frame))

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
df.feat    <- readRDS('Data/df_NUTS2_FR.rds')
df.bDeaths <- readRDS('Results/df_bDeaths.rds')

df <- df.bDeaths %>% 
  left_join(df.feat, 
            by = c('Date', 'ISOYear', 'ISOWeek', 'Region')) 

# Create lagged features
vars <- c('w_avg_tg_anom', 'w_avg_Tind95', 'w_avg_Tind5', 
          'ia100', 'ia100q', 'Nhosp', 'Nhospq')

for(v in vars){
  df <- df %>% 
    group_by(Region, Age) %>%
    mutate(!!paste0(v,'_l1')   := lag(!!sym(v), 1),
           !!paste0(v,'_l2')   := lag(!!sym(v), 2),
           !!paste0(v,'_l0.1') := (lag(!!sym(v), 0) + lag(!!sym(v), 1))/2,
           !!paste0(v,'_l2.3') := (lag(!!sym(v), 2) + lag(!!sym(v), 3))/2) %>% 
    ungroup()
}


# State restrictions (from state 0 -> 1)
df$penalty01 <- ifelse(df$w_avg_Tind95 > 0, 0, -1000)

# Remove missing observations (weeks 1-3 of year 2019)
df <- df %>% na.omit()

# Arrange data frame by (Region, Age, Time)
df <- df %>% 
  dplyr::arrange(Region, Age, Date) %>%
  mutate(Region = as.factor(Region))

##### 3) ICAR process for spatial effect ----

# Adjacency matrix
adjmat             <- nb2mat(poly2nb(shapef$geometry), zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat)   <- list(shapef$NUTS_ID, shapef$NUTS_ID)
diag(adjmat)       <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))

# ICAR set-up
D   <- diag(diag(adjmat))
W   <- - (adjmat - D)
tau <- 10
Q   <- tau * (D - W)
cov <- MASS::ginv(Q)

# Function to calculate Moore-Penroze generalized determinant
det.mp <- function(x) {
  sigma <- zapsmall(svd(x)$d)  
  prod(sigma[sigma != 0])       # Product of the nonzero singular values
}

##### 4) Design matrices -----

# Define parameters
R   <- length(unique(df$Region))
nA  <- length(unique(df$Age))
nT  <- nrow(df)/(nA*R)

# Group consecutive age groups (to limit model complexity)
df$Age0 <- df$Age
df$Age  <- plyr::mapvalues(df$Age0, from = c('Y65-69','Y70-74','Y75-79',
                                             'Y80-84','Y85-89','Y_GE90'),
                           to = c('Y65-74','Y65-74','Y75-84','Y75-84',
                                  'Y_GE85','Y_GE85'))
df$Age <- factor(df$Age, levels = c('Y65-74','Y75-84','Y_GE85'))

# State-specific Poisson model specifications
form1.0 <- as.formula(Deaths ~ -1 + offset(log(bDeaths)) + w_avg_tg_anom:Age + 
                        w_avg_tg_anom_l1:Age + w_avg_tg_anom_l2:Age)
m10     <- model.matrix(form1.0, data = df)
x1      <- simplify2array(lapply(1:R, function(j) m10[(nA*nT*(j-1) + 1):(nA*nT*j),]))

form2.0 <- as.formula(Deaths ~ -1 + offset(log(bDeaths)) + 
                        w_avg_Tind5_l0.1:Age +  w_avg_Tind5_l2.3:Age + 
                        ia100q_l0.1:Age + ia100q_l2.3:Age + 
                        Nhospq_l0.1:Age + Nhospq_l2.3:Age)
m20     <- model.matrix(form2.0, data = df)
x2      <- simplify2array(lapply(1:R, function(j) m20[(nA*nT*(j-1) + 1):(nA*nT*j),]))


# Transition model matrices (age-independent)
m01 <- model.matrix( ~ offset(log(bDeaths)) + w_avg_Tind95,
                    data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx01 <- simplify2array(lapply(1:R, function(j) 
  matrix(m01[(nT*(j-1) + 1):(nT*j),], ncol = ncol(m01))))

m11 <- model.matrix( ~ offset(log(bDeaths)) + 
                      w_avg_Tind95 + w_avg_Tind95_l1 + w_avg_Tind95_l2,
                    data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx11 <- simplify2array(lapply(1:R, function(j) m11[(nT*(j-1) + 1):(nT*j),]))

m02 <- model.matrix( ~ offset(log(bDeaths)) + 
                      ia100q_l0.1 + Nhospq_l0.1,
                    data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx02 <- simplify2array(lapply(1:R, function(j) m02[(nT*(j-1) + 1):(nT*j),]))

m22 <- model.matrix( ~ offset(log(bDeaths)) + 
                      ia100q_l0.1 + Nhospq_l0.1 +
                      ia100q_l2.3 + Nhospq_l2.3,
                    data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx22 <- simplify2array(lapply(1:R, function(j) m22[(nT*(j-1) + 1):(nT*j),]))

# Initial parameter estimates
alpha.1 <- rep(0, ncol(x1))
alpha.2 <- rep(0, ncol(x2))
beta.01 <- rep(0, dim(xx01)[2])
beta.02 <- rep(0, dim(xx02)[2])
beta.0  <- c(beta.01, beta.02)
beta.11 <- rep(0, dim(xx11)[2])
beta.22 <- rep(0, dim(xx22)[2])

# Initial marginal state probabilities
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

# Objective functions to be minimized in EM-algorithm
obj.alpha.1 <- function(alpha, msprobs, xmat1, dtr, btr){
  Pt1  <- do.call('rbind', (replicate(6, msprobs, simplify = FALSE)))
  x1a1 <- sapply(1:dim(xmat1)[3], function(j) xmat1[,,j] %*% alpha)
  Et1  <- exp(x1a1)
  
  - sum(Pt1 * (-btr*Et1 + dtr*x1a1 + dtr*log(btr) - lgamma(dtr+1)))
}
obj.alpha.2 <- function(alpha, msprobs, xmat2, dtr = dt, btr = bt){
  Pt2  <- do.call('rbind', (replicate(6, msprobs, simplify = FALSE)))
  x2a2 <- sapply(1:dim(xmat2)[3], function(j) xmat2[,,j] %*% alpha)
  Et2  <- exp(x2a2)
  
  - sum(Pt2 * (-btr*Et2 + dtr*x2a2 + dtr*log(btr) - lgamma(dtr+1)))
}
obj.beta   <- function(beta, jsprobs, U, xmat01 = xx01, xmat02 = xx02,
                       xmat11 = xx11, xmat22 = xx22, dtr = dt, btr = bt, 
                       pena01 = pen01, pena02 = 0, pena11 = 0, 
                       pena22 = 0){
  # Dimensions
  i01   <- dim(xmat01)[2]
  i11   <- dim(xmat11)[2]
  i02   <- dim(xmat02)[2]
  i22   <- dim(xmat22)[2]
  
  # Parameters
  b01   <- beta[1:i01]
  b02   <- beta[(i01+1):(i01+i02)]
  b11   <- beta[(i01+i02+1):(i01+i02+i11)]
  b22   <- beta[(i01+i02+i11+1):(i01+i02+i11+i22)]
  
  # Linear predictor
  x1b01 <- pena01 + sapply(1:R, function(r) xmat01[,,r] %*% b01 + U[r])
  x2b02 <- pena02 + sapply(1:R, function(r) xmat02[,,r] %*% b02 + U[r]) 
  x1b11 <- pena11 + sapply(1:R, function(r) xmat11[,,r] %*% b11 + U[r])
  x2b22 <- pena22 + sapply(1:R, function(r) xmat22[,,r] %*% b22 + U[r]) 
  
  # Exponentiate
  exp.x1b01 <- exp(x1b01) 
  exp.x2b02 <- exp(x2b02)
  exp.x1b11 <- exp(x1b11) 
  exp.x2b22 <- exp(x2b22) 
  
  # Transition probabilities
  p00 <- 1/(1 + exp.x1b01 + exp.x2b02)
  p01 <- exp.x1b01/(1 + exp.x1b01 + exp.x2b02)
  p02 <- exp.x2b02/(1 + exp.x1b01 + exp.x2b02)
  
  p10 <- 1/(1 + exp.x1b11)
  p11 <- exp.x1b11/(1 + exp.x1b11)
  
  p20 <- 1/(1 + exp.x2b22)
  p22 <- exp.x2b22/(1 + exp.x2b22)
  
  # Avoid numerical issues
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
  
  
  # Log-likelihood
  -sum(jsprobs[-1,1:R] * log(p00[-1,]) + 
         jsprobs[-1,(R+1):(2*R)] * log(p01[-1,]) +
         jsprobs[-1,(2*R+1):(3*R)] * log(p02[-1,]) + 
         jsprobs[-1,(3*R+1):(4*R)] * log(p10[-1,]) + 
         jsprobs[-1,(4*R+1):(5*R)] * log(p11[-1,]) +
         jsprobs[-1,(6*R+1):(7*R)] * log(p20[-1,]) + 
         jsprobs[-1,(8*R+1):(9*R)] * log(p22[-1,])) +
    0.5*log(det.mp(H1+H2))
}
obj.gamma <- function(gamma, msprobs, dens1, dens2, dens3){
  g1   <- gamma
  g2   <- 0
  g3   <- 1 - g1 - g2
  gvec <- rep(c(g1, g2, g3), each = R)
  
  ind.kp  <- which(! gvec == 0)
  term.g0 <- - sum(msprobs[ind.kp] * log(gvec[ind.kp])) 
  term.g0
}

# Laplace-approximated log-likelihood (function of spatial effect vector)
Uopt<- function(U){
  # Specs
  nA <- length(unique(df$Age0))
  nT <- length(unique(df$Date))
  R  <- length(U) + 1
  
  # Make sure random effect sums to zero
  U <- c(U, -sum(U))
  
  # Gamma
  id  <- 1 + nT*(0:(nA-1))
  g1 <- gamma[1]
  g2 <- 0 # do not start in state 1 (only in summer)
  g3 <- 1 - g1 - g2
  gamma <- rep(c(g1, g2, g3), each = R)
  
  term.g0 <- - sum(ms.prob[1,-c((R+1):(2*R))] * log(gamma[-c((R+1):(2*R))])) -
    sum(ms.prob[1,1:R] * log(fdts1[id,])) - 
    sum(ms.prob[1,(R+1):(2*R)] * log(fdts2[id,])) - 
    sum(ms.prob[1,(2*R+1):(3*R)] * log(fdts3[id,]))
  
  # Alpha 0
  Pt0     <- do.call(rbind, replicate(nA, ms.prob[,1:R], simplify = FALSE))
  term.a0 <- - sum(Pt0[-id,] * (-bt[-id,] + dt[-id,]*log(bt[-id,]) - lgamma(dt[-id,] + 1)))
  
  # Alpha 1
  term.a1 <- obj.alpha.1(alpha.1, ms.prob[-1,(R+1):(2*R)], xmat1 = x1[-id,,],
                         dtr = dt[-id,], btr = bt[-id,])
  
  # Alpha 2
  term.a2 <- obj.alpha.2(alpha.2, ms.prob[-1,(2*R+1):(3*R)], xmat2 = x2[-id,,],
                         dtr = dt[-id,], btr = bt[-id,])
  
  # Beta
  term.b <- obj.beta(beta = c(beta.0, beta.11, beta.22), jsprobs = js.prob, U = U,
                     xmat01 = xx01, xmat02 = xx02, xmat11 = xx11, xmat22 = xx22, 
                     dtr = dt, btr = bt, pena01 = pen01, pena02 = 0,
                     pena11 = 0, pena22 = 0)
  
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

# Create index list (see loop) 
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
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% alpha.1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% alpha.2)
  fdts1 <- matrix(dpois(dt, lambda = bt), nrow = nT*nA, ncol = R)
  fdts2 <- matrix(dpois(dt, lambda = bt * exp(x1a1)), nrow = nT*nA, ncol = R)
  fdts3 <- matrix(dpois(dt, lambda = bt * exp(x2a2)), nrow = nT*nA, ncol = R)
  fdts  <- cbind(fdts1, fdts2, fdts3)
  
  # Change 0 values by very small positive value to avoid numerical issues
  fdts[which(fdts == 0)] <- 10^(-30)
  
  # Transition model matrices
  x1b01 <- sapply(1:dim(xx01)[3], function(j) xx01[,,j] %*% beta.01 + U.new[j])
  x2b02 <- sapply(1:dim(xx02)[3], function(j) xx02[,,j] %*% beta.02 + U.new[j])
  x1b11 <- sapply(1:dim(xx11)[3], function(j) xx11[,,j] %*% beta.11 + U.new[j])
  x2b22 <- sapply(1:dim(xx22)[3], function(j) xx22[,,j] %*% beta.22 + U.new[j])
  
  # Calculate the transition probabilities
  pen01  <- pen01[1:nT,]
  pen02  <- x2b02*0
  pen11  <- x1b11*0
  pen22  <- x2b22*0
  pt.01 <- exp(pen01 + x1b01)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
  pt.02 <- exp(pen02 + x2b02)/(1 + exp(pen01 + x1b01) + exp(pen02 + x2b02))
  pt.00 <- 1 - pt.01 - pt.02
  
  pt.11 <- exp(pen11 + x1b11)/(1 + exp(pen11 + x1b11))
  pt.12 <- matrix(0, nrow = nT, ncol = R)
  pt.10 <- 1 - pt.11
  
  pt.21 <- matrix(0, nrow = nT, ncol = R)
  pt.22 <- exp(pen22 + x2b22)/(1 + exp(pen22 + x2b22))
  pt.20 <- 1 - pt.22
  
  P <- unname(cbind(pt.00, pt.01, pt.02, pt.10, pt.11, 
                    pt.12, pt.20, pt.21, pt.22))
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
    b[b == 0] <- 10^(-30)
    c <- a/rep(b, times = 9)
    
    if(any(is.na(c))) stop('na')
    
    lst[[t]][['a']][1:(9*R)] <- a
    lst[[t]][['b']][1:R]     <- b
    lst[[t]][['c']][1:(9*R)] <- c
  }
  
  # From joint to marginal
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
  term.a0 <- - sum(Pt0[-id,] * (-bt[-id,] + dt[-id,]*log(bt[-id,]) - lgamma(dt[-id,] + 1)))
  
  # Gamma
  fit0 <- optim(par = gamma[1], fn = obj.gamma, msprobs = ms.prob[1,],
                dens1 = fdts1, dens2 = fdts2, dens3 = fdts3,
                control = list(reltol = 10^(-16)), method = 'Brent',
                lower = 0.001, upper = 0.999)
  g1 <- fit0$par
  g2 <- 0 
  g3 <- 1 - g1 - g2
  gamma <- rep(c(g1, g2, g3), each = R)
  term.g0 <- - sum(ms.prob[1,-c((R+1):(2*R))] * log(gamma[-c((R+1):(2*R))])) -
    sum(ms.prob[1,1:R] * log(fdts1[id,])) - 
    sum(ms.prob[1,(R+1):(2*R)] * log(fdts2[id,])) - 
    sum(ms.prob[1,(2*R+1):(3*R)] * log(fdts3[id,]))
  
  # Alpha 1
  fit1 <- glm(form1.0,
              family = poisson(link = 'log'), data = df[-c(1 + nT*(0:(nA*R-1))),],
              start = alpha.1, control = list(epsilon = 10^(-16)),
              weights = as.vector(do.call('rbind', replicate(nA, ms.prob[-1,(R+1):(2*R)], 
                                                             simplify = FALSE))))
  alpha.1 <- fit1$coefficients
  
  # Alpha 2
  fit2 <- glm(form2.0,
              family = poisson(link = 'log'), data = df[-c(1 + nT*(0:(nA*R-1))),],
              start = alpha.2, control = list(epsilon = 10^(-16)),
              weights = as.vector(do.call('rbind', replicate(nA, ms.prob[-1,(2*R+1):(3*R)], 
                                                             simplify = FALSE))))
  alpha.2 <- fit2$coefficients
  
  # Beta's
  fit3 <- optim(par = c(beta.0, beta.11, beta.22), fn = obj.beta, U = U.new,
                jsprobs = js.prob, xmat01 = xx01, xmat02 = xx02, xmat11 = xx11, 
                xmat22 = xx22, dtr = dt, btr = bt, pena01 = pen01[1:nT,],
                pena02 = pen02[1:nT,], pena11 = pen11[1:nT,], pena22 = pen22[1:nT,],
                method = 'BFGS', control = list(reltol = 10^(-6)))
  
  beta.0  <- fit3$par[1:length(beta.0)]
  beta.01 <- fit3$par[1:length(beta.01)]
  beta.02 <- fit3$par[(length(beta.01)+1):(length(beta.01) + length(beta.02))]
  beta.11 <- fit3$par[(length(beta.01)+length(beta.02)+1):(length(beta.01) + length(beta.02) + length(beta.11))]
  beta.22 <- fit3$par[(length(beta.01)+length(beta.02)+length(beta.11)+1):(length(beta.01) + length(beta.02) + length(beta.11) + length(beta.22))]
  
  # Update spatial effect
  if(abs(Uopt(U.new[-R]) - ll0)/ll0 > 10^(88)){
    U.fit <- optim(par = U.new[-R], fn = Uopt, method = 'BFGS', 
                   control = list(reltol = 10^(-12)))
    U.new <- c(U.fit$par, -sum(U.fit$par))
    ll1   <- U.fit$value
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

# Save output + global environment
saveRDS(object = obj, file = paste0('Results/EMoutput_tau', tau, '.rds'))
save(list = ls(.GlobalEnv), file = paste0("Results/EMresult_tau", tau, ".Rdata"))

##### 7) Visualisations ----

# Load global environment
load(paste0("Results/EMresult_tau", tau, ".Rdata"))

# Summarize results
tab1 <- data.frame()
tab1 <- data.frame('Date'    = df$Date, 
                   'Deaths'  = df$Deaths,
                   'bDeaths' = df$bDeaths,
                   'Region'  = df$Region,
                   'Age'     = df$Age0,
                   'Prob0'   = as.vector(apply(ms.prob[,1:R], 1, mean)),
                   'Prob1'   = as.vector(apply(ms.prob[,(R+1):(2*R)], 1, mean)), 
                   'Prob2'   = as.vector(apply(ms.prob[,(2*R+1):(3*R)], 1, mean)),
                   'Tind95'  = df$w_avg_Tind95, 
                   'Flu'     = (df$ia100 - min(df$ia100))/(max(df$ia100) - min(df$ia100)),
                   'State'   = as.vector(sapply(1:R, function(j) 
                     apply(ms.prob[,c(j,j+R,j+2*R)],1, which.max) - 1)))

### Plot death counts  + baseline trend
p1 <- ggplot(tab1) + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Date, y = Deaths), col = 'black', linewidth = 1, alpha = 0.1) + 
  geom_line(aes(x = Date, y = bDeaths), col = ared, linewidth = 0.9)

ggsave(plot = p1,
       filename = paste0('C:/Users/u0131219/OneDrive - KU Leuven/Documents/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/weeklydeaths.pdf'),
       width = 6, height = 4)  

### Plot conditional probabilities
p21 <- ggplot(tab1) + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Date, y = Prob1, col = ablue), linewidth = 0.7, alpha = 0.75) + 
  geom_line(aes(x = Date, y = Tind95, col = ared), linewidth = 0.7) + 
  ylab('Smoothed marginal state prob.') + 
  scale_color_manual(values = c(ablue, ared), name = '',
                     labels = c(bquote('P(S'[t]~'=1 | D'['t']~', X'['t']~')'),'w_avg_Tind95')) + 
  theme(legend.position = 'bottom')

p22 <- ggplot(tab1) + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Date, y = Prob2, col = ablue), linewidth = 0.7, alpha = 0.75) + 
  geom_line(aes(x = Date, y = Flu, col = ared), linewidth = 0.7) + 
  ylab('Smoothed marginal state prob.') + 
  scale_color_manual(values = c(ablue, ared), name = '',
                     labels = c(bquote('P(S'[t]~'=2 | D'['t']~', X'['t']~')'),'Flu')) + 
  theme(legend.position = 'bottom')

p2 <- ggpubr::ggarrange(p21, p22, nrow = 1)

ggsave(plot = p2,
       filename = paste0('C:/Users/u0131219/OneDrive - KU Leuven/Documents/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/mstateprob.pdf'),
       width = 8, height = 4)  

### States on death count time series
tab1r <- tab1 %>% dplyr::filter(Region == regions[19] & Age == 'Y_GE90')
p3 <- ggplot(tab1r) + 
  theme_bw(base_size = 15) + ylab('Death counts') + 
  geom_rect(aes(xmin = Date - 3.5, xmax = Date + 3.5, ymin = min(Deaths)-100, 
                ymax = max(Deaths) + 100, fill = as.factor(State)), alpha = 0.4) + 
  scale_fill_manual(values = brewer.pal(3,'Dark2'), name = 'State') +
  geom_line(aes(x = Date, y = bDeaths), col = 'red') + 
  geom_line(aes(x = Date, y = Deaths), col = 'gray20') +
  theme(legend.position = 'bottom') + 
  coord_cartesian(ylim = c(min(tab1r$Deaths), max(tab1r$Deaths)))
p3

ggsave(plot = p3,
       filename = paste0('C:/Users/u0131219/OneDrive - KU Leuven/Documents/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/mstateprobdeaths.pdf'),
       width = 7, height = 4.5)  

