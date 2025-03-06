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
              'progress', 'mgcv')
invisible(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load
load("Results/EMresult_tau10.RDATA")

# rm(list=setdiff(ls(), c('obj','R','nA','nT','dt','bt','df','pen1','pen2',
#                         'x1','x2','xx1','xx2', 'indl0', 'indl1', 'indl2',
#                         'fit1','fit2', 'Q', 'cov', 'pt.00', 'pt.01', 'pt.02',
#                         'pt.10', 'pt.11', 'pt.20', 'pt.22', 'obj.alpha.1',
#                         'obj.alpha.2','obj.beta', 'ms.prob')))


##### 1) Fisher information -----

### Object
obj[[length(obj)]] <- NULL
iteropt <- length(obj)

### Iteration results
alpha1.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.1))
alpha2.iter <- do.call('rbind', lapply(obj, function(x) x$alpha.2))
beta0.iter  <- do.call('rbind', lapply(obj, function(x) x$beta.0))
beta11.iter <- do.call('rbind', lapply(obj, function(x) x$beta.11))
beta22.iter <- do.call('rbind', lapply(obj, function(x) x$beta.22))
gamma.iter  <- do.call('rbind', lapply(obj, function(x) x$gamma))
U.iter      <- do.call('rbind', lapply(obj, function(x) x$U)) 

theta.iter <- cbind(alpha1.iter, alpha2.iter, beta0.iter, 
                    beta11.iter, beta22.iter, gamma.iter)

### Optimal parameter vector

# Length vector
d <- ncol(alpha1.iter) + ncol(alpha2.iter) + ncol(beta0.iter) + 
  ncol(beta11.iter) + ncol(beta22.iter) + ncol(gamma.iter) - 2

# Optimal parameter vector
alpha1.opt <- alpha1.iter[iteropt,]
alpha2.opt <- alpha2.iter[iteropt,]
beta0.opt  <- beta0.iter[iteropt,]
beta01.opt <- beta0.opt[1:ncol(xx01)]
beta02.opt <- beta0.opt[-c(1:ncol(xx01))]
beta11.opt <- beta11.iter[iteropt,]
beta22.opt <- beta22.iter[iteropt,]
gamma.opt  <- gamma.iter[iteropt,1] # g2 set at zero (g3 st g1 + g3 = 0)
U.opt      <- U.iter[iteropt,]

alpha.opt  <- c(alpha1.opt, alpha2.opt)
beta.opt   <- c(beta0.opt, beta11.opt, beta22.opt)
theta.opt  <- unname(c(alpha1.opt, alpha2.opt, beta0.opt, beta11.opt,
                       beta22.opt))

##### Results SE's
B   <- 25000
obj <- readRDS(file = 'Results/FI.rds')

FI.alpha <- Reduce('+',lapply(obj, function(x) x$alpha))/B
FI.beta  <- Reduce('+',lapply(obj, function(x) x$beta))/B

COV.alpha <- ginv(FI.alpha)
COV.beta  <- ginv(FI.beta)

SD.alpha <- sqrt(diag(COV.alpha))
SD.beta  <- sqrt(diag(COV.beta))

SD.alpha
SD.beta

##### 2) Predicted base-line fit  ----

### Calibration and forecasting period
t0 <- as.Date('2012-12-31', format = '%Y-%m-%d')
t1 <- as.Date('2024-07-01', format = '%Y-%m-%d') # End of 26-th ISO-week
t2 <- as.Date('2029-12-30', format = '%Y-%m-%d')

t01 <- seq(t0, t1 - 1, by="days")
t12 <- seq(t1, t2, by="days")
df.t01 <- data.frame('Date' = t01, 'ISOWeek' = isoweek(t01), 
                     'ISOYear' = isoyear(t01))
df.t12 <- data.frame('Date' = t12, 'ISOWeek' = isoweek(t12), 
                     'ISOYear' = isoyear(t12))

### Prepare dataset for forecasting
ages <- sort(unique(df$Age0))

dff <- df.t12 %>% 
  dplyr::select(c(ISOYear, ISOWeek)) %>% 
  unique() %>%
  mutate('Date' = ISOweek::ISOweek2date(sprintf("%d-W%02d-%d", ISOYear, ISOWeek, 1)))

dff <- dff %>%
  mutate('fsin52' = sin(2*pi*ISOWeek/52.1775),
         'fcos52' = cos(2*pi*ISOWeek/52.1775),
         'fsin26' = sin(4*pi*ISOWeek/52.1775),
         'fcos26' = cos(4*pi*ISOWeek/52.1775))
dff <- do.call('rbind', lapply(1:length(ages),
                               function(a) data.frame(dff, 'Age' = ages[a])))
dff <- do.call('rbind', lapply(1:length(regions), 
                               function(r) data.frame(dff, 'Region' = regions[r])))
dff$Region <- factor(dff$Region, levels = sort(unique(dff$Region)))
dff$Age    <- factor(dff$Age, levels = sort(unique(dff$Age)))
dff$Time   <- dff$ISOYear  - isoyear(t0)

dff <- dff %>%
  dplyr::arrange(Region, Age, Date)

# Formulas
formula0 <- ~ -1 + Region:Age
formula1 <- ~ -1 + Time:Region:Age
formula2 <- ~ -1 + fsin52:Region:Age
formula3 <- ~ -1 + fcos52:Region:Age
formula4 <- ~ -1 + fsin26:Region:Age
formula5 <- ~ -1 + fcos26:Region:Age

# Model matrices
X0 <- model.matrix(formula0, data = dff)
X1 <- model.matrix(formula1, data = dff)
X2 <- model.matrix(formula2, data = dff)
X3 <- model.matrix(formula3, data = dff)
X4 <- model.matrix(formula4, data = dff)
X5 <- model.matrix(formula5, data = dff)

# Adjacency matrices
adjmat <- nb2mat(poly2nb(shapef[which(shapef$NUTS_ID %in% regions),'geometry']), 
                 zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat) <- list(shapef[which(shapef$NUTS_ID %in% regions),]$NUTS_ID, 
                         shapef[which(shapef$NUTS_ID %in% regions),]$NUTS_ID)
diag(adjmat) <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))

# With penalty
nR  <- nrow(adjmat)
nA  <- length(unique(dff$Age))
ADJ <- matrix(0, nrow = nR*nA, ncol = nR*nA)
for(a in 1:nA){
  ADJ[((a-1)*nR + 1):(a*nR),((a-1)*nR + 1):(a*nR)] <- adjmat
}

# Add model matrices
dff$M0 <- X0
dff$M1 <- X1
dff$M2 <- X2
dff$M3 <- X3
dff$M4 <- X4
dff$M5 <- X5

# Baseline fit
fitb.train  <- readRDS('Results/GAM.rds')
dgamb.train <- readRDS('Results/df_GAM.rds') %>%
  dplyr::filter(Date < as.Date('2024-07-01'))

# Forecasting dataset
dff <- dff %>% 
  mutate('Expo' = 1000)
dgamb.forecast <- dff %>%
  dplyr::select(c('Expo','M0','M1','M2','M3','M4','M5'))

# Baseline death counts in matrix
nT.forecast <- length(unique(dff$Date))
bt.forecast <- matrix(predict(fitb.train, newdata = dgamb.forecast, 
                              offset = dgamb.forecast$Expo, type = 'response'),
                      nrow = nT.forecast*nA, ncol = R) * dgamb.forecast$Expo

# Extract parameter estimates
Mu     <- fitb.train$coefficients %>% unname()
Sigma  <- vcov(fitb.train) %>% unname()
Xd     <- dgamb.forecast %>% dplyr::select(-c('Expo')) %>% SparseM::as.matrix()

# Remove some things 
rm(fitb.train)
gc()

##### 3) Using scenarios for the covariate information ----

### Temperature anomalies and extreme indicators for forecasting period ----

# Left join temperature data
df01 <- df.t01 %>% 
  left_join(readRDS(file = 'Data/E-OBS/FR_NUTS2_tg_daily.rds'), by = 'Date') %>%
  mutate(Method = 'Obs')  %>%
  mutate(Day = yday(Date), .before = 'ISOWeek') 
df12.rcp26 <- df.t12 %>% 
  left_join(readRDS(file = 'Data/CDS/FR_NUTS2_tg_daily_rcp26_2024_2030.rds'), by = c('Date')) %>%
  mutate(Method = 'RCP 2.6')
df12.rcp45 <- df.t12 %>%
  left_join(readRDS(file = 'Data/CDS/FR_NUTS2_tg_daily_rcp45_2024_2030.rds'), by = c('Date')) %>%
  mutate(Method = 'RCP 4.5')
df12.rcp85 <- df.t12 %>%
  left_join(readRDS(file = 'Data/CDS/FR_NUTS2_tg_daily_rcp85_2024_2030.rds'), by = c('Date')) %>%
  mutate(Method = 'RCP 8.5')
df12 <- do.call('rbind', list(df12.rcp26, df12.rcp45, df12.rcp85)) %>%
  mutate(Day = yday(Date), .before = 'ISOWeek') 

# Plot
plot_rcp <- function(r){
  df02  <- rbind(df01, df12) 
  df02r <- df02 %>% dplyr::filter(Region == r)
  df02r$Obs <- as.character((df02r$Date <= max(df01$Date)))
  
  ggplot(df02r) + 
    xlab('Date') + ylab('Daily average temperature') + 
    ggtitle(shapef$NUTS_NAME[shapef$NUTS_ID == r]) + 
    theme_bw(base_size = 15) +
    geom_line(aes(x = Date, y = tg, col = Method, group = Method, linetype = Method)) + 
    scale_color_manual(name = '', 
                       values = c('Obs' = 'black', 'RCP 2.6' = agreen,
                                  'RCP 4.5' = ablue, 'RCP 8.5' = ared),
                       labels = c('Observed', 'RCP 2.6', 'RCP 4.5', 'RCP 8.5')) + 
    scale_linetype_manual(name = '', 
                          values = c('Obs' = 1, 'RCP 2.6' = 1,
                                     'RCP 4.5' = 2, 'RCP 8.5' = 3),
                          labels = c('Observed', 'RCP 2.6', 'RCP 4.5', 'RCP 8.5')) + 
    guides(color = guide_legend(override.aes = list(linewidth = 1.25, size = 10))) +
    theme(legend.position = 'bottom')
}

# Function to construct anomalies
excess_temp <- function(r){
  # Data from region r
  df01r <- df01 %>% dplyr::filter(Region == r)
  df12r <- df12 %>% dplyr::filter(Region == r)
  
  # Baseline fits
  fitr.tg <- robust::lmRob(tg ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25),
                           data = df01r)
  # Predictions
  pred.tg <- predict(fitr.tg, df12r)
  
  # Quantiles
  q01 <- quantile(df01r$tg, c(0.05, 0.95), na.rm = TRUE)
  
  # Information
  dd <- data.frame(Date    = df12r$Date, 
                   Region  = df12r$Region,
                   tg_anom = df12r$tg - pred.tg,
                   Tind95  = (df12r$tg > q01[2])*1,
                   Tind5   = (df12r$tg < q01[1])*1,
                   Method  = df12r$Method)
  dd
}

# Apply function
list.exc.t <- lapply(sort(unique(df12$Region)), excess_temp)
df.exc.t   <- do.call('rbind', list.exc.t)
df.t       <- df12 %>% left_join(df.exc.t, by = c('Date','Region', 'Method'))

# Transform to weekly - on ISO week basis
df.t.w <- df.t %>% group_by(ISOYear, ISOWeek, Region, Method) %>%
  summarize(w_avg_tg      = mean(tg), 
            w_avg_tg_anom = mean(tg_anom), 
            w_avg_Tind95  = mean(Tind95), 
            w_avg_Tind5   = mean(Tind5)) %>%
  ungroup() 


### Respiratory activity scenarios  ----

### Baseline death counts
df.bDeaths01 <- readRDS('GAM_Results_bDeathsx.rds')
df.bDeaths12 <- dff %>% 
  dplyr::select(c('Date','ISOYear','ISOWeek','Region','Age','Expo')) %>%
  mutate(bDeaths = as.vector(bt.forecast))


### Influenza and ARI
df.epi01 <- readRDS('df_NUTS2_FR.rds') %>%
  dplyr::select(c('Date','ISOYear', 'ISOWeek', 'Region','ia100')) %>%
  mutate(Type = 'Obs')

# Influenza
list.ia_scenario <- lapply(1:length(regions), function(r) {
  # Read JAGS model
  if1 <- readRDS(paste0('SIRS models/SIR_', regions[r], '.rds'))
  
  # Extract posterior samples
  I_pred_samples <- if1$BUGSoutput$sims.list$I_pred
  
  # Compute quantiles
  l1 <- sapply(1:ncol(I_pred_samples), function(x) quantile(I_pred_samples[,x], 0.50))
  l2 <- sapply(1:ncol(I_pred_samples), function(x) quantile(I_pred_samples[,x], 0.75))
  l3 <- sapply(1:ncol(I_pred_samples), function(x) quantile(I_pred_samples[,x], 0.95))
  
  # Make data frame to summarize results
  time.proj   <- seq(max(df.epi01$Date) + 7, 
                     max(df.epi01$Date) + 7*52*5 + 7*27, by = 'week')
  
  dd <- data.frame('Date'    = rep(time.proj, times = 3),
                   'ISOYear' = rep(isoyear(time.proj), times = 3),
                   'ISOWeek' = rep(isoweek(time.proj), times = 3),
                   'Region'  = regions[r], 
                   'ia100'   = c(tail(l1, nT.forecast), tail(l2, nT.forecast),
                                 tail(l3, nT.forecast)),
                   'Type' = rep(c('M', 'H', 'S'), each = nT.forecast))
  dd
})
df.ia12 <- do.call('rbind', list.ia_scenario)

# Combine
df.epi12 <- df.ia12 

# Plot
plot_ia <- function(r){
  dd0 <- df.epi01 %>% dplyr::filter(Region == r)
  dd1 <- df.epi12 %>% dplyr::filter(Region == r)
  
  dd01 <- do.call('rbind', list(dd0, dd1))
  dd01$Type <- factor(dd01$Type, levels = c('Obs', 'S', 'H', 'M'))
  
  ggplot(dd01) + 
    xlab('Date') + ylab('IA incidence (per 1 000 inh.)') + 
    ggtitle(shapef$NUTS_NAME[shapef$NUTS_ID == r]) + 
    theme_bw(base_size = 15) +
    geom_line(aes(x = Date, y = ia100/1000, col = Type, group = Type, linetype = Type)) +
    scale_color_manual(name = '', 
                       labels = c('Obs' = 'Observed', 'S' = 'Severe', 'H' = 'High', 'M' = 'Moderate'),
                       values = c('Obs' = 'black', 'S' = ared, 'H' = ablue, 'M' = agreen),
                       breaks = c('Obs', 'M', 'H', 'S')) +
    scale_linetype_manual(name = '',
                          values = c('Obs' = 1, 'S' = 3, 'H' = 2, 'M' = 1),
                          breaks = c('Obs', 'M', 'H', 'S'),
                          labels = c('Obs' = 'Observed', 'S' = 'Severe', 'H' = 'High', 'M' = 'Moderate')) + 
    guides(color = guide_legend(override.aes = list(linewidth = 1.25, size = 10))) +
    theme(legend.position = 'bottom') 
  
}

reg1 <- regions[12]
reg2 <- regions[21]
p.rcp <- plot_rcp(r = reg1) + 
  theme(legend.position = "none") + 
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
q.rcp <- plot_rcp(r = reg2)
p.ia  <- plot_ia(r = reg1) + 
  theme(legend.position = "none") + 
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
q.ia <- plot_ia(r = reg2) 

pscenario <- ggpubr::ggarrange(plotlist = list(p.rcp, p.ia, q.rcp, q.ia),
                               align = 'hv', nrow = 2, ncol = 2, vjust = 1000)
# ggsave(plot = pscenario,
#        filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/scenarios.pdf'),
#        width = 15, height = 11)

# for(j in 1:length(regions)){
#   p.rcp <- plot_rcp(r = regions[j])
#   q.ia  <- plot_ia(r  = regions[j])
# 
#   pq.rcp.ia <- ggpubr::ggarrange(plotlist = list(p.rcp, q.ia),
#                                  align = 'hv', nrow = 1, ncol = 2)
# 
#   ggsave(plot = pq.rcp.ia,
#          filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/',
#                            'Reserving phd/Papers/Paper 5 - climate and influenza shocks/',
#                            'Figures/Scenarios/', regions[j], '.pdf'),
#          width = 12, height = 5)
# }

### Hospital admissions
df.hosp01 <- readRDS('Hosp_NUTS2.rds') %>%
  dplyr::filter(Date <= t1 - 1, Date >= t0) 
df.epi12$Nhosp <- 0

### Final data-set
df <- df.bDeaths12 %>%
  dplyr::left_join(df.t.w, by = c('ISOYear', 'ISOWeek', 'Region'),
                   relationship = 'many-to-many') %>%
  dplyr::left_join(df.epi12,  by = c('Date', 'ISOYear', 'ISOWeek', 'Region'),
                   relationship = 'many-to-many') 

# Per 1000 inhabitants
df.epi01$ia100  <- df.epi01$ia100/1000
df$ia100        <- df$ia100/1000

# Seasonal
ia.seasonal <- df.epi01 %>%
  dplyr::group_by(Region, ISOWeek) %>%
  summarize('avg_ia100' = mean(ia100)) 

df <- df %>%
  dplyr::left_join(ia.seasonal, by = c('ISOWeek', 'Region')) %>%
  mutate(ia100 = ia100 - avg_ia100)

# Quantile exceedances
qvec <- quantile(df.epi01$ia100,  0.75)

# Quantiles
df$ia100q  <- pmax(df$ia100  - qvec, 0)
df$Nhospq  <- 0

# Lags
vars <- c("Date", "ISOYear", "ISOWeek", "Region",  "Age", "Expo", "bDeaths", 
          "w_avg_tg", "w_avg_tg_anom", "w_avg_Tind95", "w_avg_Tind5", "ia100", 
          "Nhosp", "ia100q", "Nhospq", 'Type', 'Method')

df <- df %>% 
  dplyr::select(all_of(vars))

# Select Method and Type
TYPE   <- 'M'
METHOD <- 'RCP 2.6'
df     <- df %>% 
  dplyr::filter(Type == TYPE, Method == METHOD) %>%
  dplyr::select(-c('Type', 'Method'))

# Extend with last 3 weeks (to compute lags)
df.prev <- readRDS('dfRS.rds') %>%
  dplyr::select(-c('Deaths', 'geo_name','Region.y.y','geometry', 'avg_ia100', 
                   'ari100', 'ari100q')) %>%
  dplyr::filter(ISOYear == year(t1) & ISOWeek %in% 24:26)

df <- rbind(df.prev, df) %>%
  arrange(Region, Age, Date)

vars <- c("w_avg_tg_anom", "w_avg_Tind95", "w_avg_Tind5", "ia100", 
          "Nhosp", "ia100q", "Nhospq")

### Compute lags ----
for(v in vars){
  df <- df %>% 
    group_by(Region, Age) %>%
    mutate(!!paste0(v,'_l1')   := lag(!!sym(v), 1),
           !!paste0(v,'_l2')   := lag(!!sym(v), 2),
           !!paste0(v,'_l0.1') := (lag(!!sym(v), 0) + lag(!!sym(v), 1))/2,
           !!paste0(v,'_l2.3') := (lag(!!sym(v), 2) + lag(!!sym(v), 3))/2) %>% 
    ungroup()
}

### Transition restrictions ----

# States 1-2
df$penalty01 <- ifelse(df$w_avg_Tind95 > 0, 0, -1000)
df$penalty02 <- 0

### Remove missing observations
df <- df %>% na.omit()

### Arrange data frame by (Region, Time)
df <- df %>% 
  dplyr::arrange(Region, Age, ISOYear, ISOWeek) %>%
  mutate(Region = as.factor(Region),
         Time = ISOYear - min(ISOYear))

# In matrix form
pen01  <- matrix(df$penalty01, nrow = nT.forecast*nA, ncol = R, dimnames = list(NULL, regions))
pen02  <- matrix(df$penalty02, nrow = nT.forecast*nA, ncol = R, dimnames = list(NULL, regions))
pen01  <- pen01[1:nT.forecast,]
pen02  <- pen02[1:nT.forecast,]
pen11  <- pen01*0
pen22  <- pen02*0

### New design matrices ----

# Age grouping
df$Age0 <- df$Age
df$Age  <- plyr::mapvalues(df$Age0, from = c('Y65-69','Y70-74','Y75-79',
                                             'Y80-84','Y85-89','Y_GE90'),
                           to = c('Y65-74','Y65-74','Y75-84','Y75-84',
                                  'Y_GE85','Y_GE85'))
df$Age <- factor(df$Age, levels = c('Y65-74','Y75-84','Y_GE85'))

# State 1
form1.0 <- as.formula( ~ -1 + offset(log(bDeaths)) + 
                         w_avg_Tind95:Age +
                         w_avg_Tind95_l1:Age + 
                         w_avg_Tind95_l2:Age +                       
                         w_avg_tg_anom:Age +
                         w_avg_tg_anom_l1:Age + 
                         w_avg_tg_anom_l2:Age)
m10     <- model.matrix(form1.0, data = df)
x1      <- simplify2array(lapply(1:R, function(j) 
  m10[(nA*nT.forecast*(j-1) + 1):(nA*nT.forecast*j),]))

# State 2
form2.0 <- as.formula( ~ -1 + offset(log(bDeaths)) +
                         w_avg_Tind5_l0.1:Age +
                         w_avg_Tind5_l2.3:Age + 
                         ia100q_l0.1:Age +
                         ia100q_l2.3:Age +
                         Nhospq_l0.1:Age + 
                         Nhospq_l2.3:Age)
m20     <- model.matrix(form2.0, data = df)
x2      <- simplify2array(lapply(1:R, function(j) 
  m20[(nA*nT.forecast*(j-1) + 1):(nA*nT.forecast*j),]))

# Trans probs
m01 <- model.matrix( ~ offset(log(bDeaths)) + w_avg_Tind95,
                     data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx01 <- simplify2array(lapply(1:R, function(j) 
  matrix(m01[(nT.forecast*(j-1) + 1):(nT.forecast*j),], ncol = ncol(m01))))

m11 <- model.matrix( ~ offset(log(bDeaths)) + 
                       w_avg_Tind95 + w_avg_Tind95_l1 + w_avg_Tind95_l2,
                     data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx11 <- simplify2array(lapply(1:R, function(j) 
  m11[(nT.forecast*(j-1) + 1):(nT.forecast*j),]))

m02 <- model.matrix( ~ offset(log(bDeaths)) + 
                       ia100q_l0.1 + Nhospq_l0.1,
                     data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx02 <- simplify2array(lapply(1:R, function(j) 
  m02[(nT.forecast*(j-1) + 1):(nT.forecast*j),]))

m22 <- model.matrix( ~ offset(log(bDeaths)) + 
                       ia100q_l0.1 + Nhospq_l0.1 +
                       ia100q_l2.3 + Nhospq_l2.3,
                     data = df %>% dplyr::filter(Age0 == 'Y_GE90'))
xx22 <- simplify2array(lapply(1:R, function(j) 
  m22[(nT.forecast*(j-1) + 1):(nT.forecast*j),]))

### CI death counts: param unc (base + RS param) + RS state unc + Poisson noise ----

B <- 25000

cl <- snow::makeCluster(detectCores() - 1) 
registerDoSNOW(cl)
snow::clusterEvalQ(cl, c(library(dplyr), library(MASS)))
snow::clusterExport(cl, list = c('Mu','Sigma','Xd'))   

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
  Dxtb  <- exp(Xd %*% t(Mub))*1000
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
  
  # matrices
  x1a1  <- sapply(1:dim(x1)[3], function(j) x1[,,j] %*% a1)
  x2a2  <- sapply(1:dim(x2)[3], function(j) x2[,,j] %*% a2)
  x1b01 <- pen01 + sapply(1:R,  function(r) xx01[,,r] %*% b01 + Usim[r])
  x2b02 <- pen02 + sapply(1:R,  function(r) xx02[,,r] %*% b02 + Usim[r]) 
  x1b11 <- pen11 + sapply(1:R,  function(r) xx11[,,r] %*% b11 + Usim[r])
  x2b22 <- pen22 + sapply(1:R,  function(r) xx22[,,r] %*% b22 + Usim[r]) 
  
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
  St <- matrix(NA, nrow = nT.forecast, ncol = R)
  St[1,] <- 1
  
  # Sample paths regime states
  for(t in 2:nT.forecast){
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
  
  Dt  <- matrix(NA, nrow = nT.forecast*nA, ncol = R)
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

df.p1 <- cbind(df, qmat)

# saveRDS(df.p1, file = 'Forecast_RCP26v2_M.rds')
df.pL <- readRDS('Forecast_RCP26v2_M.rds') %>%
  mutate('Scenario' = 'Low')
df.pM <- readRDS('Forecast_RCP45v2_H.rds') %>%
  mutate('Scenario' = 'Medium')
df.pH <- readRDS('Forecast_RCP85v2_S.rds') %>%
  mutate('Scenario' = 'High')
df.p1 <- do.call(rbind, list(df.pL, df.pM, df.pH))
df.p1$Scenario <- factor(df.p1$Scenario, levels = c('High', 'Medium', 'Low'))

df.p1 <- df.p1 %>%
  dplyr::select(c(Date, ISOYear, ISOWeek, Region, Age0, 
                  `2.5%`, `50%`, `97.5%`, Scenario)) %>%
  dplyr::mutate(Deaths = NA,
                bDeaths = NA,
                Expo = NA,
                bDeaths1000 = NA,
                Deaths1000 = NA, .after = 'Age0')
colnames(df.p1) <- c('Date', 'ISOYear', 'ISOWeek', 'Region', 'Age', 'Deaths', 
                     'bDeaths', 'Expo', 'bDeaths1000', 'Deaths1000', 
                     '2.5%', '50%', '97.5%', 'Scenario')

# Baseline death counts (forecasted)
fitb.train  <- readRDS('GAMx.rds')
pred <- predict(fitb.train, newdata = dgamb.forecast, 
                offset = dgamb.forecast$Expo, type = 'response')* dgamb.forecast$Expo
dff$bDeathsPred <- pred
dff.add <- dff %>% 
  dplyr::select(c('Date', 'ISOYear', 'ISOWeek', 'Age', 'Region', 'bDeathsPred'))

df.p1 <- df.p1 %>%
  dplyr::left_join(dff.add, by = join_by(Date, ISOYear, ISOWeek, Region, Age)) %>%
  mutate(bDeaths1000 = bDeathsPred, .after = 'Expo') %>%
  dplyr::select(-c(bDeathsPred))

df.bDeaths01 <- df.bDeaths01 %>%
  dplyr::select(c(Date, ISOYear, ISOWeek, Region, Age, Deaths, bDeaths, Expo)) %>%
  mutate('2.5%' = NA,
         '50%' = NA, 
         '97.5%' = NA)

plot_scenmort <- function(region, age, agegr){
  dfplotf <- df.p1 %>% dplyr::filter(Age == age, Region == region)
  dfobsf  <- df.bDeaths01 %>% 
    dplyr::filter(Age == unique(dfplotf$Age),  
                  Region == unique(dfplotf$Region),
                  Date < as.Date('2024-07-01')) %>%
    mutate(bDeaths1000 = bDeaths/Expo * 1000, 
           Deaths1000 = Deaths/Expo * 1000, .after = 'Expo') %>%
    mutate(Scenario = 'Observed')
  dfra   <- do.call('rbind', list(dfplotf, dfobsf))
  
  
  p1 <- ggplot(dfra) + 
    theme_bw(base_size = 15) + xlab('Date') + ylab('Death rate') + 
    geom_line(data = dfra %>% dplyr::filter(Date < as.Date('2024-07-01')) %>%
                mutate(Scenario = factor(Scenario, levels = c('Observed'))), 
              aes(x = Date, y = Deaths1000, col = Scenario)) +
    geom_ribbon(data = dfra %>% dplyr::filter(Date >= as.Date('2024-07-01')),
                aes(x = Date, ymin = `2.5%`, ymax  = `97.5%`, group = Scenario, 
                    fill = Scenario), alpha = 0.75) +
    scale_fill_manual(name = '', values = c('Low' = agreen, 'Medium' = ablue, 'High' = ared),
                      labels = c('Low' = 'Scenario 1', 'Medium' = 'Scenario 2', 'High' = 'Scenario 3'),
                      breaks = c('Low','Medium', 'High')) +
    scale_color_manual(name = '', values = c('Observed' = 'black')) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.25, size = 10)),
           fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.position = 'bottom') + 
    ggtitle(paste0(shapef[shapef$NUTS_ID == dfra$Region[1],'NUTS_NAME'],' - Age group ',agegr)) 
  
  p1
}

p1 <- plot_scenmort(region = regions[2],  age = 'Y_GE90', agegr = '90+')
p2 <- plot_scenmort(region = regions[12], age = 'Y_GE90', agegr = '90+')
p3 <- plot_scenmort(region = regions[17], age = 'Y_GE90', agegr = '90+')
p4 <- plot_scenmort(region = regions[21], age = 'Y_GE90', agegr = '90+')

p1234 <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4),
                           common.legend = TRUE, nrow = 2,ncol = 2,
                           align = 'hv', legend = 'bottom')

# ggsave(plot = p1234,
#        filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/scenpredddd.pdf'),
#        width = 14, height = 10)

# for(j in 1:R){
#   df.p1r  <- df.p1 %>% dplyr::filter(Region == regions)
#   bounds  <- range(c(df.p1r$`2.5%`, df.p1r$`97.5%`))
# 
#   a1 <- plot_scenmort(regions[j], 'Y65-69', '65-69')
#   a2 <- plot_scenmort(regions[j], 'Y70-74', '70-74')
#   a3 <- plot_scenmort(regions[j], 'Y75-79', '75-79')
#   a4 <- plot_scenmort(regions[j], 'Y80-84', '80-84')
#   a5 <- plot_scenmort(regions[j], 'Y85-89', '85-89')
#   a6 <- plot_scenmort(regions[j], 'Y_GE90', '90+')
# 
#   a123456 <- ggpubr::ggarrange(plotlist = list(a1, a2, a3, a4, a5, a6),
#                                common.legend = TRUE, align = 'hv',
#                                nrow = 2, ncol = 3, legend = 'bottom')
# 
#   ggsave(plot = a123456,
#          filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/',
#                            'Reserving phd/Papers/Paper 5 - climate and influenza shocks',
#                            '/Figures/MortalityForecasts/', regions[j], '.pdf'),
#          width = 14, height = 10)
# 
# }


### Table for excess deaths

# Baseline death counts (forecasted)
fitb.train  <- readRDS('GAMx.rds')
pred <- predict(fitb.train, newdata = dgamb.forecast, 
                offset = dgamb.forecast$Expo, type = 'response')* dgamb.forecast$Expo
dff$bDeathsPred <- pred
dff.add <- dff %>% 
  dplyr::select(c('Date', 'ISOYear', 'ISOWeek', 'Age', 'Region', 'bDeathsPred'))

dt <- df.p1 %>%
  dplyr::left_join(dff.add, by = join_by(Date, ISOYear, ISOWeek, Region, Age)) %>%
  mutate(bDeaths1000 = bDeathsPred) %>%
  dplyr::select(-c(bDeathsPred))

summary <- dt %>%
  group_by(Region, Age, Scenario) %>%
  reframe('Excess_M' = c(sum(`50%`) - sum(bDeaths1000))/sum(bDeaths1000)) %>%
  dplyr::left_join(shapef, by = c('Region' = 'NUTS_ID'))

sL.90 <- summary %>% dplyr::filter(Age == 'Y_GE90', Scenario == 'Low')
sM.90 <- summary %>% dplyr::filter(Age == 'Y_GE90', Scenario == 'Medium')
sH.90 <- summary %>% dplyr::filter(Age == 'Y_GE90', Scenario == 'High')

pL <- ggplot() +
  geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
  geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
  geom_sf(data = sL.90$geometry, aes(fill = as.numeric(sL.90$Excess_M))) +
  theme_bw(base_size = 10) + xlab('Long') + 
  theme(legend.position = 'bottom') + 
  scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                       name = '') +
  geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
  geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
  ggtitle(paste0('Scenario 1'))
pM <- ggplot() +
  geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
  geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
  geom_sf(data = sM.90$geometry, aes(fill = as.numeric(sM.90$Excess_M))) +
  theme_bw(base_size = 10) + xlab('Long') + 
  theme(legend.position = 'bottom') + 
  scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                       name = '') +
  geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
  geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
  ggtitle(paste0('Scenario 2'))
pH <- ggplot() +
  geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
  geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
  geom_sf(data = sH.90$geometry, aes(fill = as.numeric(sH.90$Excess_M))) +
  theme_bw(base_size = 10) + xlab('Long') + 
  theme(legend.position = 'bottom') + 
  scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                       name = '') +
  geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
  geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
  ggtitle(paste0('Scenario 3'))

pLMH <- ggpubr::ggarrange(plotlist = list(pL, pM, pH), align = 'hv',
                          legend = 'bottom', nrow = 1)

# ggsave(plot = pLMH,
#        filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/Reserving phd/Papers/Paper 5 - climate and influenza shocks/Figures/scenED.pdf'),
#        width = 10, height = 4.75)

plot_excess_deaths <- function(age, agegr){
  sL.90 <- summary %>% dplyr::filter(Age == age, Scenario == 'Low')
  sM.90 <- summary %>% dplyr::filter(Age == age, Scenario == 'Medium')
  sH.90 <- summary %>% dplyr::filter(Age == age, Scenario == 'High')
  
  pL <- ggplot() +
    geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
    geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
    geom_sf(data = sL.90$geometry, aes(fill = as.numeric(sL.90$Excess_M))) +
    theme_bw(base_size = 10) + xlab('Long') + 
    theme(legend.position = 'bottom') + 
    scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                         name = '') +
    geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
    geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
    ggtitle(paste0('Scenario 1'))
  pM <- ggplot() +
    geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
    geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
    geom_sf(data = sM.90$geometry, aes(fill = as.numeric(sM.90$Excess_M))) +
    theme_bw(base_size = 10) + xlab('Long') + 
    theme(legend.position = 'bottom') + 
    scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                         name = '') +
    geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
    geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
    ggtitle(paste0('Scenario 2'))
  pH <- ggplot() +
    geom_sf(data = shapef$geometry, col = 'white', show.legend = FALSE) + 
    geom_sf(data = shapef %>% st_union(), col = 'black', size = 2) +
    geom_sf(data = sH.90$geometry, aes(fill = as.numeric(sH.90$Excess_M))) +
    theme_bw(base_size = 10) + xlab('Long') + 
    theme(legend.position = 'bottom') + 
    scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                         name = '') +
    geom_sf(data = shapef$geometry, col = 'white', fill = NA, show.legend = FALSE) + 
    geom_sf(data = shapef$geometry %>% st_union(), col = 'white', fill = NA, show.legend = FALSE, linewidth = 1) +
    ggtitle(paste0('Scenario 3'))
  
  pLMH <- ggpubr::ggarrange(plotlist = list(pL, pM, pH), align = 'hv',
                            legend = 'bottom', nrow = 1)
  
  pcomba <- ggpubr::annotate_figure(pLMH, top = 
                                      ggpubr::text_grob(paste0("     Age group ", agegr),
                                                        face = "bold", size = 20))
  
  ggsave(plot = pLMH,
         filename = paste0('C:/Users/u0131219/OneDrive - UvA/Documenten/', 
                           'Reserving phd/Papers/Paper 5 - climate and influenza shocks/',
                           'Figures/ExcessDeaths/',age,'.pdf'),
         width = 10, height = 4.75)
}

for(a in unique(df.pL$Age0)){
  age   <- a
  agegr <- substr(a, 2, nchar(a))
  
  if(age == "Y_GE90") agegr <- '90+'
  
  plot_excess_deaths(age = age, agegr = agegr)
}

