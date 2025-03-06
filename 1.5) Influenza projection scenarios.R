##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'stats',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr','httr',
              'utils', 'R2jags')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source function
source('Functions/F. Diverse.R')
source('Functions/F. Population Weights v2.R')

# Country/countries of interest
ctry.spec <- c("FR")

##### 1) Load shape file NUTS 2 European regions ----

# Read shape file - downloaded from Eurostat
shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp')

# Extract shapefile for NUTS 2 regions in the countries of interest 
nuts.spec <- 2
shapef    <- shapefile[shapefile$CNTR_CODE %in% ctry.spec & 
                         shapefile$LEVL_CODE == nuts.spec,]

# Remove overseas areas
overseas  <- c(paste0('FRY',1:5), 'FRM0')
ind.rm    <- unlist(sapply(overseas, function(x) 
  which(grepl(x, shapef$NUTS_ID, fixed = TRUE))))
shapef    <- shapef[-c(ind.rm), ] %>%
  dplyr::arrange(NUTS_ID)

# Box long lat coordinates where country falls into
coords_ctry <- st_coordinates(st_transform(shapef, 4326)) %>% 
  as.data.frame() %>% dplyr::select('X','Y')
long.range <- range(coords_ctry[,1])
lat.range  <- range(coords_ctry[,2])

##### 2) Download Sentinelles data ----

# URL to French Sentinel influenza data
url_ia <- 'https://www.sentiweb.fr/datasets/all/inc-3-REG.csv'

# Read (no data for overseas region 'OUTRE-MER)
df_ia <- read.csv(url_ia, header = TRUE, skip = 1) %>%
  mutate('ISOYear' = as.integer(substr(week, 1, 4)),
         'ISOWeek' = as.integer(substr(week, 5, 6)),
         'Date'    = ISOweek2date(sprintf("%d-W%02d-%d", ISOYear, ISOWeek, 1)),
         'ia100'   = inc100,
         .before   = 1) %>%
  dplyr::filter(! geo_name %in% c('OUTRE-MER','CORSE')) %>% 
  dplyr::select(c(ISOYear, ISOWeek, Date, ia100, geo_name)) %>%
  dplyr::arrange(Date, geo_name)

# Match with NUTS 2 regions
df_ia$Region <- plyr::mapvalues(df_ia$geo_name, 
                                from = sort(unique(df_ia$geo_name)),
                                to   = sort(unique(shapef$NUTS_NAME)))
df_ia <- df_ia %>%
  dplyr::left_join(shapef[,c('NUTS_ID','NUTS_NAME')], 
                   by = c('Region' = 'NUTS_NAME')) %>%
  arrange(ISOYear, ISOWeek) %>%
  dplyr::filter(Date >= as.Date('2012-12-31') & 
                  Date <= as.Date('2024-06-30')) %>%
  mutate(ia100 = as.integer(ia100))


##### 3) Create forecast scenarios using Bayesian SIRS model ----
# Forecasting period
time.proj   <- seq(max(df_ia$Date) + 7, max(df_ia$Date) + 7*52*5 + 7*27, by = 'week')
nT          <- length(unique(df_ia$Date))
nT_forecast <- length(unique(time.proj))  # Number of time steps to predict
nT_total    <- nT + nT_forecast

# JAGS Model Implementation for the Bayesian DT-SIRS Model 1 (M1)
jags_model <- "model {
  ### Define likelihood 
  
  # At t = 2
  S[2] <- S[1] - S[1] * 1/N * (2 * I[1])^(alpha) * exp(phi[mvec[2]] + P[2]) + psi * R[1]
  R[2] <- (1 - psi) * R[1] + I[1]
  P[2] ~ dunif(0,d)
  
  # For t >= 3
  for (t in 3:nT) {
    # Infections
    I[t] ~ dpois(lambda[t] + 1) 
    P[t] ~ dunif(0,d)
    lambda[t] <- S[t-1] * 1/N * (I[t-1] + I[t-2])^(alpha) * exp(phi[mvec[t]] + P[t])
    
    # Observed compartments
    S[t] <- S[t-1] - S[t-1] * 1/N * (I[t-1] + I[t-2])^(alpha) * exp(phi[mvec[t]] + P[t]) + psi * R[t-1]
    R[t] <- (1 - psi) * R[t-1] + I[t-1]
  }
  
  ### In-sample posterior predictions
  
  # At t = 2
  S_pred[2] <- S_pred[1] - S_pred[1] * 1/N * (2 * I[1])^(alpha) * exp(phi[mvec[2]] + P_pred[2]) + psi * R_pred[1]
  R_pred[2] <- (1 - psi) * R_pred[1] + I[1]
  P_pred[2] ~ dunif(0,d)
  
  # For t >= 3
  for(t in 3:(nT)){
    # In-sample predictions
    I_pred[t] ~ dpois(lambda_pred[t] + 1)  
    P_pred[t] ~ dunif(0,d)
    lambda_pred[t] <- S_pred[t-1] * 1/N * (I[t-1] + I[t-2])^(alpha) * exp(phi[mvec[t]] + P_pred[2])
    
    # Compartments
    S_pred[t] <- S_pred[t-1] - S_pred[t-1] * 1/N * (I[t-1] + I[t-2])^(alpha) * exp(phi[mvec[t]] + P_pred[2]) + psi * R_pred[t-1]
    R_pred[t] <- (1 - psi) * R_pred[t-1] + I[t-1]
  }
  
  # For t >= nT + 1
  for(t in (nT+1):(nT_total)){
    # Forecasts
    I_pred[t] ~ dpois(lambda_pred[t] + 1)  
    P_pred[t] ~ dunif(0,d)
    lambda_pred[t] <- S_pred[t-1] * 1/N * (I_pred[t-1] + I_pred[t-2])^(alpha) * exp(phi[mvec[t]] + P_pred[2])
    
    # Compartments
    S_pred[t] <- S_pred[t-1] - S_pred[t-1] * 1/N * (I_pred[t-1] + I_pred[t-2])^(alpha) * exp(phi[mvec[t]] + P_pred[2]) + psi * R_pred[t-1]
    R_pred[t] <- (1 - psi) * R_pred[t-1] + I_pred[t-1]
  }
  
  # Priors model parameters
  psi ~ dunif(0.5,1)  # Proportion of infected recovering
  alpha ~ dunif(0,2)
  d = 0.6
  for(t in 1:13){
    phi[t] ~ dnorm(0,1)
  }
  
  # Initial conditions
  S[1] <- S0 - I[1]
  R[1] <- R0
  S_pred[1] <- S0 - I[1]
  R_pred[1] <- R0
  I_pred[1] <- I[1]
}"

# Save the JAGS model to a file
writeLines(jags_model, "Results/dt_sirs_model.jags")

# Parameters to monitor
params <- c("phi", "psi", "alpha", "d", "I_pred", "S_pred", "R_pred")

# Fit for every region separately
for(r in regions) {
  # Filter dataset to region r
  dfr    <- df_ia %>% dplyr::filter(NUTS_ID == r)
  
  # Projection data set
  dfr.proj    <- data.frame("Year" = c(dfr$ISOYear, isoyear(time.proj)),
                            "Week" = c(dfr$ISOWeek, isoweek(time.proj)),
                            "Date" = c(dfr$Date, time.proj),
                            "ia100" = rep(NA, nT_total),
                            "geo_name" = rep(NA, dfr$geo_name[1], nT_total),
                            "Region" = rep(NA, dfr$Region[1], nT_total),
                            "NUTS_ID" = rep(NA, dfr$NUTS_ID[1], nT_total),
                            "geometry" = rep(NA, dfr$geometry[1], nT_total))
  
  
  # Construct mvec
  mvec <- rep(NA, length(nT_total))
  for(t in 1:nT_total){
    w <- dfr.proj$Week[t]
    if(w %% 4 == 0){
      mvec[t] <- floor(w/4)
    } else {
      mvec[t] <- min(floor(w/4) + 1, 13)
    }
  }
  
  # Define the data list as argument for jags model
  data_list <- list(
    nT = nT,
    nT_total = nT_total,
    I = c(dfr$ia100, rep(NA, nT_forecast)),
    S0 = 100000,
    R0 = 0,
    N = 100000,
    mvec = mvec
  )
  
  #sink(paste0("Report.txt"), append = TRUE)
  cat("Region ",r)
  
  # Run the JAGS model
  model_fit <- jags.parallel(
    data = data_list,
    parameters.to.save = params,
    model.file = "Results/dt_sirs_model.jags",
    n.chains = 3,
    n.iter = 5000,
    n.burnin = 2500,
    n.thin = 5
  )
  
  saveRDS(model_fit, paste0('Results/SIRS models/SIR_',r,'.rds'))
}

