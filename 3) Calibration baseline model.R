##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'stats',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr','httr',
              'utils', 'sp', 'spdep')
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

##### 2) Load weekly deaths and preprocess ----

# Weekly death counts from Eurostat
d.xtw.all <- get_eurostat(id = 'demo_r_mweek3', cache = FALSE,
                          compress_file = FALSE, time_format = 'raw')

colnames(d.xtw.all) <- c('Freq', 'Unit', 'Sex', 'Age', 
                         'Region', 'Time', 'Deaths')

# Filter and pre-process
d.xtw.c <- d.xtw.all %>% 
  dplyr::filter(Region %in% shapef$NUTS_ID &
                  ! Age %in% c('UNK', 'TOTAL')) %>%
  mutate(ISOYear = as.integer(substr(Time, 1, 4)), 
         ISOWeek = as.integer(substr(Time, 7, 9)),
         Age = factor(Age, levels = c('Y_LT5', paste0('Y',seq(5,85,5), '-', 
                                                      seq(9,89,5)), 'Y_GE90')),
         .before = 1) %>%
  dplyr::select(-c('Freq', 'Unit', 'Time')) %>% 
  dplyr::arrange(Sex, ISOYear, ISOWeek, Age) %>%
  dplyr::filter(!grepl(':', Deaths),
                as.integer(Age) >= 14, 
                ISOYear >= year(t_min) & 
                  ((ISOYear < year(t_max)) | (ISOYear == year(t_max) & ISOWeek <= 26)))

# Drop levels
d.xtw.c$Age    <- droplevels(d.xtw.c$Age)

# Focus on age groups above 65 / combined data on both sexes
d.xtw.c <- d.xtw.c %>% 
  dplyr::filter(Age %in% c('Y65-69', 'Y70-74', 'Y75-79', 
                           'Y80-84', 'Y85-89', 'Y_GE90') & 
                  Sex == 'T')

# Add date
wdate   <- ISOweek2date(sprintf("%d-W%02d-%d", d.xtw.c$ISOYear, d.xtw.c$ISOWeek, 1))
d.xtw.c <- d.xtw.c %>% 
  mutate('Date' = wdate, .before = 1)

# Remove Sex column
df.d <- d.xtw.c %>% 
  dplyr::select(-c(Sex))

##### 3) Weekly population exposures at January 1, year t ----

# Annual exposures from Eurostat
P.xt.all <- get_eurostat(id = 'demo_r_d2jan', cache = FALSE,
                         compress_file = FALSE, time_format = 'raw')

colnames(P.xt.all) <- c('Freq', 'Unit', 'Sex', 
                        'Age', 'Region', 'Time', 'Pop')

# Annual exposures- countries of interest
P.xt.c <- P.xt.all %>% 
  dplyr::filter(Region %in% shapef$NUTS_ID,
                ! Age %in% c('UNK','TOTAL')) %>%
  mutate(ISOYear = as.integer(Time),
         Age = factor(Age, levels = c('Y_LT1', paste0('Y',seq(1,99,1)), 'Y_OPEN')),
         .before = 1) %>%
  dplyr::select(-c('Freq', 'Unit', 'Time')) %>% 
  arrange(Sex, ISOYear, Age) %>%
  dplyr::filter(!grepl(':', Pop),
                as.integer(Age) >= 66, ISOYear <= year(t_max) + 1, ISOYear >= year(t_min))

# Drop levels
P.xt.c$Age <- droplevels(P.xt.c$Age)

# Focus on unisex data 
P.xt.c <- P.xt.c %>% 
  dplyr::filter(Sex == "T") %>%
  dplyr::select(-c(Sex))

# Extrapolation for up to the year 2025 (January 1)
list.FR <- P.xt.c %>% 
  split(~ Region + Age) 
newdf   <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(newdf) <- c('ISOYear', 'Age', 'Region', 'Pop')

for(s in names(list.FR)){
  # Population counts for region r
  sub   <- list.FR[[s]]
  
  # LM fit
  lmr <- mgcv::gam(Pop ~ s(ISOYear), data = sub)
  
  # Missing year(s)
  years   <- year(t_min):(year(t_max)+1)
  missing <- years[!years %in% sub$ISOYear]
  
  # Prediction on missing year(s)
  pred    <- predict(lmr, newdata = data.frame(ISOYear = missing))
  
  # Add to new data frame
  newdf <- rbind(newdf, data.frame('ISOYear' = missing, 
                                   'Region'  = substr(s,1,4),
                                   'Age'     = substr(s,6,15), 
                                   'Pop'     = pred))
}

P.xt.c <- rbind(P.xt.c, newdf)

# Create weekly exposures from population count
E.xt.c <- P.xt.c %>% 
  dplyr::group_by(Age, Region) %>% 
  arrange(ISOYear) %>% 
  reframe(ISOYear, 
          'Expo' = c((Pop[-length(ISOYear)] + Pop[-1])/2, 
                     Pop[length(ISOYear)])/52.1775) %>%
  dplyr::filter(ISOYear <= year(t_max))

# Group per age group
E.xt.c$Age <- plyr::mapvalues(E.xt.c$Age, from = levels(E.xt.c$Age),
                              to = c(rep('Y65-69', 5), rep('Y70-74', 5),
                                     rep('Y75-79', 5), rep('Y80-84', 5),
                                     rep('Y85-89', 5), rep('Y_GE90', 11)))
E.xt.c <- E.xt.c %>%
  group_by(ISOYear, Age, Region) %>%
  reframe(Expo = sum(Expo)) %>%
  ungroup()

# Add to death counts data 
df.d <- df.d %>% 
  dplyr::left_join(E.xt.c, by = c('ISOYear', 'Age', 'Region')) 

##### 4) Fit penalized Poisson GLM ----

# Arrange data set
df.d <- df.d %>% 
  arrange(Date, Region, Age)
shapef <- shapef %>% 
  arrange(NUTS_ID)

# Add Fourier terms
df <- df.d %>%
  mutate('fsin52' = sin(2*pi*ISOWeek/52.1775),
         'fcos52' = cos(2*pi*ISOWeek/52.1775),
         'fsin26' = sin(4*pi*ISOWeek/52.1775),
         'fcos26' = cos(4*pi*ISOWeek/52.1775))

# Factorize region
df$Region <- factor(df$Region, levels = regions)

# Create year-indicator (t = 0,1,2,..)
df$Time <- df$ISOYear  - min(df$ISOYear)

# Formulas
formula0 <- Deaths ~ -1 + Region:Age
formula1 <- Deaths ~ -1 + Time:Region:Age
formula2 <- Deaths ~ -1 + fsin52:Region:Age
formula3 <- Deaths ~ -1 + fcos52:Region:Age
formula4 <- Deaths ~ -1 + fsin26:Region:Age
formula5 <- Deaths ~ -1 + fcos26:Region:Age

# Model matrices
X0 <- model.matrix(formula0, data = df)
X1 <- model.matrix(formula1, data = df)
X2 <- model.matrix(formula2, data = df)
X3 <- model.matrix(formula3, data = df)
X4 <- model.matrix(formula4, data = df)
X5 <- model.matrix(formula5, data = df)

# Add model matrices to dataframe
df$M0 <- X0
df$M1 <- X1
df$M2 <- X2
df$M3 <- X3
df$M4 <- X4
df$M5 <- X5

# Save data frame
saveRDS(df, file = 'Results/df_GAM.rds')

# Adjacency matrix
adjmat <- nb2mat(poly2nb(shapef[,'geometry']), zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat) <- list(shapef[which(shapef$NUTS_ID %in% regions),]$NUTS_ID, 
                         shapef[which(shapef$NUTS_ID %in% regions),]$NUTS_ID)
diag(adjmat) <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))

# Expand for different age categories
nR  <- nrow(adjmat)
nA  <- length(unique(df.d$Age))
ADJ <- matrix(0, nrow = nR*nA, ncol = nR*nA)
for(a in 1:nA){
  ADJ[((a-1)*nR + 1):(a*nR),((a-1)*nR + 1):(a*nR)] <- adjmat
}

# Formula for gam fit
formula <- Deaths ~ -1 + M0 + M1 + M2 + M3 + M4 + M5

# Remove extreme COVID-19 weeks
dfs <- df %>% dplyr::filter(! (ISOYear == 2020 & ISOWeek %in% 13:16))

# Baseline fit
fit2 <- mgcv::gam(formula, offset = log(dfs$Expo), data = dfs,
                  family = poisson(link = 'log'),
                  control = list(epsilon = 10^(-5), trace = TRUE,
                                 newton = list(conv.tol = 1e-5), nthreads = 10),
                    drop.intercept = TRUE,
                  paraPen = list(M0 = list(ADJ, rank = nA*nR, sp = -1),
                                 M1 = list(ADJ, rank = nA*nR, sp = -1),
                                 M2 = list(ADJ, rank = nA*nR, sp = -1),
                                 M3 = list(ADJ, rank = nA*nR, sp = -1),
                                 M4 = list(ADJ, rank = nA*nR, sp = -1),
                                 M5 = list(ADJ, rank = nA*nR, sp = -1)))

# Save / Read fitted GAM object
saveRDS(fit2, 'Results/GAM.rds')
fit2 <- readRDS('Results/GAM.rds')


# Save baseline deaths in data-frame
df.fitb <- data.frame(df[,c('Date', 'ISOYear', 'ISOWeek', 'Region', 'Age', 'Deaths','Expo')],
                      'bDeaths' = predict(fit2, df, type = 'response') * df$Expo)
saveRDS(df.fitb, 'Results/df_bDeaths.rds')
