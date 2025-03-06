##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'stats',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr','httr',
              'utils')
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

# Countries of interest
ctry.spec <- c("FR")

# Colour scales
ablue  <- '#56B4E9'
agreen <- '#A3D596'
ared   <- '#F8766D'

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

##### 2) Load daily temperature data + weekly transformation ----

# Read processed temperature data
df.t <- df_itime %>% 
  left_join(readRDS(file = 'Data/E-OBS/FR_NUTS2_tg_daily.rds'), by = 'Date') %>%
  mutate(Day = yday(Date), .before = 'ISOWeek') 

# Create temperature anomalies and heat/cold indicators
excess_temp <- function(r){
  # Data from region r
  dfsub.daily <- df.t %>% dplyr::filter(Region == r)
  
  # Baseline fits
  fitr.tg <- robust::lmRob(tg ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25),
                 data = dfsub.daily)
  
  # Predictions
  pred.tg <- predict(fitr.tg, dfsub.daily)
  
  # Quantiles
  qtg <- quantile(dfsub.daily$tg, c(0.05, 0.95), na.rm = TRUE)

  # Information
  dd <- data.frame(Date    = dfsub.daily$Date, Region = dfsub.daily$Region,
                   tg_anom = dfsub.daily$tg - pred.tg,
                   Tind95  = (dfsub.daily$tg > qtg[2])*1,
                   Tind5   = (dfsub.daily$tg < qtg[1])*1)
  dd
}

list.exc.t <- lapply(sort(unique(df.t$Region)), excess_temp)
df.exc.t   <- do.call('rbind', list.exc.t)
df.t       <- df.t %>% left_join(df.exc.t)

# Transform to weekly - on ISO week basis
df.t.w <- df.t %>% group_by(ISOYear, ISOWeek, Region) %>%
  summarize(w_avg_tg      = mean(tg), 
            w_avg_tg_anom = mean(tg_anom), 
            w_avg_Tind95  = mean(Tind95), 
            w_avg_Tind5   = mean(Tind5)) %>%
  ungroup() 

##### 3) Influenza activity ----

# Read influenza data
df_ia <- readRDS(file = 'Data/SENTIWEB/FR_NUTS2_ia_weekly.rds')

# Per 100 inhabitants
df_ia$ia100  <- df_ia$ia100/1000

# Influenza anomalies
ia.seasonal <- df_ia %>%
  dplyr::group_by(Region, ISOWeek) %>%
  summarize('avg_ia100' = mean(ia100))
df_ia <- df_ia %>%
  dplyr::left_join(ia.seasonal, by = c('ISOWeek', 'Region')) %>%
  mutate(ia100 = ia100 - avg_ia100)

# Quantile excess
df_ia$ia100q  <- pmax(df_ia$ia100 - quantile(df_ia$ia100, 0.75), 0)

##### 4) Hospitalization rates for COVID-19 ----

# Read data on hospitalizations for COVID-19
df_hosp <- readRDS(file = 'Data/OPENDATA/Hosp_NUTS2.rds')

# Quantile excess
df_hosp$Nhospq  <- pmax(df_hosp$Nhosp - 
                          quantile(df_hosp$Nhosp[
                            df_hosp$Date > as.Date('2020-03-16')], 0.75))

##### 5) Combine all data sets ----
df <- df.t.w %>%
  left_join(df_ia,   by = c('ISOYear', 'ISOWeek', 'Region' = 'NUTS_ID')) %>%
  left_join(df_hosp, by = c('ISOYear', 'ISOWeek', 'Date', 'Region' = 'NUTS2_ID')) 

df$Nhosp[is.na(df$Nhosp)]   <- 0
df$Nhospq[is.na(df$Nhospq)] <- 0

file <- paste0('Data/df_NUTS2_FR.rds')
saveRDS(df, file)

