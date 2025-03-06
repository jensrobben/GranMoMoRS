##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble',
              'Polychrome','tidyr','ISOweek', 'ggplot2','RColorBrewer',
              'data.table', 'viridis', 'imputeTS')
invisible(sapply(packages, require, character.only = TRUE))
`%notin%` <- Negate(`%in%`)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

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


##### 2) Download Sentinel data ----

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

# Check if data for every isoyear, isoweek
check1 <- df_ia %>% 
  group_by(ISOYear, ISOWeek) %>%
  summarize(l = length(ia100))
all(check1$l == nrow(shapef)) # 21 is number of regions we consider

# Save
file.save <- paste0("Data/SENTIWEB/", ctry.spec, '_NUTS2_ia_weekly','.rds')
saveRDS(df_ia, file = file.save)

