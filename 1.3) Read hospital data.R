##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble',
              'Polychrome','tidyr','ISOweek', 'ggplot2','RColorBrewer',
              'data.table', 'viridis', 'imputeTS', 'eurostat', 'RCurl')
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

### Read shape file with NUTS 2 codes (geo level of death counts)
shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp') %>%
  dplyr::filter(CNTR_CODE %in% ctry.spec) %>%
  dplyr::arrange(NUTS_NAME)
shapefile$NUTS_NAME[shapefile$NUTS_NAME == "Côte-d’Or"] <- "Côte-d'Or"
shapefile$NUTS_NAME[shapefile$NUTS_NAME == "Côtes-d’Armor"] <- "Côtes-d'Armor"
shapefile$NUTS_NAME[shapefile$NUTS_NAME == "Val-d’Oise"] <- "Val-d'Oise"

shapef2    <- shapefile[shapefile$LEVL_CODE == 2,]
shapef3    <- shapefile[shapefile$LEVL_CODE == 3,] 

##### 2) Read hospital data ----
# URL to OpenDataSoft for hospital data relating to COVID-19 epidemic
url_ha <- paste0('https://public.opendatasoft.com/api/explore/v2.1/catalog/',
                 'datasets/donnees-hospitalieres-covid-19-dep-france/exports/csv')

# Read csv file
df_ha <- read.csv(url_ha, sep = ";") 

# One missing department name, but corresponds to Collectivity of Saint-Martin (overseas)
unique(df_ha$nom_dep_min[which(! df_ha$nom_dep_min %in% shapef3$NUTS_NAME)])

# Go from NUTS 3 to former French NUTS 2 regions
df <- df_ha %>%
  dplyr::left_join(shapef3, by = c('nom_dep_min' = 'NUTS_NAME')) %>%
  dplyr::filter(sex == "Tous") %>%
  arrange(date, NUTS_ID) %>%
  dplyr::select(c('date','NUTS_ID', 'day_hosp_new')) %>%
  mutate(NUTS2_ID = substr(NUTS_ID, 1, 4),
         Date     = as.Date(date),
         ISOYear  = isoyear(date),
         ISOWeek  = isoweek(date))

# Group information by week
df <- df %>% 
  group_by(ISOYear, ISOWeek, NUTS2_ID) %>%
  reframe(New.Hosp = sum(day_hosp_new, na.rm = TRUE)) %>%
  mutate(Date = ISOweek2date(paste0(ISOYear,'-W', 
                                    ifelse(ISOWeek < 10, paste0('0',ISOWeek),
                                           ISOWeek), '-',1)), .before = 1)


##### 3) Population data ----
# Download
P.xt.all <- get_eurostat(id = 'demo_r_d2jan', cache = FALSE,
                         compress_file = FALSE, time_format = 'raw')

# Pre-process
P.xt.c <- P.xt.all %>% 
  dplyr::filter(geo %in% sort(unique(df$NUTS2_ID)),
                age %notin% c('UNK','TOTAL')) %>%
  dplyr::select(-c('freq','unit')) %>%
  mutate(ISOYear = as.integer(TIME_PERIOD),
         Age = factor(age, levels = c('Y_LT1', paste0('Y',seq(1,99,1)), 'Y_OPEN')),
         .before = 1) %>%
  dplyr::select(-c('TIME_PERIOD', 'age')) %>% 
  arrange(sex, ISOYear, Age) %>%
  dplyr::filter(!grepl(':', values),
                as.integer(Age) >= 66, ISOYear <= 2024, ISOYear >= 2020, sex == 'T')

colnames(P.xt.c) <- c('ISOYear', 'Age', 'Sex', 'Region', 'Pop')

P.xt.c$Age <- droplevels(P.xt.c$Age)

Pt <- P.xt.c %>% 
  ungroup() %>%
  group_by(ISOYear, Region) %>%
  reframe('Pop' = sum(as.numeric(Pop))) 

##### 4) Merge hospital and population data ----

# Number of new hospital admissions per 1000 inhabitants
df <- df %>% 
  dplyr::left_join(Pt, by = c('ISOYear' = 'ISOYear', 'NUTS2_ID' = 'Region')) %>%
  mutate(Nhosp = New.Hosp/Pop * 1000) %>%
  dplyr::filter(! NUTS2_ID %in% c('FRY1','FRY2','FRY3','FRY4','FRY5','FRM0')) %>%
  na.omit()

df.final <- df %>% 
  dplyr::select(c('Date', 'ISOYear', 'ISOWeek', 'NUTS2_ID', 'Nhosp'))

# Save
saveRDS(df.final, file = 'Data/OPENDATA/Hosp_NUTS2.rds')

