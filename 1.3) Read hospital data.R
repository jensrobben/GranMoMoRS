##### 0) Preliminary stuff ----

# Clear environment
rm(list = ls())
gc()

# Load packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble',
              'Polychrome','tidyr','ISOweek', 'ggplot2','RColorBrewer',
              'data.table', 'viridis', 'imputeTS', 'eurostat', 'RCurl',
              "eurostat")
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
file_ha <- 'Data/OPENDATA/donnees-hospitalieres-covid-19-dep-france.csv'

# Read csv file
df_ha <- read.csv(file = file_ha, sep = ';')

# One missing department name, but corresponds to Collectivity of Saint-Martin (overseas)
unique(df_ha$Nom.département[which(! df_ha$Nom.département %in% shapef3$NUTS_NAME)])

# Go from NUTS 3 to former French NUTS 2 regions
df <- df_ha %>%
  dplyr::left_join(shapef3, by = c('Nom.département' = 'NUTS_NAME')) %>%
  dplyr::filter(Sexe == "Tous") %>%
  arrange(Date, NUTS_ID) %>%
  dplyr::select(c('Date','NUTS_ID', 'Nb.actuellement.hospitalisés',
                  'Nb.actuellement.en.soins.intensifs',
                  'Nb.Quotidien.Admis.Hospitalisation')) %>%
  mutate(NUTS2_ID = substr(NUTS_ID, 1, 4),
         Date     = as.Date(Date),
         ISOYear  = isoyear(Date),
         ISOWeek  = isoweek(Date))

# Group information by week
df <- df %>% 
  group_by(ISOYear, ISOWeek, NUTS2_ID) %>%
  reframe(New.Hosp = sum(Nb.Quotidien.Admis.Hospitalisation, na.rm = TRUE)) %>%
  mutate(Date = ISOweek2date(paste0(ISOYear,'-W', 
                                    ifelse(ISOWeek < 10, paste0('0',ISOWeek),
                                           ISOWeek), '-',1)), .before = 1)


##### 3) Population data ----
# Download
P.xt.all <- fread(paste0('Data/EUROSTAT/demo_r_d2jan.tsv'), sep = '\t',
                  header = TRUE) %>% 
  separate(col = "freq,unit,sex,age,geo\\TIME_PERIOD", 
           into = c('Freq', 'Unit', 'Sex', 'Age', 'Region'), sep = ',') 

# Pre-process
P.xt.c <- P.xt.all %>% 
  dplyr::filter(Region %in% sort(unique(df$NUTS2_ID)),
                Age %notin% c('UNK','TOTAL')) %>%
  dplyr::select(-c('Freq','Unit')) %>%
  tidyr::gather(key = 'Time', value = 'Pop', - c(Sex, Age, Region)) %>%
  mutate(ISOYear = as.integer(Time),
         Age = factor(Age, levels = c('Y_LT1', paste0('Y',seq(1,99,1)), 'Y_OPEN')),
         .before = 1) %>%
  dplyr::select(-c('Time')) %>% 
  arrange(Sex, ISOYear, Age) %>%
  dplyr::filter(!grepl(':', Pop),
                as.integer(Age) >= 66, ISOYear <= 2024, ISOYear >= 2020, Sex == 'T')

# Pre-process
P.xt.c$Pop <- readr::parse_number(P.xt.c$Pop) %>% as.integer()
P.xt.c$Age <- droplevels(P.xt.c$Age)

Pt <- P.xt.c %>% ungroup() %>%
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

