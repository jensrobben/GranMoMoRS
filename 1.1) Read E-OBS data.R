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

# Source function
source('Functions/F. Diverse.R')
source('Functions/F. Population Weights.R')

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

##### 2) Download E-OBS data (CDS) ----

# Select weather variable
var      <- 'tg'
var_name <- 'mean_temperature'

# Login details - create account on CDS and fill in your API token 
wf_set_key(key = '....')

# Request information
request <- list(
  product_type = "ensemble_mean",
  variable = var_name,
  grid_resolution = "0_1deg",
  period = "2011_2024",
  version = "30_0e",
  format = "tgz",
  dataset_short_name = "insitu-gridded-observations-europe",
  target = paste0(var,'.tar.gz'))

# Download (about 15 minutes)
wf_request(request = request, path = paste0(getwd(),'/Data/E-OBS'))

# Extract files from tar archive
lfile <- list.files('Data/E-OBS')
untar(paste0('Data/E-OBS/',lfile), exdir = "Data/E-OBS")

# Delete original files
file.remove(paste0('Data/E-OBS/',lfile))

# NC files
lfile <- paste0(var, "_ens_mean_0.1deg_reg_2011-2024_v30.0e.nc")

##### 3) Download population count data (SEDAC) ----

# Create account on EarthData (https://urs.earthdata.nasa.gov/) 

# Replace with your EarthData credentials
username <- "...."
password <- "...."

# URL of the data
url <- paste0("https://sedac.ciesin.columbia.edu/downloads",
              "/data/gpw-v4/gpw-v4-population-count-rev11/",
              "gpw-v4-population-count-rev11_totpop_2pt5_min_nc.zip")

# Send GET request with authentication (approx 10 min.)
response <- GET(
  url,
  authenticate(username, password, type = "basic"),
  write_disk("Data/SEDAC/PopulationCount.zip", overwrite = TRUE)
)

# Check if download succeeded
if (status_code(response) == 200) {
  message("Data downloaded successfully!")
} else {
  message("Download failed. Check credentials or URL.")
}

# Unzip nc file
unzip(zipfile = "Data/SEDAC/PopulationCount.zip", 
      files = "gpw_v4_population_count_rev11_2pt5_min.nc",
      exdir = "Data/SEDAC", overwrite = FALSE)

##### 4) Pre-process E-OBS data and save ----
### 4.0) Weather variable ----
filename <- lfile
var      <- substr(filename,1,2)
message(var)

### 4.1) Extract data ----

# Open the nc file
nc_ds <- nc_open(paste0('Data/E-OBS/',filename))

# Extract longitude, latitude coordinates and time object
long <- ncvar_get(nc_ds, "longitude") 
lat  <- ncvar_get(nc_ds, "latitude")
time <- ncvar_get(nc_ds, "time")

# Put as date object
t_units <- ncatt_get(nc_ds, "time", "units")
t_ustr  <- strsplit(t_units$value, " ")
t_dstr  <- strsplit(unlist(t_ustr)[3], "-")
date    <- ymd(t_dstr) + ddays(time)

# Filter spatial data frames 
id.long <- which(nc_ds$dim$longitude$vals >= long.range[1] & 
                   nc_ds$dim$longitude$vals <= long.range[2])
id.lat  <- which(nc_ds$dim$latitude$vals >= lat.range[1] &
                   nc_ds$dim$latitude$vals <= lat.range[2])
id.time <- which(nc_ds$dim$time$vals >=  as.Date('2013-01-01') - 
                   as.Date('1950-01-01') & nc_ds$dim$time$vals <= 
                   as.Date('2024-12-31') - as.Date('1950-01-01'))
long    <- long[id.long]
lat     <- lat[id.lat]
date    <- date[id.time]

# Store the data in a 3-dimensional array 
var.array <- ncvar_get(nc_ds) 
var.array <- var.array[id.long, id.lat, id.time]
dimnames(var.array) <- list(long, lat, as.character(date))

### 4.2) Intersection NUTS 2 regions with E-OBS grid ----

# All grid points in E-OBS data 
dat <- expand.grid(long, lat)
colnames(dat) <- c('Longitude', 'Latitude')

# Put in same CRS
eobs_sf  <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
nuts2_sf <- st_transform(shapef, 4326)     

# Make intersection
ints <- st_intersects(nuts2_sf, eobs_sf, sparse = TRUE)
for(k in 1:length(ints)){
  dat[ints[[k]],'region'] <- nuts2_sf$NUTS_ID[k]
}

# Keep E-OBS grid points in NUTS 2 regions
coords.ctry <- dat %>% na.omit() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 

### 4.3) Calculate population weights ----

# Remove missing data
df <- dat %>% na.omit()

# Calculate population weights
wpop <- get_pop_weights(ctry.spec = ctry.spec, long.ext = long, lat.ext = lat)
df   <- df %>% left_join(wpop, by = c('Longitude' = 'Long.ext', 
                                      'Latitude'  = 'Lat.ext', 
                                      'region'    = 'Region'))

# Check for population weights to sum up to 1 for each region
test <- df %>% dplyr::group_by(region) %>%
  summarize(sum(WPC2000, na.rm = TRUE), sum(WPC2005, na.rm = TRUE), 
            sum(WPC2010, na.rm = TRUE), sum(WPC2015, na.rm = TRUE), 
            sum(WPC2020, na.rm = TRUE))

### 4.4) Average climate feature over grid points in each NUTS region ----

# Split coordinates by region
coord.region <- df %>% split(f = df$region)

# Function to calculate population weighted average over each region
avg.clim.region <- function(df.sub){
  lo.ind <- match(df.sub$Longitude, long)
  la.ind <- match(df.sub$Latitude, lat)
  obj <- sapply(1:nrow(df.sub), function(s) 
    var.array[lo.ind[s],la.ind[s],])
  
  weight.matrix <- matrix(rep(df.sub$WPC2015, times = nrow(obj)), 
                          nrow = nrow(obj), ncol = ncol(obj), byrow = TRUE)
  mm <- rowSums((!is.na(obj)) * weight.matrix)
  
  objw <- obj * weight.matrix
  
  v.objw <- rowSums(objw, na.rm = TRUE)/mm
  names(objw) <- dimnames(var.array)[[3]]
  v.objw
}

df.region.CLI <- sapply(coord.region, function(sub) 
  avg.clim.region(sub))


### 4.5) Function to transform matrix to long data frame r (see ribiosUtils) ----
df.region.CLI.daily  <- matrix2longdf(df.region.CLI,
                                      longdf.colnames = c('Date', 'Region', var))  

df.region.CLI.daily$Date <- as.Date(df.region.CLI.daily$Date)

# Save
file.save <- paste0("Data/E-OBS/",
                    paste0(ctry.spec, collapse = '_'),'_NUTS2_',var,'_daily','.rds')
saveRDS(df.region.CLI.daily, file = file.save)




