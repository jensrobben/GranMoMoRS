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

##### 2) Download climate projections from CDS ----

# Select weather variable
var      <- 'tg'
var_name <- '2m_air_temperature'

# Login details - create account on CDS and fill in your API token 
wf_set_key(key = '....')

# Request information
request <- list(
  product_type = "essential_climate_variables",
  variable = var_name,
  processing_type = "bias_corrected",
  variable_type = "absolute_values",
  time_aggregation = "daily",
  horizontal_resolution = "5_km",
  experiment = "rcp_2_6",
  rcm = "racmo22e",
  gcm = "hadgem2_es",
  ensemble_member = "r1i1p1",
  period = as.character(seq(2024,2030)),
  format = "tgz",
  dataset_short_name = "sis-hydrology-meteorology-derived-projections",
  target = paste0(var,'.tar.gz'))

# Download (about 15 minutes)
wf_request(request = request, path = paste0(getwd(),'/Data/CDS'))

# Extract files from tar archive
lfile <- list.files('Data/CDS')
untar(paste0('Data/CDS/',lfile[2]), exdir = "Data/CDS")

# Delete original files
file.remove(paste0('Data/CDS/',lfile))

# NC files
nc_files <- list.files('Data/CDS')

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

### 4.1) Extract data ----

# Open the nc file
files  <- paste0('Data/CDS/', nc_files)
nc_ds  <- lapply(files, function(f) nc_open(f))

# France bounding box (lon_min, lon_max, lat_min, lat_max)
bbox <- c(long.range, lat.range)

# Read variables (it's important that all files have same lon, lat coordinates)
lon    <- ncvar_get(nc_ds[[1]], "lon")  # 2D array [x, y]
lat    <- ncvar_get(nc_ds[[1]], "lat")  # 2D array [x, y]
tas    <- lapply(nc_ds, function(nc) ncvar_get(nc, "tasAdjust"))  # 3D array [x, y, time]
tas    <- abind::abind(tas)
time   <- unlist(lapply(nc_ds, function(nc) ncvar_get(nc, "time")))

# Clear memory 
gc()

# Put as date object
t_units <- ncatt_get(nc_ds[[1]], "time", "units")
t_ustr  <- strsplit(t_units$value, " ")
t_dstr  <- strsplit(unlist(t_ustr)[3], "-")
date    <- ymd(t_dstr) + dseconds(time)

# Create a mask for grid points within the bounding box
mask <- (lon >= bbox[1]) & (lon <= bbox[2]) & 
  (lat >= bbox[3]) & (lat <= bbox[4])

# Find indices of x and y where mask is TRUE
indices   <- which(mask, arr.ind = TRUE)
x_indices <- indices[, 1]
y_indices <- indices[, 2]

# Get slice range for x and y
x_slice <- range(x_indices)
y_slice <- range(y_indices)

# Subset the data
tas_subset <- tas[x_slice[1]:x_slice[2], y_slice[1]:y_slice[2], ]
lon_subset <- lon[x_slice[1]:x_slice[2], y_slice[1]:y_slice[2]]
lat_subset <- lat[x_slice[1]:x_slice[2], y_slice[1]:y_slice[2]]

# Close the original NetCDF files
lapply(nc_ds, function(nc) nc_close(nc))

### 4.2) Intersection NUTS 2 regions with grid ----

# Redefine
long <- lon_subset
lat  <- lat_subset
arr  <- tas_subset

# All grid points in E-OBS data 
dat <- data.frame(as.vector(long), as.vector(lat))
colnames(dat) <- c('Longitude', 'Latitude')

# Put in same CRS
eobs_sf  <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
nuts3_sf <- st_transform(shapef, 4326)     

# Make intersection
ints <- st_intersects(nuts3_sf, eobs_sf, sparse = TRUE)
for(k in 1:length(ints)){
  dat[ints[[k]],'region'] <- nuts3_sf$NUTS_ID[k]
}

# Keep E-OBS grid points in NUTS 3 regions
coords.ctry <- dat %>% na.omit() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 


### 4.3) Calculate population weights ----

# Remove missing data
df <- dat %>% na.omit()

# Calculate population weights
wpop <- get_pop_weights(ctry.spec = ctry.spec, long.ext = long, lat.ext = lat)
df   <- df %>% dplyr::left_join(wpop, by = c('Longitude' = 'Long.ext', 
                                             'Latitude'  = 'Lat.ext', 
                                             'region'    = 'Region'))

# If NA - zero population weight
df[which(is.na(df), arr.ind = TRUE)] <- 0

# Small check for population weights to sum up to 1 for each region
test <- df %>% dplyr::group_by(region) %>%
  summarize(sum(WPC2000, na.rm = TRUE), sum(WPC2005, na.rm = TRUE), 
            sum(WPC2010, na.rm = TRUE), sum(WPC2015, na.rm = TRUE), 
            sum(WPC2020, na.rm = TRUE))


### 4.4) Average climate feature over grid points in each NUTS region ----

# Split coordinates by region
coord.region <- df %>% split(f = df$region)

# Function to calculate population weighted average over each region
avg.clim.region <- function(df.sub){
  message(df.sub$region[1])
  lola  <- data.frame('Longitude' = as.vector(long), 'Latitude' = as.vector(lat))
  indd <- sapply(1:nrow(df.sub), function(j) which(df.sub$Longitude[j] == lola$Longitude & 
                                                     df.sub$Latitude[j] == lola$Latitude))
  rowInd <- indd %% nrow(long)
  colInd <- floor(indd/nrow(long)) + 1
  obj <- sapply(1:nrow(df.sub), function(s) 
    arr[rowInd[s],colInd[s],])
  
  weight.matrix <- matrix(rep(df.sub$WPC2015, times = nrow(obj)), 
                          nrow = nrow(obj), ncol = ncol(obj), byrow = TRUE)
  mm <- rowSums((!is.na(obj)) * weight.matrix)
  
  objw <- obj * weight.matrix
  
  v.objw <- rowSums(objw, na.rm = TRUE)/mm
  v.objw - 273.15
}

df.region.CLI <- sapply(coord.region, function(sub) 
  avg.clim.region(sub))
rownames(df.region.CLI) <- as.character(substr(date, 1, 10))


### 4.5) Function to transform matrix to long data frame r (see ribiosUtils) ----
matrix2longdf <- function(mat,
                          row.names, col.names,
                          longdf.colnames=c("row","column","value")) {
  if(missing(row.names)) row.names <- rownames(mat)
  if(missing(col.names)) col.names <- colnames(mat)
  
  if(is.null(row.names)) row.names <- 1:nrow(mat)
  if(is.null(col.names)) col.names <- 1:ncol(mat)
  
  value <- as.vector(mat)
  if(length(row.names)!=nrow(mat))
    warning("row.names is inconsistent with the matrix dim")
  if(length(col.names)!=ncol(mat))
    warning("col.names is inconsistent with the matrix dim")
  
  rn <- rep(row.names, ncol(mat))
  cn <- rep(col.names, each=nrow(mat))
  res <- data.frame(row=rn,
                    column=cn,
                    value=value)
  colnames(res) <- longdf.colnames
  return(res)
}

df.region.CLI.daily  <- matrix2longdf(df.region.CLI,
                                      longdf.colnames = c('Date', 'Region', var))  

df.region.CLI.daily$Date <- as.Date(df.region.CLI.daily$Date)


# Save
file.save <- paste0("Data/CDS/FR_NUTS2_",var,'_daily_rcp26_2024_2030.rds')
saveRDS(df.region.CLI.daily, file = file.save)

