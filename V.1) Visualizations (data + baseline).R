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


##### 2) Load data sets -----

# Features
df.feat    <- readRDS('Data/df_NUTS2_FR.rds')
df.bDeaths <- readRDS('Results/df_bDeaths.rds')

df <- df.bDeaths %>% 
  left_join(df.feat, 
            by = c('Date', 'ISOYear', 'ISOWeek', 'Region')) 

##### 3) Visualisations ----

### 3.1) Visualisation model (tikz) ----
regions <- shapef$NUTS_ID
dfsubar <- df %>% dplyr::filter(Region == regions[12], Age == "Y_GE90")
fits <- glm(Deaths ~ ISOYear + sin(2*pi*ISOWeek/52) + cos(2*pi*ISOWeek/52), data = dfsubar,
            family = poisson(link = 'log'))
p3.1 <- ggplot(dfsubar) +
  theme_minimal() + 
  geom_line(aes(x = Date, y = Deaths), col = 'gray30') + 
  geom_line(aes(x = Date, y = fits$fitted.values), col = 'red', linewidth = 0.8) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p3.1

### 3.2) Aggregated death counts/exposures per week ----

# Plot death counts by age group
dff1 <- df %>%
  group_by(Region, Date) %>%
  summarize(Deaths = sum(Deaths),
            Expo = sum(Expo))

regg <- dff1 %>%
  group_by(Region) %>%
  summarize('TD' = sum(Deaths)) %>%
  arrange(TD) %>%
  pull(Region)

dff1$Region <- factor(dff1$Region, levels = regg)
soft_colors <- rev(c(
  "#A6CEE3", "#7FC97F", "#BEAED4", "#FAA086", "#FFFF99", 
  "#386CB0", "#F28E2B", "#FF9DA7", "#9FC17D", "#86BCB6", 
  "#E15759", "#FFC988", "#B07AA1", "#D37295", "#9D7660", 
  "#BAB0AC", "#79706E", "#D4A6C8", "#6D7660", "#C0CBE8", 
  "#F1CE63"
)) 
names(soft_colors) <- unique(shapef$NUTS_ID)


p1 <- ggplot(dff1) + 
  theme_bw(base_size = 15) + 
  geom_bar(aes(x = Date, y = Deaths, group = Region, fill = Region),
           stat = 'identity', width = 7) + 
  scale_fill_manual(name = 'NUTS 2',
                    values = soft_colors,
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  scale_x_continuous(breaks = c(as.Date('2013-01-01'),
                                as.Date('2015-01-01'),
                                as.Date('2017-01-01'),
                                as.Date('2019-01-01'),
                                as.Date('2021-01-01'),
                                as.Date('2023-01-01')),
                     labels = seq(2013,2023,2))

p2 <- ggplot(dff1) + 
  theme_bw(base_size = 15) + ylab('Population exposure') + 
  geom_bar(aes(x = Date, y = Expo, group = Region, fill = Region),
           stat = 'identity', width = 7) + 
  scale_fill_manual(name = 'NUTS 2',
                    values = soft_colors,
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  scale_x_continuous(breaks = c(as.Date('2013-01-01'),
                                as.Date('2015-01-01'),
                                as.Date('2017-01-01'),
                                as.Date('2019-01-01'),
                                as.Date('2021-01-01'),
                                as.Date('2023-01-01')),
                     labels = seq(2013,2023,2))



p3.2 <- ggpubr::ggarrange(plotlist = list(p1, p2),
                          common.legend = TRUE, 
                          legend = 'bottom',
                          align = 'hv')
p3.2

### 3.3) Aggregated death counts/exposures per age group ----
dff2 <- df %>%
  group_by(Age, Region) %>%
  summarize(Deaths = sum(Deaths),
            Expo = sum(Expo))

regg <- dff2 %>%
  group_by(Region) %>%
  summarize('TD' = sum(Deaths)) %>%
  arrange(TD) %>%
  pull(Region)

dff2$Region <- factor(dff2$Region, levels = regg)
soft_colors <- rev(c(
  "#A6CEE3", "#7FC97F", "#BEAED4", "#FAA086", "#FFFF99", 
  "#386CB0", "#F28E2B", "#FF9DA7", "#9FC17D", "#86BCB6", 
  "#E15759", "#FFC988", "#B07AA1", "#D37295", "#9D7660", 
  "#BAB0AC", "#79706E", "#D4A6C8", "#6D7660", "#C0CBE8", 
  "#F1CE63"
))
names(soft_colors) <- unique(shapef$NUTS_ID)


p1 <- ggplot(dff2) + 
  theme_bw(base_size = 15) + 
  geom_bar(aes(x = Age, y = Deaths, group = Region, fill = Region),
           stat = 'identity') + 
  scale_fill_manual(name = 'NUTS 2',
                    values = soft_colors,
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 

p2 <- ggplot(dff2) + 
  theme_bw(base_size = 15) + ylab('Population exposure') + 
  geom_bar(aes(x = Age, y = Expo, group = Region, fill = Region),
           stat = 'identity') + 
  scale_fill_manual(name = 'NUTS 2',
                    values = soft_colors,
                    breaks = rev(regg)) + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 


p3.3 <- ggpubr::ggarrange(plotlist = list(p1, p2),
                          common.legend = TRUE, 
                          legend = 'bottom',
                          align = 'hv')
p3.3


### 3.4) NUTS 2 regions ----
cols <- c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Set2'), brewer.pal(3, 'Accent')[1])

p3.4 <- ggplot(shapef) + 
  theme_bw(base_size = 15) + ggtitle('French NUTS 2 regions') + 
  geom_sf(aes(geometry = geometry, fill = NUTS_ID, group = NUTS_ID), show.legend = FALSE) + 
  geom_sf_text(aes(geometry = geometry,  label = NUTS_ID), size = 3,
               fun.geometry = sf::st_centroid) + 
  scale_fill_manual(values = cols) + xlab('') + ylab('')

p3.4

### 3.5) Hot- and cold-week indicator ----
pal <- colorRampPalette(brewer.pal(9,"Blues")[-1])(100)
borders <- shapef %>% summarize('geom1' = st_union(geometry))
coords <- data.frame('long' = c(-4,9), 'lat' = c(41.5, 51.3)) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = 3035) %>%
  st_coordinates()

p3.5.1 <- ggplot(df %>% dplyr::filter(Date == '2022-08-01')) + 
  theme_bw(base_size = 15) +
  ggtitle(bquote('Hot-week index (2022-08-01)')) + 
  geom_sf(aes(geometry = geometry, fill = w_avg_Tind95), col = 'gray80') + 
  geom_sf(data = borders$geom1, aes(), col = 'gray80', fill = NA) + 
  scale_fill_viridis_c(option = 'rocket', begin = 0.2, end = 1, direction = -1,
                       name = '', guide = guide_colorbar(barwidth = 10)) +  
  theme(legend.title = element_text(margin = margin(t = -0.75, unit = "lines"),
                                    size = 15),
        legend.position = 'bottom')  + 
  coord_sf(xlim = coords[,'X'], ylim = coords[,'Y']) 
p3.5.1

p3.5.2 <- ggplot(df %>% dplyr::filter(Date == '2013-01-14')) + 
  theme_bw(base_size = 15) +
  ggtitle(bquote('Cold-week index (2013-01-14)')) + 
  geom_sf(aes(geometry = geometry, fill = w_avg_Tind5), col = 'gray80') + 
  geom_sf(data = borders$geom1, aes(), col = 'gray80', fill = NA) + 
  scale_fill_gradientn(colours = pal, name = '', limits = NULL,
                       guide = guide_colorbar(barwidth = 10)) + 
  theme(legend.title = element_text(margin = margin(t = -0.75, unit = "lines"),
                                    size = 15),
        legend.position = 'bottom')  + 
  coord_sf(xlim = coords[,'X'], ylim = coords[,'Y']) 
p3.5.2

### 3.6) Influenza activity ----
REG  <- c('FR10', 'FRL0', 'FRH0')
NAME <- shapef$NUTS_NAME[which(shapef$NUTS_ID %in% REG)]

p3.6 <- ggplot(df %>% dplyr::filter(Region %in% REG, Age == 'Y_GE90')) + 
  geom_line(aes(x = Date, y = ia100 + avg_ia100, group = Region, colour = Region)) + 
  theme_bw(base_size = 15) + ylab('Incidence rate (per 100 inh.)') + 
  scale_color_manual(values = c(ablue, ared, agreen)) +
  ggtitle('Influenza-like illness') + 
  theme(legend.position = 'bottom')
p3.6

### 3.7) New hospitalizations for COVID-19 ----
REG  <- c('FR10', 'FRL0', 'FRH0')

p3.7 <- ggplot(df %>% dplyr::filter(Region %in% REG, Age == 'Y_GE90')) + 
  geom_line(aes(x = Date, y = Nhosp, group = Region, colour = Region)) + 
  theme_bw(base_size = 15) + ylab('Hospitalizations (per 1 000 inh.)') + 
  scale_color_manual(values = c(ablue, ared, agreen)) +
  ggtitle('Hospital data related to COVID-19') + 
  theme(legend.position = 'bottom') + 
  coord_cartesian(xlim = c(as.Date('2020-03-10'),as.Date('2023-02-28'))) 
p3.7


