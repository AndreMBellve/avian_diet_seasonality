#This script processes the species ranges maps into individual raster files (./output/season_rasters/) and produces a .csv of the species that have temporally variable ranges for us to reference individual months for their seasonal ranges.

#Author: André M. Bellvé

# Libraries ---------------------------------------------------------------
#Time keeping and notifications
library(tictoc)
library(beepr)
library(janitor)

#Basic data manipulation
library(dplyr)
library(dtplyr)
library(stringr)

#Visualisation
library(ggplot2)
library(viridis)

#Geospatial
library(sf)
library(terra)
library(tidyterra)


# Reading in data ---------------------------------------------------------

#Reading in all the data and selecting only relevant columns to manage model size
tic("All")#benchmarking - ~ 820 seconds for all features. PC showed 20 Gb of data
range_sf <- st_read("./data/range_maps/BOTW.gdb",
                    options = "METHOD=SKIP",
                    query = "SELECT SISID, SCI_NAME, PRESENCE, ORIGIN, SEASONAL FROM \"All_Species\"")

toc()
beep()

#Vector of species names from the original sf object to loop over
spp_vec <- unique(range_sf$SCI_NAME)

#Writing individual shapefiles because vect() seems to be having an issue with converting the entire sf range file
tic("Writing shapefiles")
for(i in seq_along(spp_vec)){
  
  sp_name <- str_replace(i, " ", "_")
  
  #Reading in the shape files using terra
  range_sf %>%
    filter(SCI_NAME == spp_vec[i]) %>% 
    #Converting them to a raster with a given target
    st_write(paste0("./output/shapefiles/",
                    spp_vec[i], ".shp"),
             append = FALSE, 
             quite = TRUE)
  
  print(paste0(i, "/", length(spp_vec)))
}
toc()
beep()

#2589.29 seconds to write shapefiles


# Rasterize ---------------------------------------------------------------
#The model raster will need to be masked, so that cells in the ocean have the value NA and otherwise the cells are 0 to indicate the species is absent in the cells on land
target_raster <- rast(extent = c(-180, 180, -90, 90), 
                      res = 0.9) 

tic("rasterising")
for(i in seq_along(spp_vec)){
  
  sp_name <- str_replace(i, " ", "_")
  
  #Reading in the shape files using terra
  vect(paste0("./output/shapefiles/",
                 spp_vec[i], ".shp")) %>%
    #Converting them to a raster with a given target
    rasterize(target_raster,
              field = "SEASONAL",
              touches = FALSE,
              fun = "min",
              background = 0,
              overwrite = TRUE,
              filename = paste0("./output/season_rasters/",
                                spp_vec[i], ".tif"))
  
  print(paste0(i, "/", length(spp_vec)))
}
toc()
beep()

# Loading prepped raster data ----------------------------------------------------
#Reading in species rasters
spp_rast <- list.files("./output/season_rasters/",
                       pattern = ".tif",
                       full.names = TRUE) %>% 
  rast()

#Assigning species names to each raster
names(spp_rast) <- list.files("./output/season_rasters/",
           pattern = ".tif") %>% 
  str_remove(".tif") %>%
  str_replace(" ", "_")

# Raster summaries --------------------------------------------------------
#Creating latitudinal masking layers to check the distribution of points in the N & S hemispheres
lat_rast <- init(spp_rast[[1]], "y")

north_mask <- ifel(lat_rast > 30, TRUE, FALSE)
south_mask <- ifel(lat_rast < -30, TRUE, FALSE)

#Masking rasters to be just each hemisphere and a 60 degree band around the equator
n_hemi_rast <- mask(spp_rast,
                    mask = north_mask,
                    maskvalues = FALSE)

s_hemi_rast <- mask(spp_rast,
                    mask = south_mask,
                    maskvalues = FALSE)

#Creating an "equatorial band as per the marious paper
equator <- ifel(lat_rast >= -30 & lat_rast <= 30, 
                     TRUE, FALSE)

equator_rast <- mask(spp_rast,
                     mask = equator,
                     maskvalues = FALSE)


#150 seconds
#Initialising storage lists
season_n_ls <- list()
season_s_ls <- list()
season_e_ls <- list()

n_hemi_mat <- as.matrix(n_hemi_rast)
s_hemi_mat <- as.matrix(s_hemi_rast)
equator_mat <- as.matrix(equator_rast)

#Loop to calculate various summary statistics
tic("season_count")
for(i in 1:nlyr(spp_rast)){
  
  #Pulling vector for each matrix the unique values by hemisphere and by equator
  season_n_ls[[i]] <- unique(n_hemi_mat[,i])
  
  season_s_ls[[i]] <- unique(s_hemi_mat[,i])
  
  season_e_ls[[i]] <- unique(equator_mat[,i])
  
}
toc()
beep()

#Creating dataframe to store values
spp_season_df <- data.frame(species = names(spp_rast), 
                            n_hemi_r = NA,
                            n_hemi_b = NA,
                            n_hemi_n = NA,
                            n_hemi_p = NA,
                            n_present = NA,
                            n_count = NA,
                            
                            s_hemi_r = NA,
                            s_hemi_b = NA,
                            s_hemi_n = NA,
                            s_hemi_p = NA,
                            s_present = NA,
                            s_count = NA,
                            
                            equa_r = NA,
                            equa_b = NA,
                            equa_n = NA,
                            equa_p = NA,
                            e_present = NA,
                            e_count = NA,
                            
                            resident_only = NA)

for(i in seq_along(season_e_ls)){
  
  spp_season_df[i,]$n_hemi_r <- ifelse(1 %in% season_n_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$n_hemi_b <- ifelse(2 %in% season_n_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$n_hemi_n <- ifelse(3 %in% season_n_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$n_hemi_p <- ifelse(4 %in% season_n_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$n_present <- ifelse(any(c(spp_season_df[i,]$n_hemi_r,
                                              spp_season_df[i,]$n_hemi_b,
                                              spp_season_df[i,]$n_hemi_n,
                                              spp_season_df[i,]$n_hemi_p)),
                                        TRUE, FALSE)
  spp_season_df[i,]$n_count <- sum(c(spp_season_df[i,]$n_hemi_r,
                                      spp_season_df[i,]$n_hemi_b,
                                      spp_season_df[i,]$n_hemi_n,
                                      spp_season_df[i,]$n_hemi_p))
  
  
  spp_season_df[i,]$s_hemi_r <- ifelse(1 %in% season_s_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$s_hemi_b <- ifelse(2 %in% season_s_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$s_hemi_n <- ifelse(3 %in% season_s_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$s_hemi_p <- ifelse(4 %in% season_s_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$s_present <- ifelse(any(c(spp_season_df[i,]$s_hemi_r,
                                              spp_season_df[i,]$s_hemi_b,
                                              spp_season_df[i,]$s_hemi_n,
                                              spp_season_df[i,]$s_hemi_p)),
                                        TRUE, FALSE)
  spp_season_df[i,]$s_count <- sum(c(spp_season_df[i,]$s_hemi_r,
                                      spp_season_df[i,]$s_hemi_b,
                                      spp_season_df[i,]$s_hemi_n,
                                      spp_season_df[i,]$s_hemi_p))
  
  spp_season_df[i,]$equa_r <- ifelse(1 %in% season_e_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$equa_b <- ifelse(2 %in% season_e_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$equa_n <- ifelse(3 %in% season_e_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$equa_p <- ifelse(4 %in% season_e_ls[[i]], TRUE, FALSE)
  spp_season_df[i,]$e_present <- ifelse(any(c(spp_season_df[i,]$equa_r,
                                              spp_season_df[i,]$equa_b,
                                              spp_season_df[i,]$equa_n,
                                              spp_season_df[i,]$equa_p)),
                                        TRUE, FALSE)
  spp_season_df[i,]$e_count <- sum(c(spp_season_df[i,]$equa_r,
                                     spp_season_df[i,]$equa_b,
                                     spp_season_df[i,]$equa_n,
                                     spp_season_df[i,]$equa_p))
  
  spp_season_df[i,]$resident_only <- ifelse(!any(spp_season_df[i,]$n_hemi_n, 
                                                 spp_season_df[i,]$n_hemi_b, 
                                                 spp_season_df[i,]$n_hemi_p,
                                                 
                                                 spp_season_df[i,]$s_hemi_n, 
                                                 spp_season_df[i,]$s_hemi_b, 
                                                 spp_season_df[i,]$s_hemi_p), 
                                            TRUE,
                                            FALSE)
}

spp_season_df <- spp_season_df %>% 
  mutate(n_hemi_only = ifelse(n_present & !e_present & !s_present,
                              TRUE, FALSE),
         equa_only = ifelse(!n_present & e_present & !s_present,
                            TRUE, FALSE),
         s_hemi_only = ifelse(!n_present & !e_present & s_present,
                              TRUE, FALSE),
         high_lat_only = ifelse(n_present & !e_present & s_present,
                                TRUE, FALSE),
         global = ifelse(n_present & e_present & s_present,
                         TRUE, FALSE))

#Saving the output file
write.csv(spp_season_df, 
          "./output/species_range_summary_incorrect_names.csv",
          row.names = FALSE)
#There was several names that had inconsistencies with other databases that were manually corrected - "./output/species_range_summary.csv"

#Creating summary to check the data
spp_season_summ <- spp_season_df %>% 
  summarise(north_only = sum(n_hemi_only),
            equa_only = sum(equa_only),
            south_only = sum(s_hemi_only),
            high_lat_only = sum(high_lat_only),
            global_spp = sum(global))

beep()

#Summary checks
spp_season_df %>%
  filter(!resident_only) %>% 
  summarise(n_breeding_ct =  sum(n_hemi_b),
            s_breeding_ct =  sum(s_hemi_b))

#Summary checks
spp_season_summ <- spp_season_df %>% 
  filter(!resident_only) %>% 
  summarise(north_only = sum(n_hemi_only),
            equa_only = sum(equa_only),
            south_only = sum(s_hemi_only),
            high_lat_only = sum(high_lat_only),
            global_spp = sum(global))

#Summary checks
spp_season_df %>% 
  filter(!resident_only) %>% 
  summarise(north_b = sum(n_hemi_b),
            equa_b = sum(equa_b),
            south_b = sum(s_hemi_b))

#Creating a df of species that have a non-resident range and saving a copy of it
spp_season_df %>%
  filter(!resident_only) %>% 
  select(species, n_present, s_present, e_present) %>% 
  write.csv("./output/non_resident_species.csv")
#We manually found data for each of these species on BOTW and then listed which data (resident, breeding, non-breeding, passage) was used for each month. 

#Loop to calculate various summary statistics
tic("season_count")
for(i in seq_along(nlyr(spp_rast))){
  #Counting the N hemisphere seasons
  season_n_ls[[i]] <- unique(n_hemi_rast[[i]])
  season_s_ls[[i]] <- unique(s_hemi_rast[[i]])
  season_e_ls[[i]] <- unique(equator_rast[[i]])
  
}
toc()
beep()

names(season_n_ls) <- names(spp_rast)
names(season_s_ls) <- names(spp_rast)
#Correcting counts and determining presences

#Initialising dataframe to fill with species counts and distributions
spp_dist_df <- data.frame(species = names(spp_rast),
                          north_seasons_n = NA,
                          south_seasons_n = NA)



################################################################################
# This script contains code for calculating predictability for temperature, precipitation, and gross primary productivity (gpp).

# Author: Reymond J. Miyajima

# libraries ---------------------------------------------------------------

# Geospatial
library(terra)
library(tidyterra)

# Data manipulation
library(tidyverse)
library(reshape2)

# Stats
library(WaveletComp)

# Misc
library(tictoc)
library(beepr)

# Data --------------------------------------------------------------------

# Load in Chelsea monthly precipitation and temperature data

# Temperature
temp_rast <- list.files("./data/tas",
                        pattern = ".tif",
                        full.names = TRUE) %>% 
  rast() 

# Clean up the names of the raster layers
clean_temp_names <- names(temp_rast) %>%
  str_remove("CHELSA_tas_") %>%
  str_remove("_V.2.1")

# Reassign raster layer names
names(temp_rast) <- clean_temp_names
rm(clean_temp_names)

# Precipitation
precip_rast <- list.files("./data/pet",
                          pattern = ".tif",
                          full.names = TRUE) %>% 
  rast() 

# Clean up the names of the raster layers
clean_precip_names <- names(precip_rast) %>%
  str_remove("CHELSA_pet_penman_") %>%
  str_remove("_V.2.1")

# Reassign raster layer names
names(precip_rast) <- clean_precip_names
rm(clean_precip_names)

# gpp layer

gpp_rast <- rast("./data/gpp_mean_5000_8_day.tif")

# Clean up the names of the raster layers
gpp_rast_names <- names(gpp_rast)%>% 
  str_remove("_Gpp")

# Reassign raster layer names
names(gpp_rast) <- gpp_rast_names
rm(gpp_rast_names)

# Read in model raster/masking layer
mask_lyr <- rast("./data/ts_mask_ts1_pc1.tif")


# Raster prep -------------------------------------------------------------

# These data are at a 1km resolution and need to be coarsened to 100km to match our biodiversity data

# Create a list to store new temperature rasters in
temp_list <- list() 

tic();for (i in 1:nlyr(temp_rast)) {
  
  # Since our world raster has a different spatial extent and resolution than the original Chelsa rasters, we need to aggregate the Chelsa rasters to a 100 km resolution
  
  temp_list[[i]] <- aggregate(temp_rast[[i]],
                              fact = 100,
                              fun = 'mean',
                              na.rm = TRUE) %>%
    project(mask_lyr)
  
  # Make a stack of the rasters in the list
  new_temp_rast <- rast(temp_list)
  
};beep();toc() # takes about ~2 hrs

# mask oceans 
new_temp_rast_masked <- mask(new_temp_rast,mask_lyr)

# Save this new raster because I don't want to wait 2 hrs again if my computer crashes
# writeRaster(new_temp_rast_masked,
#             "./data/new_temp_rast.tif",
#             overwrite = TRUE)

new_temp_rast <-  rast("./data/new_temp_rast.tif")

# Repeat this for precipitation data

# create a list to store new precipitation rasters in
precip_list <- list()

tic(); for (i in 1:nlyr(precip_rast)) {
  
  precip_list[[i]] <- aggregate(precip_rast[[i]],
                                fact = 100,
                                fun = 'mean',
                                na.rm = TRUE) %>%
    project(mask_lyr)
  
  new_precip_rast <- rast(precip_list)
}; beep(); toc() 
# Takes ~ 2 hrs

# plot check
plot(new_precip_rast_masked[[1]])

# Save new raster
# writeRaster(new_precip_rast,
#             "./data/new_precip_rast.tif",
#             overwrite = TRUE)

new_precip_rast <-  rast("./data/new_precip_rast.tif")

# GPP layer has the correct spatial dimensions however, it's at an a day interval. Aggerate site data to a monthly resolution. 


# time series goes from 2021 to 2023
months <- c(paste0("2021_",sprintf("%02d", 1:12)),
            paste0("2022_",sprintf("%02d", 1:12)),
            paste0("2023_",sprintf("%02d", 1:12)))

# convert raster into a dataframe.
gpp_df <- terra::as.data.frame(gpp_rast, xy = TRUE)
latlong <- gpp_df %>% 
  select(x,
         y)

# store each month for a given year into a list
month_ls <- list()

for (i in months) {
  month_i <- gpp_df %>% 
    select(starts_with(i))
  month_ls[[i]] <- rowMeans(month_i)
}

monthly_gpp <- do.call(cbind, month_ls) %>% 
  cbind(latlong)

# Climate data prep -------------------------------------------------------

# Convert rasters into a dataframe
temp_df <- terra::as.data.frame(new_temp_rast,
                                xy = FALSE)


# Get the years 
years <- seq(1980,2019,by=1)

year_list <- list()

for (i in years) {
  
  year <- temp_df %>% 
    select(contains(paste0("_",i)))
  
  year_list[[as.character(i)]] <- year
  
}

names(year_list) <- NULL

temp_df_final <- do.call(cbind, year_list) 

# Save data
saveRDS(temp_df_final, "./data/temp_df_final.rds")


# Do this for precipitation 
precip_df <- terra::as.data.frame(new_precip_rast,
                                  xy = FALSE) 


# Get the years 
years <- seq(1980,2018,by=1)

year_list <- list()

for (i in years) {
  
  year <- precip_df %>% 
    select(contains(paste0("_",i)))
  
  year_list[[as.character(i)]] <- year
  
}

names(year_list) <- NULL

precip_df_final <- do.call(cbind, year_list) 

# Save data
saveRDS(precip_df_final, "./data/precip_df_final.rds")


# Do this for GPP

# Get the years 
years <- seq(2021,2023,by=1)

year_list <- list()

for (i in years) {
  
  year <- monthly_gpp %>% 
    select(contains(paste0(i)))
  
  year_list[[as.character(i)]] <- year
  
}

names(year_list) <- NULL

gpp_final <- do.call(cbind, year_list) 

# convert nas to zero.
gpp_df_final <- gpp_final %>% 
  mutate_all(~replace(., is.na(.), 0))

# Save data
saveRDS(gpp_df_final, "./data/gpp_df_final.rds")

# Predictability calculations ---------------------------------------------

# Predictability is the reliable recurrence of a phenomena over time. For example, it's always hot in the jungle in June every year.

# Because we do not have daily time-series data, we cannot use Colwell's indices, and therefore we conducted a wavelet analysis. This is a more flexible approach as it will allow us to quantify predictability at any time scale. Here, following Tonkin et al 2017, we define predictability as the average proportion of wavelet power that is significant at the 12-month interval across the entire time series. More simply, predictability is how much of the total variability in the data is explained by the recurring 12-month pattern. 

# create list for grid cell sites
site_list <- list()

# calulate predictability for each site
for (i in 1:nrow(temp_df_final)) {
  
  # Pull out each site one by one
  site <- temp_df_final[i,] 
  
  # Transform data for wavelet analysis
  site <- site %>%
    melt()%>%
    rename('date' = 'variable',
           'temp' = 'value')
  
  # calculate predictability 
  site_wv <- analyze.wavelet(site,
                             loess.span = 0, # There is no need to detrend this series
                             my.series = 2,
                             dt = 1,
                             method = 'white.noise')
  
  site_predict <- site_wv[["Power.avg"]]
  
  site_list[[i]] <- site_predict[12] 
  
}

site_predict_temp_mat <- do.call(rbind,site_list)

# repeat this for precipitation  
site_list <- list()

for (i in 1:nrow(precip_df_final)) {
  
  # Pull out each site one by one
  site <- precip_df_final[i,] 
  
  # Transform data for wavelet analysis
  site <- site %>%
    melt()%>%
    rename('date' = 'variable',
           'temp' = 'value')
  
  # calculate predictability 
  site_wv <- analyze.wavelet(site,
                             loess.span = 0, # There is no need to detrend this series
                             my.series = 2,
                             dt = 1,
                             method = 'white.noise')
  
  site_predict <- site_wv[["Power.avg"]]
  
  site_list[[i]] <- site_predict[12] 
  
}

site_predict_precip_mat <- do.call(rbind,site_list)

# repeat this for gpp 

site_list <- list()

for (i in 1:nrow(gpp_df_final)) {
  
  # Pull out each site one by one
  site <- gpp_df_final[i,] 
  
  # Transform data for wavelet analysis
  site <- site %>%
    melt()%>%
    rename('date' = 'variable',
           'temp' = 'value')
  
  # calculate predictability 
  site_wv <- analyze.wavelet(site,
                             loess.span = 0, # There is no need to detrend this series
                             my.series = 2,
                             dt = 1,
                             method = 'white.noise')
  
  site_predict <- site_wv[["Power.avg"]]
  
  site_list[[i]] <- site_predict[12] 
  
}

site_predict_gpp_mat <- do.call(rbind,site_list)


# Creating layers ------------------------------------------

# Get lat long coordinates for each temperature site
latlong <- terra::as.data.frame(new_temp_rast,
                                xy = TRUE) %>% 
  select(x,
         y) %>% 
  rename(long = x,
         lat = y) %>% 
  as.matrix()

# read in temperature predictability data
temp_predict_mat <- readRDS('./data/site_predict_temp_mat.rds')

# add lat long coordinates and rasterize and reproject to have same extent as our original model raster
temp_predict_rast <- temp_predict_mat %>% 
  cbind(latlong) %>% 
  as.data.frame() %>% 
  rename(temp_predict = 'V1') %>% 
  as_spatraster(xycols = 2:3,
                crs = 'EPSG:4326') %>% 
  project(mask_lyr)

# mask
temp_predict_rast_masked <- mask(temp_predict_rast,
                                 mask_lyr)

plot(temp_predict_rast_masked)

# Repeat for precipitation
latlong <- terra::as.data.frame(new_precip_rast,
                                xy = TRUE) %>% 
  select(x,
         y) %>% 
  rename(long = x,
         lat = y) %>% 
  as.matrix()

precip_predict_mat <- readRDS('./data/site_predict_precip_mat.rds')

precip_predict_rast <- precip_predict_mat %>% 
  cbind(latlong) %>% 
  as.data.frame() %>% 
  rename(precip_predict = 'V1') %>% 
  as_spatraster(xycols = 2:3,
                crs = 'EPSG:4326') %>% 
  project(mask_lyr)

# mask
precip_predict_rast_masked <- mask(precip_predict_rast,
                                   mask_lyr)

plot(precip_predict_rast_masked)


# gpp
gpp_predict_mat <- readRDS('./data/site_predict_gpp_mat.rds')

gpp_predict_rast <- gpp_predict_mat %>% 
  cbind(latlong) %>% 
  as.data.frame() %>% 
  rename(gpp_predict = '.') %>% 
  as_spatraster(xycols = 2:3,
                crs = 'EPSG:4326') %>% 
  project(mask_lyr)

# mask
gpp_predict_rast_masked <- mask(gpp_predict_rast,
                                mask_lyr)

plot(gpp_predict_rast_masked)

# save outputs.
writeRaster(temp_predict_rast_masked,
            "./data/temp_predict_rast_masked.tif",
            overwrite = TRUE)
writeRaster(precip_predict_rast_masked,
            "./data/precip_predict_rast_masked.tif",
            overwrite = TRUE)
writeRaster(gpp_predict_rast_masked,
            "./data/gpp_predict_rast_masked.tif",
            overwrite = TRUE)



