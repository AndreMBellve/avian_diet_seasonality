#Regression of assemblage diet variability by environmental covariates of seasonality and it's predictability.

#Author: André M. Bellvé

# Libraries ---------------------------------------------------------------
#Data manipulation
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(matrixStats)

#Geospatial processing
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(geodata)

#Visualisation
library(ggplot2)
library(khroma)
library(patchwork)

#Modelling
library(lqr)

# Data prep ---------------------------------------------------------------

#Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
pca_geography_matrix <- readRDS("./output/pca_geography_matrix.rds")

#geography pca has 80000*6 chunks The first 80000 is pc1
ms_pca_mat_t <- t(ms_pca_matrix)

#Dataframe of row indices to iterate through
row_index_df <- data.frame(begin = c(1, 80001, 160001,
                                     240001, 320001, 400001),
                           end = c(80000, 160000, 240000, 
                                   320000, 400000, 480000))

#Iterating through to build the rasters for each of our original six PCA axes - these are split into  80k chunks
for(pc in 1:6){
  
  cur_rast_vals <- ms_pca_mat_t[seq(from = row_index_df[pc,"begin"],
                                    to = row_index_df[pc,"end"],
                                    by = 1),]
  
  #load in raster and make a model raster
  model_rast <- rast("./output/model_raster.tif")
  
  #create 3 replicate layers 
  pc_rast <- rep(model_rast, 2)
  
  #create a list to store matrix monthly PC data
  pc_list <- list()
  
  for(c in 1:ncol(cur_rast_vals)){
    
    pc_list[[c]] <- cur_rast_vals[,c] #add data to list
    
    #change raster values to the pc values for the corresponding 
    values(pc_rast[[c]]) <- pc_list[[c]]
    
  }
  
  #relabel names
  pc_names <- paste0(c('TS1','TS2'), "-PC", pc)
  
  #Naming the rasters for clarity
  names(pc_rast) <- pc_names
  
  #Either initialising a spatraster storage object or 
  if(pc == 1){
    ts_rast <- pc_rast
  }else{
    ts_rast <- c(ts_rast, pc_rast)
  }
}

ts_rast

#Downloading the bioclim data
bioclim_data <- worldclim_global(var = "bio",
                                 res = 5,
                                 path = "./data/env_rasters/")
#Renaming for ease
names(bioclim_data) <- names(bioclim_data) %>% 
  str_remove("wc2.1_5m_")

#Selecting seasonality variables and renaming
season_rast <- bioclim_data[[c("bio_4", "bio_15")]]
names(season_rast) <- c("temp_seasonality", "precip_seasonality")

season_rast <- project(season_rast,
                       ts_rast)

#Reprojecting the seasonality rasters to match the dimensions of the ts_rast
ts_mask_rast <- mask(ts_rast,
                     season_rast)

#Writing a model raster to upload to GEE to reproject the MODIS GPP layers
writeRaster(ts_mask_rast[[1]],
            "./output/ts_pca_raster/ts_mask_ts1_pc1.tif",
            filetype = "GTiff",
            overwrite = TRUE)

#This one has all the zones outside of bird ranges masked, which the other model raster did not
model_rast <- rast("./output/ts_pca_raster/ts_mask_ts1_pc1.tif")


# GPP CV creation ---------------------------------------------------------
#Downloaded the GEE summarised 8 day resolution GPP means (based on 5000 cells in each 100 × 100 km pixel)
gpp_rast <- rast("./data/env_rasters/modis_gpp_gee/gpp_mean_5000_8_day.tif")

#The three years of data that we aggregate over
year_vec <- c("2021", "2022", "2023")

#For loop to iterate through the 8 day layer for each year and calculate the CV
for(i in seq_along(year_vec)){
  
  gpp_yr_rast <- tidyterra::select(gpp_rast, 
                                   starts_with(year_vec[i]))
  
  gpp_yr_sd_rast <- gpp_yr_rast %>% 
    app(fun = "sd",
        na.rm = TRUE)
  
  gpp_yr_mean_rast <-gpp_yr_rast %>% 
    app(fun = "mean",
        na.rm = TRUE)
  
  gpp_yr_cv_rast <- (gpp_yr_sd_rast / gpp_yr_mean_rast) * 100 
  
  names(gpp_yr_cv_rast) <- paste0(year_vec[i],"_gpp_cv")
  
  if(i == 1){
    gpp_cv_rast <- gpp_yr_cv_rast
  }else{
    gpp_cv_rast <- c(gpp_cv_rast, gpp_yr_cv_rast)
  }
}
plot(gpp_cv_rast)

#Summarising the mean CV across our three years of data
gpp_mean_cv_rast <- app(gpp_cv_rast,
                        fun = "mean",
                        na.rm = TRUE)

#Matching the dimension/projection of the season and ts rasters and overwriting the missing values with NAs
gpp_mean_cv_reproj_rast <- project(gpp_mean_cv_rast,
                                   ts_rast) %>% 
  #This one remove the values that the GPP layer has that the others do not
  mask(x = ., 
       mask = season_rast,
       updatevalue = NA)

#Renaming for clarity
names(gpp_mean_cv_reproj_rast[[1]]) <- "gpp_seasonality" 

#Adding this to the SpatRaster object for processing
season_rast <- c(season_rast, gpp_mean_cv_reproj_rast[[1]])

# Difference rasters ------------------------------------------------------

#Listing all the summary matrices to iterate through
summary_matrices <- list.files("./output/summary_matrices/",
                               full.names = TRUE,
                               pattern = "mean_")

#For loop to iteratively create rasters containing the max difference for each PC
for(i in seq_along(summary_matrices)){
  
  #Iteratively reading in the PC mean matrices
  pc_monthly_means <- readRDS(summary_matrices[i])
  
  #Creating a dataframe to pull from
  pc_diff_df <- data.frame(max = rowMaxs(pc_monthly_means, na.rm =  TRUE),
                           min = rowMins(pc_monthly_means, na.rm =  TRUE)) %>% 
    mutate(id = 1:nrow(.),
           pc_id = paste0("PC", i),
           diff = max - min,
           diff = replace(diff, is.infinite(diff), NA)) 
  
  #Initialising a raster with the correct dimensions/CRS
  pc_diff_rast <- ts_rast[[1]]
  
  #Making sure there are no errant values carried over
  values(pc_diff_rast) <- NA
  
  #Assinging the values
  values(pc_diff_rast) <- pc_diff_df$diff
  
  #Coherent naming
  names(pc_diff_rast) <- paste0("pc", i, "_diff")
  
  #Masking to ensure that we do not have any mismatched data
  pc_diff_rast <- mask(x = pc_diff_rast,
                       mask = season_rast[[1]]) %>% 
    mask(x = .,
         mask = ts_rast[[1]])
  
  #Initialising the list to store all of the PC diff rasters or adding to an existing one
  if(i == 1){
    all_diff_rast <- pc_diff_rast
  }else{
    all_diff_rast <- c(all_diff_rast, pc_diff_rast)
  }
  
  #Iteration tracker
  cat("\nPC", i, " completed",
      sep = "")
  
}

#Saving output
writeRaster(all_diff_rast,
            filename = paste0("./output/annual_diff_rasters/",
                              names(all_diff_rast),
                              ".tif"),
            overwrite = TRUE)

# Predictability rasters --------------------------------------------------

#Reading in our predictability rasters
season_pred_rast <- list.files("./data/env_rasters/predictability/",
                               full.names = TRUE,
                               pattern = ".tif") %>% 
  rast()

#Creating visualisations of predictability to inspect trends


# Get the continent and ocean shapefile
globe_shp <- ne_download(scale = "medium", 
                         type = "land",
                         category = "physical", 
                         returnclass = "sf") 
ocean_shp <- ne_download(scale = "medium", 
                         type = "ocean",
                         category = "physical", 
                         returnclass = "sf")

#Create a grid to overlay the earth
grid_full <- list(
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-180, -120, -60, 0, 60, 120, 180),
                                  function(x) cbind(x, seq(-90, 90, 1)))), 
    crs = 'WGS84')),
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-90, -60, -30, 0, 30, 60, 90), function(x) {
      cbind(seq(-180, 180, 1), x)
    })), crs = 'WGS84'))) %>% 
  bind_rows()

#Maps of predictability
ggplot() + 
  
    geom_sf(data = ocean_shp, 
          fill = "grey70",
          colour = "transparent",
          alpha = 1, 
          size = 0.2) +
  
  geom_sf(data = grid_full,
          linewidth = 0.2,
          colour = "grey40") +
  
  geom_spatraster(data = season_pred_rast,
                  maxcell = 80000) + 
  
  facet_wrap(~lyr,
             ncol = 1,
             labeller = as_labeller(c(gpp_predict = "GPP", 
                                      precip_predict = "Precipitation", 
                                      temp_predict = "Temperature"))) + 
  
  scale_fill_viridis_c(option = "viridis",
                       na.value = "transparent",
                       name = "Predictability\n",
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1)) + 
  
  coord_sf(crs = "ESRI:54030",
           expand = FALSE) +
  
  theme_void() + 
  
  theme(legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.key.width = unit(dev.size()[1] / 15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", 
                                       color = NA),
        panel.background = element_rect(fill = "transparent", 
                                        color = NA))


#Saving map for supplementary materials
ggsave("./graphs/seasonality_regression/predictability_maps.png",
       width = 12, height = 14.5)


#Maps of seasonality
#TEMPERATURE MAP
temp_map_gg <- ggplot() + 
  
  geom_sf(data = ocean_shp, 
          fill = "grey70",
          colour = "transparent",
          alpha = 1, 
          size = 0.2) +
  
  geom_sf(data = grid_full,
          linewidth = 0.2,
          colour = "grey40") +
  
  geom_spatraster(data = season_rast[["temp_seasonality"]],
                  maxcell = 80000) + 
  
  scale_fill_viridis_c(option = "magma",
                       na.value = "transparent",
                       name = "Temperature\nSeasonality\n") + 
  
  coord_sf(crs = "ESRI:54030",
           expand = FALSE) +
  
  theme_void() + 
  
  theme(legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.key.width = unit(dev.size()[1] / 15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", 
                                       color = NA),
        panel.background = element_rect(fill = "transparent", 
                                        color = NA));temp_map_gg

#PRECIPITATION MAP
precip_map_gg <- ggplot() + 
  
  geom_sf(data = ocean_shp, 
          fill = "grey70",
          colour = "transparent",
          alpha = 1, 
          size = 0.2) +
  
  geom_sf(data = grid_full,
          linewidth = 0.2,
          colour = "grey40") +
  
  geom_spatraster(data = season_rast[["precip_seasonality"]],
                  maxcell = 80000) + 
  
  scale_fill_viridis_c(option = "cividis",
                       na.value = "transparent",
                       name = "Precipitation\nSeasonality\n", 
                       direction = 1) + 
  
  coord_sf(crs = "ESRI:54030",
           expand = FALSE) +
  
  theme_void() + 
  
  theme(legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.key.width = unit(dev.size()[1] / 15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", 
                                       color = NA),
        panel.background = element_rect(fill = "transparent", 
                                        color = NA));precip_map_gg

#GPP MAP
gpp_map_gg <- ggplot() + 
  
  geom_sf(data = ocean_shp, 
          fill = "grey70",
          colour = "transparent",
          alpha = 1, 
          size = 0.2) +
  
  geom_sf(data = grid_full,
          linewidth = 0.2,
          colour = "grey40") +
  
  geom_spatraster(data = season_rast[["gpp_seasonality"]],
                  maxcell = 80000) + 
  
  scale_fill_viridis_c(option = "viridis",
                       na.value = "transparent",
                       name = "GPP\nSeasonality\n", 
                       direction = 1, 
                       breaks = c(0, 100, 200, 300),
                       limits = c(0,300)) + 
  
  coord_sf(crs = "ESRI:54030",
           expand = FALSE) +
  
  theme_void() + 
  
  theme(legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.key.width = unit(dev.size()[1] / 15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", 
                                       color = NA),
        panel.background = element_rect(fill = "transparent", 
                                        color = NA));gpp_map_gg


#Binding all the maps into one
gpp_map_gg + precip_map_gg + temp_map_gg + 
  plot_layout(ncol = 1)
ggsave("./graphs/seasonality_regression/seasonality_maps.png",
       width = 12, height = 14.5)

#Combining values to have all env variables in a single object
season_rast <- c(season_rast, season_pred_rast)

# Difference df -----------------------------------------------------------

#Reading in the variation raster as all diff to test the output
all_diff_rast <- list.files("./output/annual_diff_rasters/",
                            pattern = ".tif",
                            full.names = TRUE) %>% 
  rast()

#Creating a dataframe with all differences and the seasonality correlates
diff_season_df <- values(c(all_diff_rast, season_rast), 
                         dataframe = TRUE,
                         na.rm = FALSE) %>% 
  filter(!is.na(pc1_diff)) %>% 
  bind_cols(cell_id = cells(all_diff_rast)) %>% 
  clean_names() 

#Pulling the coordinates so I can calculate their latitudinal position
coords <- xyFromCell(all_diff_rast,
                     diff_season_df$cell_id) %>% 
  data.frame(latitude = .[,"y"],
             longitude = .[,"x"])

#Reading in the eigenvalues from our diet LRA to weight our scores by
eigenval_df <- readRDS("./output/diet_pca/diet_pca.rds")$sdev %>% 
  data.frame(eigen_prop = round((.^2 / sum(.^2)), digits = 4),
             PC = paste0("pc", 1:7)) %>% 
  dplyr::select(-c(`.`))

#Min-max feature scaling function
mm_feature_scale <- function(x){
  ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }

#Creating the difference dataframe and deriving the additional information for graphing and modelling
diff_season_long_df <- diff_season_df %>% 
  
  #Binding on the coordinates
  bind_cols(dplyr::select(coords, c(latitude, longitude))) %>% 
  
  #Creating latitudinal data
  mutate(
    #Latitudinal bands at 10 degree increments
    lat_band = round(latitude, digits = -1),
    lat_band = replace(lat_band, lat_band == -90, -80),
    lat_band = factor(as.character(lat_band), 
                      levels = c(as.character(seq(from = 80,
                                                  to = -80, 
                                                  by = -10)))),
    #Three band split
    triband = case_when(
      latitude >= 30 ~ "North [90° – 30°]",
      latitude < 30 & latitude > -30 ~ "Equatorial (30° – -30°)",
      latitude <= - 30 ~ "South [-30° – -90°]"),
    triband = factor(triband, levels = c("North [90° – 30°]", 
                                         "Equatorial (30° – -30°)", 
                                         "South [-30° – -90°]")),
    
     lat_30_band = case_when(
        latitude > 60 ~ "> 60 N",
        latitude > 30 & latitude <= 60 ~ "30 – 60 N",
        latitude > 0 & latitude <= 30 ~ "0 – 30 N",
        latitude < 0 & latitude >= -30 ~ "-30 – 0 S",
        latitude < -30 & latitude >= -60 ~ "-60 – -30 S",
        latitude < -60 & latitude >= -90 ~ "< -60 S"),
    lat_30_band = factor(lat_30_band,
                         levels = c("> 60 N", "30 – 60 N",
                                    "0 – 30 N", "-30 – 0 S", 
                                    "-60 – -30 S", "< -60 S")),
    
    #Hemisphere designation
    hemisphere = ifelse(latitude >= 0, "North", "South")) %>% 
  
  #Normalising all the PC data to be between 0 and 1 so that all values are commensurate using min-max feature scaling
  mutate(across(ends_with("_diff"), mm_feature_scale)) %>% 
  
  #Pivoting to long format to assign the eigenvalue proportions
  pivot_longer(cols = ends_with("_diff"),
               names_to = "PC",
               values_to = "diff") %>% 
  
  #Modifying names
  mutate(PC = str_remove(PC, "_diff")) %>% 
  
  #PC score summary
  left_join(eigenval_df,
            join_by("PC")) %>% 
    
  #Multiplying all PC's by the eigen proportion to have their weighting reflected
  mutate(eigen_wgt_diff = diff * eigen_prop,
         PC = toupper(PC)) 


# Summed Variation --------------------------------------------------------

#Summarising by summing the eigen-weighted differences of all the PCs
diff_season_summ_df <- diff_season_long_df %>% 
  group_by(cell_id, latitude, longitude, 
           lat_band, triband, hemisphere, lat_30_band,
           precip_seasonality, temp_seasonality, gpp_seasonality,
           precip_predict, temp_predict, gpp_predict) %>% 
  summarise(summed_wgt_diff = sum(eigen_wgt_diff))


# Summed diff raster ------------------------------------------------------

summed_diff_df <- data.frame(cell_id = 1:80000) %>% 
  left_join(diff_season_summ_df[,c("cell_id", "summed_wgt_diff")],
            by = "cell_id")

#Initialising a raster with the correct dimensions/CRS
summed_diff_rast <- ts_rast[[1]]

#Making sure there are no errant values carried over
values(summed_diff_rast) <- NA
values(summed_diff_rast) <- summed_diff_df$summed_wgt_diff
names(summed_diff_rast) <- "summed_diff"

# Get the continent and ocean shapefile
globe_shp <- ne_download(scale = "medium", 
                         type = "land",
                         category = "physical", 
                         returnclass = "sf") 
ocean_shp <- ne_download(scale = "medium", 
                         type = "ocean",
                         category = "physical", 
                         returnclass = "sf")

#Create a grid to overlay the earth
grid_full <- list(
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-180, -120, -60, 0, 60, 120, 180),
                                  function(x) cbind(x, seq(-90, 90, 1)))), 
    crs = 'WGS84')),
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-90, -60, -30, 0, 30, 60, 90), function(x) {
      cbind(seq(-180, 180, 1), x)
    })), crs = 'WGS84'))) %>% 
  bind_rows()

#Map of the summed eigen weighted change in diet variability
ggplot() + 
  
  geom_sf(data = ocean_shp, 
          fill = "grey70",
          colour = "transparent",
          alpha = 1, 
          size = 0.2) +
  
  geom_sf(data = grid_full,
          linewidth = 0.2,
          colour = "grey40") +
  
  geom_spatraster(data = summed_diff_rast,
                  maxcell = 80000) + 
  
  scale_fill_viridis_c(option = "magma",
                       na.value = "transparent",
                       name = "Assemblage\nDiet Variability  ",
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1)) + 

  coord_sf(crs = "ESRI:54030",
           expand = FALSE) +
  
  theme_void() + 
  
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "left",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.width = unit(dev.size()[1] / 15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", 
                                       color = NA),
        panel.background = element_rect(fill = "transparent", 
                                        color = NA))

ggsave("./graphs/diet_variability/summed_eigen_weighted_diet_variability.png",
       width = 16, height = 9.7)


# Quantile Regression -----------------------------------------------------

#GPP
gpp_0.5_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(gpp_seasonality, 2) * lat_30_band,
                                data = diff_season_summ_df %>% filter(!is.na(gpp_seasonality)) %>% droplevels(),
                                p = 0.5, b = 0.76)
#Saving model output
gpp_0.5_band_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/gpp_qr_coefficients.csv")

#Individual level models for quantile regression to avoid discarding missing data
#Precip
precip_0.5_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(precip_seasonality, 2) * lat_30_band,
                                data = diff_season_summ_df %>% filter(!is.na(precip_seasonality)),
                                p = 0.5, b = 0.76)
#Saving model output
precip_0.5_band_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/precip_qr_coefficients.csv")

#Temp
temp_0.5_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(temp_seasonality, 2) * lat_30_band,
                                data = diff_season_summ_df %>% filter(!is.na(temp_seasonality)),
                                p = 0.5, b = 0.76)
#Saving model output
temp_0.5_band_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/temp_qr_coefficients.csv")

# Graphing ----------------------------------------------------------------


## Summed variation --------------------------------------------------------

#Precip band quantile ribbons
precip_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(precip_seasonality, 2) * lat_30_band,
                            data = diff_season_summ_df %>% filter(!is.na(precip_seasonality)) %>% droplevels(),
                            p = c(0.05, 0.5, 0.95), b = 0.76)


precip_rqs <- bind_cols(precip_0.05 = precip_band_qr[[1]]$fitted.values,
                     precip_0.5 = precip_band_qr[[2]]$fitted.values,
                     precip_0.95 = precip_band_qr[[3]]$fitted.values,
                     diff_season_summ_df %>% 
                       filter(!is.na(precip_seasonality)))


precip_summed_diff_gg <- ggplot(data = diff_season_summ_df,
                                aes(y = summed_wgt_diff,
                                    x = precip_seasonality,
                                    colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_line(data = precip_rqs,
            aes(y = precip_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1)) +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(y = "Assemblage Diet Variability",
       x = "Precipitation Seasonality\n(CV)",
       linetype = "Latitudinal\nBand (\u00B0)") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),,
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"));precip_summed_diff_gg



#Temperature
temp_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(temp_seasonality, 2) * lat_30_band,
                            data = diff_season_summ_df %>% filter(!is.na(temp_seasonality)) %>% droplevels(),
                            p = c(0.05, 0.5, 0.95), b = 0.76)


temp_rqs <- bind_cols(temp_0.05 = temp_band_qr[[1]]$fitted.values,
                     temp_0.5 = temp_band_qr[[2]]$fitted.values,
                     temp_0.95 = temp_band_qr[[3]]$fitted.values,
                     diff_season_summ_df %>% 
                       filter(!is.na(temp_seasonality)))


temp_summed_diff_gg <-  ggplot(data = diff_season_summ_df,
                               aes(y = summed_wgt_diff,
                                   x =  temp_seasonality,
                                   colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_line(data = temp_rqs,
            aes(y = temp_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1)) +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(linetype = "Latitudinal\nBand (\u00B0)",
       x = "Temperature Seasonality\n(SD × 100)") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),,
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()); temp_summed_diff_gg


#GPP
gpp_band_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(gpp_seasonality, 2) * lat_30_band,
                                data = diff_season_summ_df %>% filter(!is.na(gpp_seasonality)) %>% droplevels(),
                                p = c(0.05, 0.5, 0.95), b = 0.76)


gpp_rqs <- bind_cols(gpp_0.05 = gpp_band_qr[[1]]$fitted.values,
                     gpp_0.5 = gpp_band_qr[[2]]$fitted.values,
                     gpp_0.95 = gpp_band_qr[[3]]$fitted.values,
                     diff_season_summ_df %>% 
                       filter(!is.na(gpp_seasonality)))


gpp_summed_diff_gg <- ggplot(data = diff_season_summ_df,
                             aes(y = summed_wgt_diff,
                                 x = gpp_seasonality,
                                 colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_line(data = gpp_rqs,
            aes(y = gpp_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1),
                        guide = "none") +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +

  scale_x_continuous(limits = c(0, 300)) +
  
  labs(y = "Assemblage Diet Variability",
       x = "GPP Seasonality\n(CV)",
       linetype = "Latitudinal\nBand (\u00B0)") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()); gpp_summed_diff_gg


triband_hemi_summed_gg <- precip_summed_diff_gg +
  temp_summed_diff_gg +
  gpp_summed_diff_gg +
  plot_layout(guides = "collect")

triband_hemi_summed_gg

ggsave("./graphs/seasonality_regression/three_pane_seasonality_summed_diff.png",
       width = 22.5, height = 17)



## SM quantile ribbon -----------------------------------------------------

#Precip band quantile ribbons
precip_summed_diff_gg <- ggplot(data = diff_season_summ_df,
                                aes(y = summed_wgt_diff,
                                    x = precip_seasonality,
                                    colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +

  geom_ribbon(data = precip_rqs,
              aes(ymin = precip_0.05,
                  ymax = precip_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "transparent",
              fill = "grey",
              alpha = 0.5,
              linewidth = 1) +
  
  geom_ribbon(data = precip_rqs,
              aes(ymin = precip_0.05,
                  ymax = precip_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "#111111",
              fill = "transparent",
              linewidth = 1) +
  
  geom_line(data = precip_rqs,
            aes(y = precip_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
                colour = "black",
                linewidth = 1.5) +
  
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1)) +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(y = "Assemblage Diet Variability",
       x = "Precipitation Seasonality\n(CV)",
       linetype = "Latitudinal\nBand (\u00B0)") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),,
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"));precip_summed_diff_gg



#Temperature
temp_summed_diff_gg <-  ggplot(data = diff_season_summ_df,
                               aes(y = summed_wgt_diff,
                                   x =  temp_seasonality,
                                   colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_ribbon(data = temp_rqs,
              aes(ymin = temp_0.05,
                  ymax = temp_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "transparent",
              fill = "grey",
              alpha = 0.5,
              size = 1) +
  
  geom_ribbon(data = temp_rqs,
              aes(ymin = temp_0.05,
                  ymax = temp_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "#111111",
              fill = "transparent",
              size = 1) +
  
  geom_line(data = temp_rqs,
            aes(y = temp_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1)) +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(linetype = "Latitudinal\nBand (\u00B0)",
       x = "Temperature Seasonality\n(SD × 100)",
       y = "") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),,
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank()); temp_summed_diff_gg


#GPP
gpp_summed_diff_gg <- ggplot(data = diff_season_summ_df,
                             aes(y = summed_wgt_diff,
                                 x = gpp_seasonality,
                                 colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_ribbon(data = gpp_rqs,
              aes(ymin = gpp_0.05,
                  ymax = gpp_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "transparent",
              fill = "grey",
              alpha = 0.5,
              size = 1) +
  
  geom_ribbon(data = gpp_rqs,
              aes(ymin = gpp_0.05,
                  ymax = gpp_0.95,
                  group = lat_30_band,
                  linetype = lat_30_band),
              colour = "#111111",
              fill = "transparent",
              size = 1) +
  
  geom_line(data = gpp_rqs,
            aes(y = gpp_0.5,
                group = lat_30_band,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  facet_wrap(~triband,
             ncol = 1) + 
  
  scale_fill_nightfall(name = "Latitude (\u00B0)",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude (\u00B0)",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1, 4, 5, 
                                   3, 4, 1),
                        guide = "none") +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  scale_x_continuous(limits = c(0, 300)) +
  
  labs(y = "",
       x = "GPP Seasonality\n(CV)",
       linetype = "Latitudinal\nBand (\u00B0)") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 20),
        axis.text = element_text(colour = "black",size = 20),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()); gpp_summed_diff_gg

#Patching the plots together
triband_hemi_summed_gg <- precip_summed_diff_gg +
  temp_summed_diff_gg +
  gpp_summed_diff_gg +
  plot_layout(guides = "collect")

triband_hemi_summed_gg

ggsave("./graphs/seasonality_regression/three_pane_seasonality_summed_diff_ribbon.png",
       width = 22.5, height = 17)


# Precictability ----------------------------------------------------------

#Precip pred
#All data
precip_pred_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(precip_predict, 3),
                                    b = 0.76,
                                    data = diff_season_summ_df %>% filter(!is.na(precip_predict)) %>% droplevels(),
                                    p = 0.5)

#Saving model output
precip_pred_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/precip_pred_qr_coefficients.csv")

#Pulling fitted values for all of the points
precip_pred_rqs <- bind_cols(precip_0.5 = precip_pred_qr$fitted.values,,
                             diff_season_summ_df %>% 
                               filter(!is.na(precip_predict)) %>% droplevels())

#Graph
precip_pred_summed_diff_gg <- ggplot(data = diff_season_summ_df,
                                aes(y = summed_wgt_diff,
                                    x = precip_predict,
                                    colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_line(data = precip_pred_rqs,
            aes(y = precip_0.5),
            colour = "black",
            linewidth = 1.5) +

  scale_fill_nightfall(name = "Latitude",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(y = "Assemblage Diet Variability",
       x = "Precipitation Predictability",
       linetype = "Latitudinal\nBand") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 22),
        axis.text = element_text(colour = "black",size = 22),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black")); precip_pred_summed_diff_gg



#Temperature
temp_pred_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(temp_predict, 3),
                             data = diff_season_summ_df %>% filter(!is.na(temp_predict)),
                             b = 0.76,
                             p = 0.5)

#Saving model output
temp_pred_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/temp_pred_qr_coefficients.csv")

#Pulling fitted values for all of the points
temp_pred_rqs <- bind_cols(temp_0.5 = temp_pred_qr$fitted.values,
                           diff_season_summ_df %>% 
                             filter(!is.na(temp_predict)))

temp_summed_predict_diff_gg <- ggplot(data = diff_season_summ_df,
                                      aes(y = summed_wgt_diff,
                                          x =  temp_predict,
                                          colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +

  geom_line(data = temp_pred_rqs,
            aes(y = temp_0.5,
                linetype = lat_30_band),
            colour = "black",
            linewidth = 1.5) +
  
  scale_fill_nightfall(name = "Latitude",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_linetype_manual(values = c(1,1,1,1,1,1)) +
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(linetype = "Latitudinal\nBand",
       x = "Temperature Predictability") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 22),
        axis.text = element_text(colour = "black",size = 22),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()); temp_summed_predict_diff_gg

#GPP
gpp_pred_qr <- lqr::Log.lqr(summed_wgt_diff ~ poly(gpp_predict, 3),
                                 data = diff_season_summ_df %>% filter(!is.na(gpp_predict)) %>% droplevels(),
                            p = 0.5,
                            b = 0.76)

#Saving model output
gpp_pred_qr$table[,-c(5)] %>% 
  as.data.frame() %>% 
  mutate_all(~round(., digits = 3)) %>% 
  write.csv("./output/seasonality_regression/gpp_pred_qr_coefficients.csv")

#Pulling fitted values for all of the points
gpp_pred_rqs <- bind_cols(gpp_0.5 = gpp_pred_qr$fitted.values,
                          diff_season_summ_df %>% 
                            filter(!is.na(gpp_predict)))

gpp_summed_predict_diff_gg <- ggplot(data = diff_season_summ_df,
                                      aes(y = summed_wgt_diff,
                                          x =  gpp_predict,
                                          colour = latitude)) + 
  
  geom_point(aes(fill = latitude),
             size = 4,
             alpha = 0.5,
             colour = "black",
             shape = 21) +
  
  geom_line(data = gpp_pred_rqs,
            aes(y = gpp_0.5),
            colour = "black",
            linewidth = 1.5) +
  
  scale_fill_nightfall(name = "Latitude",
                       midpoint = 0,
                       discrete = FALSE) + 
  
  scale_color_nightfall(name = "Latitude",
                        midpoint = 0,
                        reverse = TRUE,
                        discrete = TRUE) + 
  
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                     limits = c(0, 1)) +
  
  labs(x = "GPP Predictability") + 
  
  theme_classic() +
  
  theme(axis.title = element_text(colour = "black",size = 22),
        axis.text = element_text(colour = "black",size = 22),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.key.width = unit(2,"cm"),
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "grey", colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()); gpp_summed_predict_diff_gg



#Creating patchwork plot for combining both graphs
triband_hemi_pred_summed_gg <- precip_pred_summed_diff_gg +
  temp_summed_predict_diff_gg +
  gpp_summed_predict_diff_gg + 
  plot_layout(guides = "collect")

triband_hemi_pred_summed_gg

ggsave("./graphs/seasonality_regression/three_pane_predictability_summed_diff.png",
       width = 22.5, height = 6)
