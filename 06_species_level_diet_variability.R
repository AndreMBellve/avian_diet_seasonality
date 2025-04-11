#Script to quantify species-level temporal variability in diet`

#Author: Marta Jarzyna`

# Libraries ---------------------------------------------------------------
rm(list = ls())

#Data manipulation
library(dplyr)
library(terra)

#Geospatial manipulation
library(sf)
library(terra)
library(rnaturalearth)
require(cluster)
require(fastcluster)
require(fpc)
require(ggdensity)
require(tidyr)
require(tidyverse)

# Step 1: quantify diet variability
# Read in the diet values for all species (raw and PC scores)
spp_all_scores <- read.csv("./output/species_diet_scores.csv") #This is a result of pca: clr_diet_pca$x

# Quantify temporal variability as the difference between max and min
variance_weighted_sum <- spp_all_scores %>%
  group_by(species_scientific_name) %>%
  summarise(
    PC1_var = max(PC1, na.rm = TRUE) - min(PC1, na.rm = TRUE),
    PC2_var = max(PC2, na.rm = TRUE) - min(PC2, na.rm = TRUE),
    PC3_var = max(PC3, na.rm = TRUE) - min(PC3, na.rm = TRUE),
    PC4_var = max(PC4, na.rm = TRUE) - min(PC4, na.rm = TRUE),
    PC5_var = max(PC5, na.rm = TRUE) - min(PC5, na.rm = TRUE),
    PC6_var = max(PC6, na.rm = TRUE) - min(PC6, na.rm = TRUE), 
    na.rm = TRUE) 
  
#Min-max feature scaling function
mm_feature_scale <- function(x){
  ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

#Scaling the difference to be between 0 and 1 (all LCs on the same scale)
variance_weighted_sum <- variance_weighted_sum %>%
  mutate(across(ends_with("_var"), mm_feature_scale))

#Reading in the eigenvalues from our diet LRA to weight our scores by
eigenval_df <- readRDS("./output/diet_pca/diet_pca.rds")$sdev %>% 
  data.frame(eigen_prop = round((.^2 / sum(.^2)), digits = 4),
             PC = paste0("PC", 1:7)) %>% 
  dplyr::select(-c(`.`))

#weighting the values
variance_weighted_sum_FIN <- variance_weighted_sum %>%
#Pivoting to long format to assign the eigenvalue proportions
pivot_longer(cols = ends_with("_var"),
             names_to = "PC",
             values_to = "var") %>% 
  
  #Modifying names
  mutate(PC = str_remove(PC, "_var")) %>% 
  
  #PC score summary
  left_join(eigenval_df,
            join_by("PC")) %>% 
  
  #Weighting by variance explained
  #Multiplying all PC's by the eigen proportion to have their weighting reflected
  mutate(eigen_wgt_var = var * eigen_prop,
         PC = toupper(PC)) 

#sum weighted variance
variance_weighted_sum_FIN <- variance_weighted_sum_FIN %>%
  group_by(species_scientific_name) %>%
  summarise(summed_wgt_var = sum(eigen_wgt_var))

# Read in species range designations
spp_range_df <- read.csv("./output/species_range_summary_clements2022.csv") %>%
mutate(species_scientific_name = str_replace(species, "_", " ")) %>%
  mutate(species = str_replace(species_scientific_name, " ", "_"))

# Filter spp_range_df to include only species present in variance_weighted_sum_FIN
spp_range_df_filter <- spp_range_df[spp_range_df$species_scientific_name %in% variance_weighted_sum_FIN$species_scientific_name, ]

# Join the two data frames by the "species" column
merged_df <- merge(spp_range_df_filter, variance_weighted_sum_FIN, by = "species_scientific_name")
sum(merged_df$resident_only == TRUE, na.rm = TRUE)
sum(merged_df$resident_only == FALSE, na.rm = TRUE)


# Step 2: Get species mean latitude of its range
#List of all species names which will be used to index values in the loop. These match the names in the range summary, the SpatRaster object and the species present in the breeding phenology data.
spp_name <- list.files("./output/season_rasters/",
                       full.names = FALSE, 
                       pattern = ".tif") %>% 
  str_remove("[.]tif")

#Reading in taxonomy crosswalk to match the SAviTraits database
taxa_crosswalk <- read.csv("./data/taxonomy_crosswalk/SAviTraits_1-0_3.csv") %>%
  #Just the relevant taxonomies
  dplyr::select(eBird_Clements_v2022, BirdLife_v7) %>% 
  #Dropping all the rows which repeat a particular species
  distinct_all() %>% 
  #Removing species which do not have a range raster
  dplyr::filter(!is.na(BirdLife_v7)) %>% 
  #Arranging them all alphabetically
  dplyr::arrange(BirdLife_v7)

#Reading in a file of global lakes to mask them out of our range analysis
if(!dir.exists("./data/masking_files/")){
  dir.create("./data/masking_files/")
}
lakes_vect <- ne_download(category = "physical",
                          type = "lakes",
                          scale = 110,
                          returnclass = "sv",
                          destdir = "./data/masking_files/")

#Reading in all the raster files of the global distribution - this takes about a minute to run. Species that are in both hemispheres will be separately processed by splitting the raster at the equator, creating individual monthly summaries and binding them back together.
range_rast <- list.files("./output/season_rasters/",
                         full.names = TRUE, 
                         pattern = ".tif") %>%
  rast() %>% 
  
  #Masking out all the major lakes from our maps - updating these values to be zero to exclude them from our analysis without fundamentally altering the nature of the raster data and potentially messing up downstream code
  mask(lakes_vect, 
       inverse = TRUE, 
       updatevalue = 0)

#Altering names to species to check values match the file names
names(range_rast) <- spp_name

#Removing species which are not present in the SAviTraits database
range_rast <- tidyterra::select(range_rast, 
                                all_of(taxa_crosswalk$BirdLife_v7))

#Overwriting names in range_rast with the names that SAviTraits uses
names(range_rast) <- taxa_crosswalk$eBird_Clements_v2022

spp_name <- taxa_crosswalk$eBird_Clements_v2022

sp_lat <- as.data.frame(spp_name)
sp_lat$medianlat <- NA
sp_lat$lowerlat <- NA
sp_lat$upperlat <- NA
sp_lat$difflat <- NA

# Loop through each layer in the raster stack (this takes a long time)
#for (i in 1:nrow(sp_lat)) {
  for (i in 1:10) {
  cat(i)
  # Get the current raster layer
  spp <- range_rast[[i]]
  
  # Extract cell values and coordinates
  cell_data <- as.data.frame(cbind(values(spp), xyFromCell(spp, 1:ncell(spp))))
  colnames(cell_data) <- c("value", "x", "y")
  
  # Filter for cells with values 1, 2, 3, 4
  filtered_cells <- subset(cell_data, value %in% c(1, 2, 3, 4))
  
  # Calculate quantiles of latitude or return NA
  if (nrow(filtered_cells) > 0) {
    sp_lat[i,2] <- quantile(filtered_cells$y, probs = 0.5)
    sp_lat[i,3] <- quantile(filtered_cells$y, probs = 0.025)
    sp_lat[i,4] <- quantile(filtered_cells$y, probs = 0.975)
    #sp_lat[i,2] <- mean(filtered_cells$y)
    #sp_lat[i,3] <- min(filtered_cells$y)
    #sp_lat[i,4] <- max(filtered_cells$y)
    sp_lat[i,5] <- sp_lat[i,4] - sp_lat[i,3]
  } else {
    sp_lat[i,2] <- NA
    sp_lat[i,3] <- NA
    sp_lat[i,4] <- NA
    sp_lat[i,5] <- NA
  }
}

sp_lat_sub <- sp_lat %>%
  mutate(species_scientific_name = spp_name) %>%
  mutate(species = str_replace(species_scientific_name, " ", "_")) %>%
  dplyr::select(species_scientific_name, medianlat, lowerlat, upperlat, difflat)

#merge with the rest of data
merged_df_all <- merge(merged_df, sp_lat_sub, by = "species_scientific_name")

#assign short vs long distance migrant label - a bit arbitrary, but we will use 45 degrees lat as cutoff
merged_df_all$mig_type <- with(merged_df_all, 
                               ifelse(resident_only == TRUE, 0, 
                                      ifelse(resident_only == FALSE & difflat < 45, 1, 
                                             ifelse(resident_only == FALSE & difflat >= 45, 2, NA))))

write.csv(merged_df_all, file = "./output/species_diet_variability.csv", row.names = FALSE)
merged_df_all <- read.csv(file = "./output/species_diet_variability.csv")
