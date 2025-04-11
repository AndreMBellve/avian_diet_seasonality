#This scripts takes the rasters produced in the 00_data_preparation.r script and splits them into monthly maps based on their breeding phenology

#Author: Andre Bellve

# Libraries ---------------------------------------------------------------
# install.packages(c("dplyr", "rgdal", "tidyr", 
#                    "stringr", "terra", "tidyterra", 
#                    "ggplot2", "tictoc", "beepr", "rnaturalearth))
# 

#Data manipulation
library(dplyr)
library(tidyr)
library(stringr)

#Raster data handling
library(terra)
library(tidyterra)
library(rnaturalearth)

#Visualisation
library(ggplot2)

#Time keeping
library(tictoc)
library(beepr)

# Data preparation --------------------------------------------------------

#Reading in the spp_range summary produced in the previous script which summarises the species distributions
spp_range_df <- read.csv("./output/species_range_summary.csv") %>% 
  mutate(species = str_replace(species, "_", " "))

#Reading in breeding phenology data - not, this data only covers species which have distinct breeding and/or non-breeding ranges. All other birds will default to using the resident range.
breed_pheno_df <- read.csv("./data/monthly_phenology/bird_season_by_month.csv") %>% 
  #Removing the birds for which we have no data what so ever (only two species: Camptorhynchus labradorius and Moringilla nivalis)
  filter(certainty != "none")

#Data checks look good - no NA's after initial cleaning
for(i in 9:20){
cat(colnames(breed_pheno_df %>% select(any_of(i))), ": ", any(is.na(breed_pheno_df[,i])), "\n")
}

#Converting the entries to be values 1 to 4 instead of named breeding column
#CHECK THAT THESE ARE CORRECTLY DEFINED!
for(i in c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")){

  #Successive ifelse statements to redefine values
    breed_pheno_df[,i] <- ifelse(breed_pheno_df[,i] == "Resident", 1,
                               ifelse(breed_pheno_df[,i] == "Breeding", 2,
                                     ifelse(breed_pheno_df[,i] == "Non-Breeding", 3, 4)))
#The last value (Passage) defaults to 4. 
}

# Reading in range data ---------------------------------------------------

#List of all species names which will be used to index values in the loop. These match the names in the range summary, the SpatRaster object and the species present in the breeding phenology data.
spp_name <- list.files("./output/season_rasters/",
                       full.names = FALSE, 
                       pattern = ".tif") %>% 
  str_remove("[.]tif")

#Note, these data need to be downloaded from the SAviTraits publication (DOI: doi.org/10.1111/geb.13738)
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

#Updating species names in other objects
#Names that are iterated over
spp_name <- taxa_crosswalk$eBird_Clements_v2022

#Species range summary df
spp_range_df <- spp_range_df %>% 
  dplyr::left_join(taxa_crosswalk, 
                   by = join_by(species == BirdLife_v7)) %>% 
  dplyr::mutate(species = eBird_Clements_v2022) %>% 
  filter(!is.na(species)) %>% 
  dplyr::select(-c(eBird_Clements_v2022))

#Breeding Phenology df
breed_pheno_df <- breed_pheno_df %>% 
  dplyr::left_join(taxa_crosswalk, 
                   by = join_by(species_rast == BirdLife_v7)) %>% 
  mutate(species_rast = eBird_Clements_v2022) %>% 
  filter(!is.na(species_rast)) %>% 
  dplyr::select(-c(eBird_Clements_v2022))
  
#Converting rasters to a large matrix for speed of processing - this will be iterated through to create a new matrix for each month based on the species distribution
range_mat <- as.matrix(range_rast, 
                       wide = FALSE,
                       na.rm = FALSE)

#Empty list to store the loop output in with each element being a matrix of grid cells × species
monthly_mat_ls <- as.list(rep(NA, 12))
#Naming list elements
names(monthly_mat_ls) <- c("Jan", "Feb", "Mar", "Apr", 
                           "May", "Jun", "Jul", "Aug", 
                           "Sep", "Oct", "Nov", "Dec")
#Adding lists within each months list to store the species data
for(i in 1:12){
  monthly_mat_ls[[i]] <- list()
}

# List elements -----------------------------------------------------------
#Code has been adapted to create list elements, rather than building up a matrix as it is much more efficient to run and digest the information. The output will be a list for each month that subsequently is filled with a list of values (80000 × 1) for each species containing the information for the species distribution in each month
tic("Range Phenology Processing")
#For loop to process each raster into a matrix of grid cell × species × month
for(sp in seq_along(spp_name)){
  
  #Pulling out the species in question
  sp_name <- spp_name[sp]
  
  #Tracking progress and problem species
  cat(paste0(sp_name," - ", sp, "/", length(spp_name), " \n"))
  
  #RESIDENT ONLY
  #For species which only have a resident range, we cut the processing down and use the same distribution for all months
  if(spp_range_df[spp_range_df$species == sp_name,]$resident_only){
    
    #Iterating through the months
    for(mon in names(monthly_mat_ls)){
      
      monthly_mat_ls[[mon]][[sp]] <- range_mat[,sp]
      
    }
    
    #ALL BIRDS WITH ANY NON-RESIDENT RANGE
  }else{
    
    #Checking to see if the species has a distinct breeding phenology
    if(spp_name[sp] %in% breed_pheno_df$species_rast){
      
      #Pulling out the monthly breakdown of the species breeding phenology
      spp_seasons <- breed_pheno_df %>% 
        filter(species_rast == spp_name[sp]) %>% 
        select(c("Jan", "Feb", "Mar", "Apr", 
                 "May", "Jun", "Jul", "Aug", 
                 "Sep", "Oct", "Nov", "Dec"))
      
      #Pulling out individual species data
      spp_mat <- range_mat[,sp]
      
      #SPECIES WITH ONLY ONE BREEDING PHENOLOGY
      if(nrow(spp_seasons) == 1){
        
        #Iterating through each month
        for(mon in colnames(spp_seasons)){
          
          #Pulling out the monthly value for the species
          mon_season <- spp_seasons[,mon]
          
          #Checking to see if the phenology value is actually present in the matrix, and if not, recoding as appropriate
          if(!mon_season %in% spp_mat){
            #When there is no passage range information, these values are instead allotted to the non-breeding range
            if(mon_season == 4){
              mon_season <- 3 #Non-Breeding code
            }
          }
          
          #Reclassification table for spp × month combination - 1 indicates presence, while 0 indicates absence
          monthly_mat_ls[[mon]][[sp]] <- ifelse(spp_mat %in% c(1, mon_season), 
                                                1, 0) 
          
        }
        #SPECIES WWITH MULTIPLE BREEDING PHENOLOGIES
      }else{
        
        #Split the raster into northern and southern hemisphere portions
        sp_rast_n <- range_rast[[spp_name[sp]]][1:100, 1:400, drop = FALSE]
        sp_rast_s <- range_rast[[spp_name[sp]]][101:200, 1:400, drop = FALSE]
        
        #Iterating through each month
        for(mon in colnames(spp_seasons)){
          
          #NORTHERN DISTRIBUTION
          #Pulling out the monthly value for the species - This assumes North is always entered first.
          mon_season <- spp_seasons[1, mon]
          
          #Checking to see if the phenology value is present in the matrix, and if not, recoding
          if(!mon_season %in% spp_mat){
            #Where there is no passage range information, these values are instead allotted to the non-breeding range
            if(mon_season == 4){
              mon_season <- 3 #Non-Breeding code
            }
          }
          
          #Reclassification table for spp × month combination - 1 indicates presence, while 0 indicates absence
          n_spp_mon_rast <- ifel(sp_rast_n %in% c(1, mon_season), 
                                 1, 0) 
          
          #SOUTHERN DISTRIBUTION
          #Pulling out the monthly value for the species
          mon_season <- spp_seasons[2, mon]
          
          #Checking to see if the phenology value is present in the matrix, and if not, recoding
          if(!mon_season %in% spp_mat){
            #Where there is no passage range information, these values are instead allotted to the non-breeding range
            if(mon_season == 4){
              mon_season <- 3 #Non-Breeding code
            }
          }
          
          #Reclassification table for spp × month combination - 1 indicates presence, while 0 indicates absence
          s_spp_mon_rast <- ifel(sp_rast_s %in% c(1, mon_season), 1, 0)  
          
          #Binding the two hemispheres back together
          monthly_mat_ls[[mon]][[sp]] <- merge(n_spp_mon_rast, s_spp_mon_rast) %>% 
            as.matrix()
          
        } #End of monthly for loop
      }
    } #End of distinct breeding phenology condition
  } #End of resident condition
  cat("Finished \n")
} toc(); beep() #End of entire loop

#Code to add species names to list elements - (format corrected) and converting them to monthy matrices
for(i in seq_along(monthly_mat_ls)){
  
  #Name for each month
  month <- names(monthly_mat_ls)[i]
  
  #Naming columns in each months species lists
  names(monthly_mat_ls[[i]]) <- paste(str_replace(spp_name, " ", "_"), month, sep = "_")
  
  #Converting the monthly lists for each species into a single matrix with species represented as individual columns and rows corresponding to grid cells
  monthly_mat_ls[[i]] <- simplify2array(monthly_mat_ls[[i]])
  
  #Saving output for each month
  saveRDS(monthly_mat_ls[[i]],
          paste0("./output/monthly_ranges/", 
                 month, ".rds"))
}


