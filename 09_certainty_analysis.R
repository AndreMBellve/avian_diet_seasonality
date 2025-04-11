#Script to determine the effect of uncertain diets on our analysis. It follows an identical workflow to the 02_diet_pca_scores script, but removes rows based on our certainty about a given species diet

#Author: André M. Bellvé

# Libraries ---------------------------------------------------------------

#Data manipulation
require(dplyr)
require(tidyr)
require(stringr)
require(janitor)

#Geospatial manipulation
#require(terra)
require(matrixStats)

#Analysis packages
require(compositions)

#Visualisation
require(ggplot2)
require(ggrepel)

#Time keepers
library(tictoc)
library(beepr)

# Uncertainty analysis ----------------------------------------------------

## 25% certainty ----------------------------------------------------------

#Reading in all the PCA matrices
pca_25_geog_ls <- list.files("./output/certainty_analysis/monthly_ranges/certainty_25/summary_matrices",
                             pattern = "mean_",
                             full.names = TRUE) %>% 
  as.list() %>% 
  lapply(readRDS) %>% 
  lapply(t) 

#Loading a model raster to overwrite values
model_rast <- rep(rast("./output/model_raster.tif"), 12)

#Iteratively looping through the matrices in the list, transforming them to rasters to mask out the lakes and oceans, putting them back into matrixes and then overwriting the values in the original list.
for(i in seq_along(pca_25_geog_ls)){
  
  values(model_rast) <- t(pca_25_geog_ls[[i]])
  
  mask_rast <- mask(model_rast, 
                    ocean_vect,
                    updatevalue = NA,
                    inverse = TRUE, 
                    touches = FALSE) %>% 
    mask(.,
         lake_vect,
         updatevalue = NA,
         inverse = TRUE, 
         touches = FALSE)
  
  pca_25_geog_ls[[i]] <- t(values(mask_rast))
  
}

#Binding them all together to create a matrix of months × PCA/geographic cell
pca_25_geog_mat <- do.call(cbind, pca_25_geog_ls) %>% 
  as.data.frame()

saveRDS(pca_25_geog_mat, 
        "./output/certainty_analysis/pca_25_geography_matrix.rds")

#Removing NA's
na_col_count <- apply(is.na(pca_25_geog_mat),2,sum)

pca_geog_25_na_free_mat <- pca_25_geog_mat[, na_col_count == 0]

ts_25_pca <- prcomp(pca_geog_25_na_free_mat,
                    retx = TRUE, 
                    scale = FALSE)

#Percent variation explained
round(((ts_25_pca$sdev^2)/sum(ts_25_pca$sdev^2)), 
      digits = 2)

saveRDS(ts_25_pca,
        "./output/certainty_analysis/month_space_pca/ts_25_pca.rds")

ts_25_geog_mat <- pca_geog_mat[1:3,]

ts_25_geog_mat[, na_col_count == 0] <- t(ts_25_pca$rotation[,1:3])

row.names(ts_25_geog_mat) <- c("PC1", "PC2", "PC3")

saveRDS(ts_25_geog_mat,
        "./output/certainty_analysis/month_space_pca/ms_pca_25_matrix.rds")


## 50% certainty ----------------------------------------------------------

#Reading in all the PCA matrices
pca_50_geog_ls <- list.files("./output/certainty_analysis/monthly_ranges/certainty_50/summary_matrices",
                             pattern = "mean_",
                             full.names = TRUE) %>% 
  as.list() %>% 
  lapply(readRDS) %>% 
  lapply(t) 

#Loading a model raster to overwrite values
model_rast <- rep(rast("./output/model_raster.tif"), 12)

#Iteratively looping through the matrices in the list, transforming them to rasters to mask out the lakes and oceans, putting them back into matrixes and then overwriting the values in the original list.
for(i in seq_along(pca_50_geog_ls)){
  
  values(model_rast) <- t(pca_50_geog_ls[[i]])
  
  mask_rast <- mask(model_rast, 
                    ocean_vect,
                    updatevalue = NA,
                    inverse = TRUE, 
                    touches = FALSE) %>% 
    mask(.,
         lake_vect,
         updatevalue = NA,
         inverse = TRUE, 
         touches = FALSE)
  
  pca_50_geog_ls[[i]] <- t(values(mask_rast))
  
}

#Binding them all together to create a matrix of months × PCA/geographic cell
pca_50_geog_mat <- do.call(cbind, pca_50_geog_ls) %>% 
  as.data.frame()

saveRDS(pca_50_geog_mat, 
        "./output/certainty_analysis/pca_50_geography_matrix.rds")

#Removing NA's
na_col_count <- apply(is.na(pca_50_geog_mat),2,sum)

pca_geog_50_na_free_mat <- pca_50_geog_mat[, na_col_count == 0]

ts_50_pca <- prcomp(pca_geog_50_na_free_mat,
                    retx = TRUE, 
                    scale = FALSE)

#Percent variation explained
round(((ts_50_pca$sdev^2)/sum(ts_50_pca$sdev^2)), 
      digits = 2)

#biplot(ts_50_pca)

saveRDS(ts_50_pca,
        "./output/certainty_analysis/month_space_pca/ts_50_pca.rds")

ts_50_geog_mat <- pca_geog_mat[1:3,]

ts_50_geog_mat[, na_col_count == 0] <- t(ts_50_pca$rotation[,1:3])

row.names(ts_50_geog_mat) <- c("PC1", "PC2", "PC3")

saveRDS(ts_50_geog_mat,
        "./output/certainty_analysis/month_space_pca/ms_pca_50_matrix.rds")



# Data preparation --------------------------------------------------------

#SAvi Traits diet certainty data
certainty_df <- read.csv("./data/diet_certainty/savi_traits_diet_certainty.csv") %>% 
  dplyr::select(Species_Scientific_Name,
                Certainty)

#SAvi Traits diet database
raw_trait_df <- read.csv('./database-files-v1.0/SAviTraits_1-0_1.csv') %>%
  #Trimming out unnecessary variables
  dplyr::select(Species_Scientific_Name, Diet_Sub_Cat,
                Jan, Feb, Mar, Apr, May, Jun,
                Jul, Aug, Sep, Oct, Nov, Dec) %>%
  left_join(certainty_df)


# #SAvi Traits diet database
# raw_trait_df <- read.csv('./data/diet_data/SAviTraits_1-0_1.csv') %>% 
#   #Trimming out unnecessary variables
#   dplyr::select(Species_Scientific_Name, Diet_Sub_Cat, 
#                 Jan, Feb, Mar, Apr, May, Jun, 
#                 Jul, Aug, Sep, Oct, Nov, Dec) %>% 
#   left_join(certainty_df)

#Reading in the necessary full data to match up the scores for the species so they align
full_trait_df <- raw_trait_df %>% 
  #Preparing the data for distance matrix and manipulation by reformatting to make diets columns and months rows
  pivot_longer(cols=c(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec),
               names_to = "Month",
               values_to = "Proportion") %>% 
  pivot_wider(names_from = "Diet_Sub_Cat",
              values_from = "Proportion") %>% 
  
  #Collapsing the vertebrate diet estimates
  mutate(Vertebrate = Ectotherms + Endotherms + Fish + Unknown) %>% 
  dplyr::select(-c(Ectotherms, Endotherms, Fish, Unknown)) %>% 
  
  #Changing species names to match the naming in our range files
  mutate(Species = str_replace(Species_Scientific_Name, " ", "_"),
         id_uniq = paste(Species, Month, sep = "_")) %>% 
  
  #Standardising naming format
  clean_names() 


# 25% certainty -----------------------------------------------------------
#Collapses the vertebrate diet groups into a single vertebrate diet class to reduce the number of traits analysed and because most of the diet proportions lie in the "unknown" vertebrate category anyway
trait_cer_25_df <- raw_trait_df %>%
  
  #Filtering down to 50% certainty or greater and then removing the column
  filter(Certainty >= 0.25) %>%
  dplyr::select(-Certainty) %>%
  
  #Preparing the data for distance matrix and manipulation by reformatting to make diets columns and months rows
  pivot_longer(cols=c(Jan, Feb, Mar, Apr, May, Jun,
                      Jul, Aug, Sep, Oct, Nov, Dec),
               names_to = "Month",
               values_to = "Proportion") %>%
  pivot_wider(names_from = "Diet_Sub_Cat",
              values_from = "Proportion") %>%
  
  #Collapsing the vertebrate diet estimates
  mutate(Vertebrate = Ectotherms + Endotherms + Fish + Unknown) %>%
  dplyr::select(-c(Ectotherms, Endotherms, Fish, Unknown)) %>%
  
  #Changing species names to match the naming in our range files
  mutate(Species = str_replace(Species_Scientific_Name, " ", "_"),
         id_uniq = paste(Species, Month, sep = "_")) %>%
  
  #Standardising naming format
  clean_names()

#Unique diets to create diet space from these
uni_diet_df <- trait_cer_25_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE)

nrow(uni_diet_df)

#Setting seed so that the multiplicative simple replacement is reproducible
set.seed(42)
#Collapsing to unique entries, effectively ignoring species traits - this will allow us to calculate positions in the ordination space, which we can then use to map onto functional space
clr_diet_df <- trait_cer_25_df %>%
  
  #Isolating the compositional columns
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) %>%
  dplyr::select(-c(species_scientific_name, species, id_uniq, month)) %>%
  as.data.frame() %>%
  
  #Multiplicative simple replacement of the zeros in our dataset
  multRepl(label = 0,
           frac = 1,
           z.warning = 0.95) %>%
  clr()

#Checking the correlations between the original and transformed diet values (just invertebrates)
ori_diet_df <- uni_diet_df %>%
  dplyr::select(names(uni_diet_df)[names(uni_diet_df) %in% names(clr_diet_df)]) %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "diet",
               values_to = "untransformed")

corr_test_df <- clr_diet_df %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "diet",
               values_to = "clr") %>%
  bind_cols(ori_diet_df %>%
              dplyr::select(untransformed))
#Generally all the values are monotonically transformed,

#Fitting a PCA to the centered log ratio transformed data
clr_diet_pca <- prcomp(clr_diet_df,
                       center = FALSE)
saveRDS(clr_diet_pca,
        "./output/certainty_analysis/diet_pca_25_cert.rds")

#Saving the PC loadings to retain variable correlations with the PC axes
clr_diet_pca$rotation %>%
  as.data.frame() %>%
  mutate_all(round, digits = 4) %>%
  mutate(variable = row.names(.)) %>%
  write.csv("./output/certainty_analysis/pca_variable_corr_cer_25.csv",
            row.names = FALSE)

#Variance explained by each PC axis
round(((clr_diet_pca$sdev^2)/sum(clr_diet_pca$sdev^2)),
      digits = 2)

## Species PC scores -------------------------------------------------------

#Reading in the original full diet PCA as there were no significant changes between any of the datasets when correcting for uncertainty
clr_diet_pca <- readRDS("./output/diet_pca/clr_diet_pca.RDS")

#Creating a df that has the unique CLR transformed diets and their PC scores to join back onto the species dataframe
diet_scores_df <- full_trait_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) %>% 
  dplyr::select(-c(species_scientific_name, species, id_uniq, month)) %>% 
  as.data.frame() %>% 
  bind_cols(clr_diet_pca$x)

#Preparing the full data frame to have species scores joined on so we can map these back onto the world through space and time
spp_scores_df <- trait_cer_25_df %>%
  left_join(diet_scores_df,
            by = c("invertebrate", "scavenger", "fruit",
                   "nectar", "seed", "other", "vertebrate"))

#Check that no rows have been missed
spp_scores_df %>%
  dplyr::select(starts_with("PC")) %>%
  drop_na() %>%
  nrow() == nrow(spp_scores_df)
#Same number of rows, so no NA's! Good to go.

#Saving the output from this script, which has all the scores for each species diet
write.csv(spp_scores_df,
          "./output/certainty_analysis/species_diet_scores_cer_25.csv",
          row.names = FALSE)

## Raster conversion ---------------------------------------------------------

#Diet scores for each species/month created in 01_diet_pca_scores.R
spp_scores_df <- read.csv("./output/certainty_analysis/species_diet_scores_cer_25.csv")

#Values to iterate through
months <- unique(spp_scores_df$month)
spp <- unique(spp_scores_df$species)

#Iterating through six of the seven axes (7th axis is excluded because it is equally correlated with all diets and accounts for very little variation)
for(pc in 1:6){
  
  #Saving the PC axis name to include in the file name
  iter_axis <- paste0("PC", pc)
  cat("\n Beginning", iter_axis, "... \n")
  
  
  for(m in seq_along(months)){
    
    #Creating empty list to fill with month by month calculations for all species
    diet_axes_ls <- list() 
    
    #Iteration tracker
    cat("\n month", m, "... \n")
    
    #Reading in the monthly matrix data
    month_spp_ls <- readRDS(paste0("./output/monthly_ranges/", 
                                   months[m], 
                                   ".rds"))
    
    for(sp in seq_along(spp)){
      
      #Subsetting down to a single spp × month combination - this is a matrix of global cell values
      sp_mat <- month_spp_ls[[spp[sp]]]
      
      #Pulling the axis scores for a given species × month to overwrite those in the matrix
      diet_axes_scores <- spp_scores_df %>%
        filter(month == months[m]) %>%
        filter(species == spp[sp])
      
      if(spp[sp] != diet_axes_scores$species){
        stop("Species filtering mismatch - check diet axes scores!")
      }
      
      #Converting the presence values to be the PC scores for that species × month to make them the species scores where they were present
      diet_axes_ls[[sp]] <- ifelse(sp_mat == 1,
                                   diet_axes_scores[,iter_axis],
                                   NA)
    }
    
    matrix_output <- "./output/certainty_analysis/monthly_ranges/certainty_25/pca_matrices/"
    
    if(!dir.exists(matrix_output)){
      dir.create(matrix_output, recursive = TRUE)
    }
    
    #Saving the output list as an .rds for later manipulation. This has a row for every cell in our raster with 0 for absence and PC* scores for presences. Each list element is a separate species
    saveRDS(diet_axes_ls,
            paste0(matrix_output,
                   months[m], "_",
                   iter_axis, ".rds"))
    #Iteration tracker
    cat("\n Completed! \n")
  }
}

## Assemblage summaries ----------------------------------------------------

#Creating vectors of unique PC axes and monthsto iterate through
pc_axes <- paste0("PC", 1:6)
months <- unique(spp_scores_df$month)

#For loop to iterate through the matrix files and create summaries (mean and SD) for each month. The output is 6 matrices (one for each PC axis), with a row for each cell in a global raster and a column for each month (12)
for(pc in seq_along(pc_axes)){
  
  #Saving the PC axis name to include in the file name
  iter_axis <- pc_axes[pc]
  cat("\n Beginning", iter_axis, "... \n")
  
  for(m in seq_along(months)){
    
    #Iteration tracker
    cat("\n month", m, "... \n")
    
    
    #Creating a vector of all the pca_matrix files to iterate through
    pca_matrix_ls <- list.files("./output/certainty_analysis/monthly_ranges/certainty_25/pca_matrices",
                                pattern = paste0(months[m],
                                                 "_",
                                                 pc_axes[pc]),
                                full.names = TRUE) %>%
      readRDS()
    
    #Collapsing from a list into a single matrix so that we have monthly matrices to sum over for all the species
    pc_mat <- do.call(cbind, pca_matrix_ls)
    
    #Replacing absences (0's) with NA to avoid skewing assemblage summaries with false PC values -should have done this in the last step, instead of putting in 0's! Not a big computational cost though... since been fixed - absence cells are now NA's!
    #pc_mat[pc_mat == 0] <- NA
    
    #Calculating the arithematic row means (collapsing the species scores to a single row) to give the mean score for a given cell for that month (averaging over the assemblage)
    pc_mean_mat <- rowMeans(pc_mat,
                            na.rm = TRUE)
    pc_sd_mat <- rowSds(pc_mat,
                        na.rm = TRUE)
    
    #Conditional to either initialise or add to the existing matrix, so that we have a single matrix for each PC containing all the averaged PC scores for each month (one column per month)
    if(m == 1){
      monthly_mean_mat <- pc_mean_mat
      monthly_sd_mat <- pc_sd_mat
    }else{
      monthly_mean_mat <- cbind(monthly_mean_mat,
                                pc_mean_mat)
      monthly_sd_mat <- cbind(monthly_sd_mat,
                              pc_sd_mat)
    }
  }
  
  #Creating output directory in case it doesn't exist
  summary_output <- "./output/certainty_analysis/monthly_ranges/certainty_25/summary_matrices/"
  if(!dir.exists(summary_output)){
    dir.create(summary_output)
  }
  
  #Saving each PC matrix containing the values for all 12 months with the average of the cells
  saveRDS(monthly_mean_mat,
          paste0(summary_output,
                 "mean_",
                 iter_axis,
                 ".rds"))
  
  #Saving each PC matrix containing the values for all 12 months with the SD of the cells
  saveRDS(monthly_sd_mat,
          paste0(summary_output,
                 "sd_",
                 iter_axis,
                 ".rds"))
  
  #Iteration tracker
  cat("\n Completed! \n")
}

# 50% certainty -----------------------------------------------------------
#Collapses the vertebrate diet groups into a single vertebrate diet class to reduce the number of traits analysed and because most of the diet proportions lie in the "unknown" vertebrate category anyway
trait_cer_50_df <- raw_trait_df %>%

  #Filtering down to 50% certainty or greater and then removing the column
  filter(Certainty >= 0.5) %>%
  dplyr::select(-Certainty) %>%

  #Preparing the data for distance matrix and manipulation by reformatting to make diets columns and months rows
  pivot_longer(cols=c(Jan, Feb, Mar, Apr, May, Jun,
                      Jul, Aug, Sep, Oct, Nov, Dec),
               names_to = "Month",
               values_to = "Proportion") %>%
  pivot_wider(names_from = "Diet_Sub_Cat",
              values_from = "Proportion") %>%

  #Collapsing the vertebrate diet estimates
  mutate(Vertebrate = Ectotherms + Endotherms + Fish + Unknown) %>%
  dplyr::select(-c(Ectotherms, Endotherms, Fish, Unknown)) %>%

  #Changing species names to match the naming in our range files
  mutate(Species = str_replace(Species_Scientific_Name, " ", "_"),
         id_uniq = paste(Species, Month, sep = "_")) %>%

  #Standardising naming format
  clean_names()

#Unique diets to create diet space from these
uni_diet_df <- trait_cer_50_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE)

nrow(uni_diet_df)

#Setting seed so that the multiplicative simple replacement is reproducible
set.seed(42)
#Collapsing to unique entries, effectively ignoring species traits - this will allow us to calculate positions in the ordination space, which we can then use to map onto functional space
clr_diet_df <- trait_cer_50_df %>%

  #Isolating the compositional columns
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) %>%
  dplyr::select(-c(species_scientific_name, species, id_uniq, month)) %>%
  as.data.frame() %>%

  #Multiplicative simple replacement of the zeros in our dataset
  multRepl(label = 0,
           frac = 1,
           z.warning = 0.95) %>%
  clr()

#Checking the correlations between the original and transformed diet values (just invertebrates)
ori_diet_df <- uni_diet_df %>%
  dplyr::select(names(uni_diet_df)[names(uni_diet_df) %in% names(clr_diet_df)]) %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "diet",
               values_to = "untransformed")

corr_test_df <- clr_diet_df %>%
  as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.),
               names_to = "diet",
               values_to = "clr") %>%
  bind_cols(ori_diet_df %>%
              dplyr::select(untransformed))
#Generally all the values are monotonically transformed,

#Fitting a PCA to the centered log ratio transformed data
clr_diet_pca <- prcomp(clr_diet_df,
                       center = FALSE)
saveRDS(clr_diet_pca,
        "./output/certainty_analysis/diet_pca_50_cert.rds")

#Saving the PC loadings to retain variable correlations with the PC axes
clr_diet_pca$rotation %>%
  as.data.frame() %>%
  mutate_all(round, digits = 4) %>%
  mutate(variable = row.names(.)) %>%
  write.csv("./output/certainty_analysis/pca_variable_corr_cer_50.csv",
            row.names = FALSE)

#Variance explained by each PC axis
round(((clr_diet_pca$sdev^2)/sum(clr_diet_pca$sdev^2)),
      digits = 2)

## Species PC scores -------------------------------------------------------

#Reading in the original full diet PCA as there were no significant changes between any of the datasets when correcting for uncertainty. Reloading as it was overwritten above by the certainty one...
clr_diet_pca <- readRDS("./output/diet_pca/clr_diet_pca.RDS")

#Creating a df that has the unique CLR transformed diets and their PC scores to join back onto the species dataframe
diet_scores_df <- full_trait_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) %>% 
  dplyr::select(-c(species_scientific_name, species, id_uniq, month)) %>% 
  as.data.frame() %>% 
  bind_cols(clr_diet_pca$x)

#Preparing the full data frame to have species scores joined on so we can map these back onto the world through space and time
spp_scores_df <- trait_cer_50_df %>%
  left_join(diet_scores_df,
            by = c("invertebrate", "scavenger", "fruit",
                   "nectar", "seed", "other", "vertebrate"))

#Check that no rows have been missed
spp_scores_df %>%
  dplyr::select(starts_with("PC")) %>%
  drop_na() %>%
  nrow() == nrow(spp_scores_df)
#Same number of rows, so no NA's! Good to go.

#Saving the output from this script, which has all the scores for each species diet
write.csv(spp_scores_df,
          "./output/certainty_analysis/species_diet_scores_cer_50.csv",
          row.names = FALSE)

## Raster conversion ---------------------------------------------------------

#Diet scores for each species/month created in 01_diet_pca_scores.R
spp_scores_df <- read.csv("./output/certainty_analysis/species_diet_scores_cer_50.csv")

#Values to iterate through
months <- unique(spp_scores_df$month)
spp <- unique(spp_scores_df$species)

#Iterating through six of the seven axes (7th axis is excluded because it is equally correlated with all diets and accounts for very little variation)
for(pc in 1:6){
  
  #Saving the PC axis name to include in the file name
  iter_axis <- paste0("PC", pc)
  cat("\n Beginning", iter_axis, "... \n")
  
  
  for(m in seq_along(months)){
    
    #Creating empty list to fill with month by month calculations for all species
    diet_axes_ls <- list() 
    
    #Iteration tracker
    cat("\n month", m, "... \n")
    
    #Reading in the monthly matrix data
    month_spp_ls <- readRDS(paste0("./output/monthly_ranges/", 
                                   months[m], 
                                   ".rds"))
    
    for(sp in seq_along(spp)){
      
      #Subsetting down to a single spp × month combination - this is a matrix of global cell values
      sp_mat <- month_spp_ls[[spp[sp]]]
      
      #Pulling the axis scores for a given species × month to overwrite those in the matrix
      diet_axes_scores <- spp_scores_df %>%
        filter(month == months[m]) %>%
        filter(species == spp[sp])
      
      if(spp[sp] != diet_axes_scores$species){
        stop("Species filtering mismatch - check diet axes scores!")
      }

      #Converting the presence values to be the PC scores for that species × month to make them the species scores where they were present
      diet_axes_ls[[sp]] <- ifelse(sp_mat == 1,
                                   diet_axes_scores[,iter_axis],
                                   NA)
    }

    matrix_output <- "./output/certainty_analysis/monthly_ranges/certainty_50/pca_matrices/"

    if(!dir.exists(matrix_output)){
    dir.create(matrix_output, recursive = TRUE)
    }

    #Saving the output list as an .rds for later manipulation. This has a row for every cell in our raster with 0 for absence and PC* scores for presences. Each list element is a separate species
    saveRDS(diet_axes_ls,
            paste0(matrix_output,
                   months[m], "_",
                   iter_axis, ".rds"))
    #Iteration tracker
    cat("\n Completed! \n")
  }
}

## Assemblage summaries ----------------------------------------------------

#Creating vectors of unique PC axes and monthsto iterate through
pc_axes <- paste0("PC", 1:6)
months <- unique(spp_scores_df$month)

#For loop to iterate through the matrix files and create summaries (mean and SD) for each month. The output is 6 matrices (one for each PC axis), with a row for each cell in a global raster and a column for each month (12)
for(pc in seq_along(pc_axes)){

  #Saving the PC axis name to include in the file name
  iter_axis <- pc_axes[pc]
  cat("\n Beginning", iter_axis, "... \n")

  for(m in seq_along(months)){

    #Iteration tracker
    cat("\n month", m, "... \n")


    #Creating a vector of all the pca_matrix files to iterate through
    pca_matrix_ls <- list.files("./output/certainty_analysis/monthly_ranges/certainty_50/pca_matrices",
                                pattern = paste0(months[m],
                                                 "_",
                                                 pc_axes[pc]),
                                full.names = TRUE) %>%
      readRDS()

    #Collapsing from a list into a single matrix so that we have monthly matrices to sum over for all the species
    pc_mat <- do.call(cbind, pca_matrix_ls)

    #Replacing absences (0's) with NA to avoid skewing assemblage summaries with false PC values -should have done this in the last step, instead of putting in 0's! Not a big computational cost though... since been fixed - absence cells are now NA's!
    #pc_mat[pc_mat == 0] <- NA

    #Calculating the arithematic row means (collapsing the species scores to a single row) to give the mean score for a given cell for that month (averaging over the assemblage)
    pc_mean_mat <- rowMeans(pc_mat,
                            na.rm = TRUE)
    pc_sd_mat <- rowSds(pc_mat,
                        na.rm = TRUE)

    #Conditional to either initialise or add to the existing matrix, so that we have a single matrix for each PC containing all the averaged PC scores for each month (one column per month)
    if(m == 1){
      monthly_mean_mat <- pc_mean_mat
      monthly_sd_mat <- pc_sd_mat
    }else{
      monthly_mean_mat <- cbind(monthly_mean_mat,
                                pc_mean_mat)
      monthly_sd_mat <- cbind(monthly_sd_mat,
                              pc_sd_mat)
    }
  }

  #Creating output directory in case it doesn't exist
  summary_output <- "./output/certainty_analysis/monthly_ranges/certainty_50/summary_matrices/"
  if(!dir.exists(summary_output)){
      dir.create(summary_output)
  }

  #Saving each PC matrix containing the values for all 12 months with the average of the cells
  saveRDS(monthly_mean_mat,
          paste0(summary_output,
                 "mean_",
                 iter_axis,
                 ".rds"))

  #Saving each PC matrix containing the values for all 12 months with the SD of the cells
  saveRDS(monthly_sd_mat,
          paste0(summary_output,
                 "sd_",
                 iter_axis,
                 ".rds"))

  #Iteration tracker
  cat("\n Completed! \n")
}

