#This script produces individual rasters for each species based on their diet for each month of the year by replacing their presence values with their LC scores from the LRA. It requires the output of all preceeding scripts.

#Author: Andre Bellve

# Libraries ---------------------------------------------------------------
#Data manipulation
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)

#Geospatial manipulation
#require(terra)
library(matrixStats)

#Visualisation
library(ggplot2)
library(ggrepel)

#Time keepers
library(tictoc)
library(beepr)
library(svMisc)

# Diet scores -------------------------------------------------------------

#Diet scores for each species/month created in 01_diet_pca_scores.R
diet_scores_df <- read.csv("./output/species_diet_scores.csv")

# Raster conversion ---------------------------------------------------------

#Values to iterate through
months <- unique(diet_scores_df$month)
spp <- unique(diet_scores_df$species)

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
      diet_axes_scores <- diet_scores_df %>%
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
    
    #Saving the output list as an .rds for later manipulation. This has a row for every cell in our raster with 0 for absence and PC* scores for presences. Each list element is a separate species
    saveRDS(diet_axes_ls,
            paste0("./output/pca_matrices/", 
                   months[m], "_", 
                   iter_axis, ".rds"))
    #Iteration tracker
    cat("\n Completed! \n")
  }
}


# Assemblage summaries ----------------------------------------------------

#Creating vectors of unique PC axes and months to iterate through
pc_axes <- paste0("PC", 1:6)
months <- unique(diet_scores_df$month)

#For loop to iterate through the matrix files and create summaries (mean and SD) for each month. The output is 6 matrices (one for each PC axis), with a row for each cell in a global raster and a column for each month (12)
for(pc in seq_along(pc_axes)){
  
  #Saving the PC axis name to include in the file name
  iter_axis <- pc_axes[pc]
  cat("\n Beginning", iter_axis, "... \n")
  
  for(m in seq_along(months)){
    
    #Iteration tracker
    cat("\n month", m, "... \n")
    
    
    #Creating a vector of all the pca_matrix files to iterate through
    pca_matrix_ls <- list.files("./output/pca_matrices/",
                                pattern = paste0(months[m], 
                                                 "_", 
                                                 pc_axes[pc]),
                                full.names = TRUE) %>% 
      readRDS()
    
    #Collapsing from a list into a single matrix so that we have monthly matrices to sum over for all the species
    pc_mat <- do.call(cbind, pca_matrix_ls)
    
    #Replacing absences (0's) with NA to avoid skewing assemblage summaries with false PC values -should have done this in the last step, instead of putting in 0's! Not a big computational cost though... since been fixed - absence cells are now NA's!
    #pc_mat[pc_mat == 0] <- NA
    
    #Calculating the arithmetic row means (collapsing the species scores to a single row) to give the mean score for a given cell for that month (averaging over the assemblage)
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
  
  #Saving each PC matrix containing the values for all 12 months with the average of the cells
  saveRDS(monthly_mean_mat,
          paste0("./output/summary_matrices/",
                 "mean_",
                 iter_axis,
                 ".rds"))
  
  #Saving each PC matrix containing the values for all 12 months with the SD of the cells
  saveRDS(monthly_sd_mat,
          paste0("./output/summary_matrices/",
                 "sd_",
                 iter_axis,
                 ".rds"))
  
  #Iteration tracker
  cat("\n Completed! \n")
}