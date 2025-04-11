#This script creates assemblage level summary rasters of for each month × LC across the globe by averaging over all species rasters. Additionally, it lays the ground work for the sensitivity analysis by following the same process but excluding species with based on the certainty in their dietary traits.

#Author: Andre Bellve

# Libraries ---------------------------------------------------------------

#Data manipulation
library(dplyr)
library(terra)
library(tictoc)

#Geospatial manipulation
library(sf)
library(terra)
library(rnaturalearth)

# Data preparation --------------------------------------------------------

#Reading in all the PCA matrices
tic()
pca_geog_ls <- list.files("./output/summary_matrices",
                          pattern = "mean_",
                          full.names = TRUE) %>% 
  as.list() %>% 
  lapply(readRDS) %>% 
  lapply(t) 
toc()

#Loading in masking layers to remove the oceans and lakes from our species rasters
ocean_vect <- ne_download(scale = 110,
                          type = "ocean",
                          category = "physical",
                          returnclass = "sv")

lake_vect <- ne_download(scale = 110,
                         type = "lakes",
                         category = "physical",
                         returnclass = "sv")


#Creating a model raster to ensure that projection and dimensions match
# model_rast <- rast("./output/season_rasters/Abeillia abeillei.tif")
# values(model_rast) <- NaN
# names(model_rast) <- "model_raster"
# writeRaster(model_rast,
#             "./output/model_raster.tif",
#             overwrite = TRUE)

model_rast <- rep(rast("./output/model_raster.tif"), 12)

#Iteratively looping through the matrices in the list, transforming them to rasters to mask out the lakes and oceans, putting them back into matrixes and then overwriting the values in the original list.
for(i in seq_along(pca_geog_ls)){
  
  values(model_rast) <- t(pca_geog_ls[[i]])
    
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
    
  pca_geog_ls[[i]] <- t(values(mask_rast))
    
}

#Binding them all together to create a matrix of months × PCA/geographic cell
pca_geog_mat <- do.call(cbind, pca_geog_ls) %>% 
  as.data.frame()

saveRDS(pca_geog_mat, 
        "./output/pca_geography_matrix.rds")

#Identifying the NaN cells in our raster's value matrix
na_col_count <- apply(is.na(pca_geog_mat), 2, sum)

#Subsetting to the non-NaN values to create a matrix that can be fed into out principle component analysis
pca_geog_na_free_mat <- pca_geog_mat[, na_col_count == 0]

#Fitting a principal component analysis with 12 rows (one for each month) and 80,000 * our six original PCA axes,
ts_pca <- prcomp(pca_geog_na_free_mat,
                 retx = TRUE, 
                 scale = FALSE)

#Checking variance explaind by our axes
round(((ts_pca$sdev^2)/sum(ts_pca$sdev^2)), 
      digits = 2)

#Saving the PCA object in full
saveRDS(ts_pca,
        "./output/month_space_pca/ts_pca.rds")

#Subsetting three rows from our original matrix to initialise a object with the correct number of dimensions - these three rows are placeholders for the first three axes from our time-space PCA
ts_geog_mat <- pca_geog_mat[1:3,]

#Overwriting the placeholders with our PCA scores from the first three axes
ts_geog_mat[, na_col_count == 0] <- t(ts_pca$rotation[,1:3])

#Renaming the rows for clarity
row.names(ts_geog_mat) <- c("PC1", "PC2", "PC3")

#Saving the matrix of tsPCA scores. 
saveRDS(ts_geog_mat,
        "./output/month_space_pca/ms_pca_matrix.rds")

