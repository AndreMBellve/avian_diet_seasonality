#Script to find clusters of congruent patterns in seasonality in dietary characteristics
#Author: Marta A. Jarzyna

# Libraries ---------------------------------------------------------------

#Data manipulation
library(dplyr)
library(terra)

#Geospatial manipulation
library(sf)
library(rnaturalearth)
require(cluster)
require(fastcluster)
require(fpc)

#Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
pca_geography_matrix <- readRDS("./output/pca_geography_matrix.rds")
ts_pca <- readRDS("./output/month_space_pca/ts_pca.rds")

# reorganize ms_pca_matrix so that it has nrows = nsites, ncolumns = dietaxis*pc
# in columns: dietaxis1 x pc1, dietaxis2 x pc1..., dietaxis1 x pc2, dietaxis2 x pc2, ...
ms_pca_matrix <- t(ms_pca_matrix)
row_index_df <- data.frame(begin = c(1, 80001, 160001, 240001, 320001, 400001),
                           end = c(80000, 160000, 240000, 320000, 400000, 480000))

ms_pca_matrix_wide <- cbind(ms_pca_matrix[row_index_df[1,1]:row_index_df[1,2],],
                            ms_pca_matrix[row_index_df[2,1]:row_index_df[2,2],],
                            ms_pca_matrix[row_index_df[3,1]:row_index_df[3,2],],
                            ms_pca_matrix[row_index_df[4,1]:row_index_df[4,2],],
                            ms_pca_matrix[row_index_df[5,1]:row_index_df[5,2],],
                            ms_pca_matrix[row_index_df[6,1]:row_index_df[6,2],])
  
coln <- c("Tr1_PC1", "Tr1_PC2", "Tr2_PC1", "Tr2_PC2", 
          "Tr3_PC1", "Tr3_PC2", "Tr4_PC1", "Tr4_PC2", 
          "Tr5_PC1", "Tr5_PC2", "Tr6_PC1", "Tr6_PC2") 

colnames(ms_pca_matrix_wide) <- coln

##################################################
##  Process Clustering - K-Means
##################################################
cluster.list <- 2:10
  
rownames(ms_pca_matrix_wide) <- c(1:80000)
ms_pca_matrix_wide_nona <- na.omit(ms_pca_matrix_wide) 
grid.names.all <- rownames(ms_pca_matrix_wide)
grid.names.sel <- rownames(ms_pca_matrix_wide_nona)

clusters <- matrix(NA,dim(ms_pca_matrix_wide_nona)[1],length(cluster.list))
colnames(clusters) <- paste("k_",cluster.list, sep="")
  
cluster_goodness <- matrix(NA,length(cluster.list), 5)
rownames(cluster_goodness) <- paste("k_",cluster.list, sep="")
colnames(cluster_goodness) <- c("k", "totss","betweenss","perc_betw","sil_width")
  
##  Calculate the k-means clustering
for (k in cluster.list) {
  cat(k)
    clust_result <- kmeansruns(ms_pca_matrix_wide_nona, krange=k, runs=10,criterion="asw")
    clusters[,k-1]=as.factor(clust_result$cluster)
    cluster_goodness[k-1,] <- c(k, clust_result$totss,clust_result$betweenss,100*(clust_result$betweenss/clust_result$totss),clust_result$crit[k])
}
  
clusters <- cbind(ms_pca_matrix_wide_nona,clusters)
saveRDS(clusters, file=paste0("./output/clustering/kmeans-clusters.rds"))
saveRDS(cluster_goodness, file=paste0("./output/clustering/kmeans-clusters-goodness.rds"))

clusters <- readRDS(file=paste0("./output/clustering/kmeans-clusters.rds"))
cluster_goodness <- readRDS(file=paste0("./output/clustering/kmeans-clusters-goodness.rds"))


