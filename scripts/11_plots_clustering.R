#Author: Marta Jarzyna

# Libraries ---------------------------------------------------------------

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
require(tidyterra)
require(scales)


#################Code part 1: mapping clusters
#Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
pca_geography_matrix <- readRDS("./output/pca_geography_matrix.rds")
ts_pca <- readRDS("./output/month_space_pca/ts_pca.rds")

# reorganize ms_pca_matrix so that it has nrows = nsites, ncolumns = index*pc
# in columns: index1 x pc1, index2 x pc1..., index1 x pc2, index2 x pc2, ...
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
rownames(ms_pca_matrix_wide) <- c(1:80000)
ms_pca_matrix_wide_nona <- na.omit(ms_pca_matrix_wide) # This is the data that's passed for clustering
grid.names.all <- rownames(ms_pca_matrix_wide)
grid.names.sel <- rownames(ms_pca_matrix_wide_nona)

# read in clustering output
clusters <- readRDS(file=paste0("./output/clustering/kmeans-clusters.rds"))
rownames(clusters) <- grid.names.sel
cluster_goodness <- readRDS(file=paste0("./output/clustering/kmeans-clusters-goodness.rds"))

cluster_goodness <- as.data.frame(cluster_goodness)

p <- ggplot(cluster_goodness, aes(k,sil_width)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Number of Clusters", breaks=seq(0,15,2)) %>%
  + scale_y_continuous(name = "Avg. Silhouette Width") %>%
  + theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "white", size=1),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=22),
    axis.text.y = element_text(colour='black',size=22),
    #axis.title.x = element_text(colour='black',size=22),
    #axis.title.y = element_text(colour='black',size=22),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(color = 'black', size=1.5),
    axis.ticks.length=unit(0.3,"cm"),
    legend.position="none",
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))
p

ggsave(p, file=paste0("./graphics/kmeans-clusters-goodness.png"), bg="transparent", width=8, height=6, dpi=600)

#compute the number of cells, and their proportion, for each cluster
# to be reported in the manuscript
clusters8 <- as.data.frame(clusters[,c(1:12,19)])
#cluster 1
sum(clusters8$k_8 == 1, na.rm = TRUE)
(sum(clusters8$k_8 == 1, na.rm = TRUE)/nrow(clusters8))*100

#cluster 2
sum(clusters8$k_8 == 2, na.rm = TRUE)
(sum(clusters8$k_8 == 2, na.rm = TRUE)/nrow(clusters8))*100

#cluster 3
sum(clusters8$k_8 == 3, na.rm = TRUE)
(sum(clusters8$k_8 == 3, na.rm = TRUE)/nrow(clusters8))*100

#cluster 4
sum(clusters8$k_8 == 4, na.rm = TRUE)
(sum(clusters8$k_8 == 4, na.rm = TRUE)/nrow(clusters8))*100

#cluster 5
sum(clusters8$k_8 == 5, na.rm = TRUE)
(sum(clusters8$k_8 == 5, na.rm = TRUE)/nrow(clusters8))*100

#cluster 6
sum(clusters8$k_8 == 6, na.rm = TRUE)
(sum(clusters8$k_8 == 6, na.rm = TRUE)/nrow(clusters8))*100

#cluster 7
sum(clusters8$k_8 == 7, na.rm = TRUE)
(sum(clusters8$k_8 == 7, na.rm = TRUE)/nrow(clusters8))*100

#cluster 8
sum(clusters8$k_8 == 8, na.rm = TRUE)
(sum(clusters8$k_8 == 8, na.rm = TRUE)/nrow(clusters8))*100


# Fill ms_pca_matrix_complete with values from clusters where rows match
ms_pca_matrix_complete <- matrix(NA, nrow = length(grid.names.all), ncol = ncol(clusters),
                                 dimnames = list(grid.names.all, colnames(clusters)))
ms_pca_matrix_complete[rownames(clusters), ] <- clusters

#plot
model_rast <- rast("./output/model_raster.tif")
clus_rast <- model_rast
values(clus_rast) <- ms_pca_matrix_complete[,19] #local maximum at k=8
clus_rast <- as.factor(clus_rast)


### Prep full grid for plotting
grid_full <- list(
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-180, -120, -60, 0, 60, 120, 180),
                                  function(x) cbind(x, seq(-90, 90, 1)))), crs = 'WGS84')),
  st_sf(geometry = st_sfc(
    st_multilinestring(x = lapply(c(-90, -60, -30, 0, 30, 60, 90), function(x) {
      cbind(seq(-180, 180, 1), x)
    })), crs = 'WGS84'))) %>% bind_rows()

ylabs <- lapply(c(-90, -60, -30, 0, 30, 60, 90 ), function(x) {
  st_sf(label = paste0(abs(x), '\u00b0',
                       ifelse(x == 0, '', ifelse(x < 0, 'S', 'N'))),
        geometry = st_sfc(st_point(c(-180, x)), crs = 'WGS84'))
}) %>% bind_rows()

# Get the continent and ocean shapefile
globe_shp <- ne_download(scale = "medium", type = "land",
                         category = "physical", returnclass = "sf")
globe_shp <- fortify(globe_shp)
ocean_shp <- ne_download(scale = "medium", type = "ocean",
                         category = "physical", returnclass = "sf")
ocean_shp <- fortify(ocean_shp)


custom_colors1 <- c("1" = "midnightblue", "2" = "deeppink3", "8" = "aquamarine2", "4" = "chartreuse1", 
                    "7" = "yellow", "3" = "darkgoldenrod1", "5"="darkorchid3", "6"="brown4")

p <- ggplot() + 
  geom_spatraster(data = clus_rast,
                  maxcell = 80000) +
  scale_fill_manual(values = custom_colors1, na.value = "transparent") + 
  #scale_fill_viridis_d(option = "turbo", na.value = "transparent", direction = 1) +
  #scale_fill_viridis_d(option = "viridis", na.value = "transparent",direction = 1) +
  #geom_sf(data = grid_full_poly, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  #geom_sf(data = grid_full, linewidth = 0.1, colour = "white", fill="grey70") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  coord_sf(crs = "ESRI:54030") 
p
ggsave(p, file=paste0("./graphics/clusters-8_legend.png"), bg = "transparent", width = 14.5, height = 16, dpi=800)

p <- ggplot() + 
  geom_spatraster(data = clus_rast,
                  maxcell = 80000) +
  scale_fill_manual(values = custom_colors1, na.value = "transparent") + 
  #scale_fill_viridis_d(option = "turbo", na.value = "transparent", direction = 1) +
  #scale_fill_viridis_d(option = "viridis", na.value = "transparent",direction = 1) +
  #geom_sf(data = grid_full_poly, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  #geom_sf(data = grid_full, linewidth = 0.1, colour = "white", fill="grey70") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),     
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        legend.position="none",  
        strip.text = element_blank()) +
  coord_sf(crs = "ESRI:54030") 
p
ggsave(p, file=paste0("./graphics/clusters-8.png"), bg = "transparent", width = 14.5, height = 16, dpi=800)



###############Code part 2: plotting quantiles for each cluster
### Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
pca_geography_matrix <- readRDS("./output/pca_geography_matrix.rds")
ts_pca <- readRDS("./output/month_space_pca/ts_pca.rds")

# reorganize ms_pca_matrix so that it has nrows = nsites, ncolumns = index*pc
# in columns: index1 x pc1, index2 x pc1..., index1 x pc2, index2 x pc2, ...
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
rownames(ms_pca_matrix_wide) <- c(1:80000)
ms_pca_matrix_wide_nona <- na.omit(ms_pca_matrix_wide) # This is the data that's passed for clustering
grid.names.all <- rownames(ms_pca_matrix_wide)
grid.names.sel <- rownames(ms_pca_matrix_wide_nona)

# read in clustering output
clusters <- readRDS(file=paste0("./output/clustering/kmeans-clusters.rds"))
rownames(clusters) <- grid.names.sel
cluster_goodness <- readRDS(file=paste0("./output/clustering/kmeans-clusters-goodness.rds")) #8 clusters is best

# Fill ms_pca_matrix_complete with values from clusters where rows match
ms_pca_matrix_complete <- matrix(NA, nrow = length(grid.names.all), ncol = ncol(clusters),
                                 dimnames = list(grid.names.all, colnames(clusters)))
ms_pca_matrix_complete[rownames(clusters), ] <- clusters


### Creating list of grid cells that fall within each cluster region & quantifying quantiles for those clusters
ms_pca_matrix_complete <- as_tibble(ms_pca_matrix_complete)
ms_pca_matrix_complete <- ms_pca_matrix_complete %>%
  mutate(ID_grid = rownames(ms_pca_matrix_complete)) %>% 
  dplyr::select(ID_grid, Tr1_PC1, Tr1_PC2, Tr2_PC1, Tr2_PC2, 
                Tr3_PC1, Tr3_PC2, Tr4_PC1, Tr4_PC2, 
                Tr5_PC1, Tr5_PC2, Tr6_PC1, Tr6_PC2, k_8)
ms_pca_matrix_nona <- na.omit(ms_pca_matrix_complete) 

ms_pca_matrix_melt <- reshape2::melt(ms_pca_matrix_nona, id = c("ID_grid","k_8"),
                                     variable.name = "Variable", value.name = "Value")

vars <- unique(ms_pca_matrix_melt$Variable)

quants.all <- list()
for (j in 1:8){
  var.j <- ms_pca_matrix_melt[ms_pca_matrix_melt[,2] == j,]
  quants.i <- list()
  for (i in 1:length(vars)){
    var.i <- var.j[var.j[,3] == vars[i],]
    quants <- quantile(var.i[,4],probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = FALSE)
    mat <- matrix(rep(quants, nrow(var.i)), nrow = nrow(var.i), byrow = TRUE)
    var.i <- cbind(var.i,mat)
    quants.i <- quants.i %>%
      bind_rows(var.i)
  }
  quants.all <- quants.all %>%
    bind_rows(quants.i)
}

colnames(quants.all) <- c("ID_grid","cluster","var_pc","loading","quant10","quant25","quant50","quant75","quant90")

#assign grouping--if (25,75) <0 group 1, (25,75)>0 group 2, (25,75) overlaps 0 group 3
# this takes a long time - read in a ready-to-go file instead
quants.all$group <- NA
for (i in 1:nrow(quants.all)){
  quants.i <- quants.all[i,]
  if (quants.all[i,6] <0 & quants.all[i,8] <0) {
    quants.all[i,10] <- 1
  } else if (quants.all[i,6] >0 & quants.all[i,8] >0) {
    quants.all[i,10] <- 2
  } else {
    quants.all[i,10] <- 3
  }
}

saveRDS(quants.all, file=paste0("./output/clustering/loadings-quants-clusters.rds"))
quants.all <- readRDS(file=paste0("./output/clustering/loadings-quants-clusters.rds"))
quants.all.lev <- quants.all
quants.all.lev$var_pc <- factor(quants.all.lev$var_pc, 
                                levels = c("Tr1_PC1","Tr2_PC1","Tr3_PC1","Tr4_PC1","Tr5_PC1","Tr6_PC1",
                                           "Tr1_PC2","Tr2_PC2","Tr3_PC2","Tr4_PC2","Tr5_PC2","Tr6_PC2"))

clust.list <- c("clus1","clus2","clus3","clus4","clus5","clus6","clus7","clus8")

# plot for each cluster separately
for (i in 1:8){
  var.i <- quants.all.lev[quants.all.lev[,2] == i,]
  
  p <- ggplot(var.i, aes(x=loading, y=factor(var_pc,levels = rev(levels(var_pc))), fill= as.factor(group))) + 
    geom_boxplot() +
    #geom_violin() +
    #coord_flip() +
    #scale_x_continuous(limits = c(-0.015, 0.015)) +  # Set x-axis limits
    scale_fill_manual(values=c('deepskyblue2','indianred1','grey80')) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size=1.5) + 
    theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  p
  
  ggsave(p, file=paste0("graphics/quants_loading_cluster_",i,".png"), width=8, height=6, dpi=600)
}


### ultimately, we need a custom x axis limits for each cluster
mat_limit <- matrix(NA,8,2)
mat_limit[1,] <- c(-0.002, 0.002)
mat_limit[2,] <- c(-0.015, 0.025)
mat_limit[3,] <- c(-0.01, 0.01)
mat_limit[4,] <- c(-0.015, 0.015)
mat_limit[5,] <- c(-0.01, 0.015)
mat_limit[6,] <- c(-0.02, 0.025)
mat_limit[7,] <- c(-0.02, 0.02)
mat_limit[8,] <- c(-0.0075, 0.01)

for (i in 1:1){ #separate for cluster 1 because there are only two groups: 1 and 3
  var.i <- quants.all.lev[quants.all.lev[,2] == i,]
  
  p <- ggplot(var.i, aes(x=loading, y=factor(var_pc,levels = rev(levels(var_pc))), fill= as.factor(group))) + 
    geom_boxplot(outlier.color = "grey") +
    #geom_violin() +
    #coord_flip() +
    scale_x_continuous(limits =  c(mat_limit[i,1],mat_limit[i,2]), oob=squish) +  # Set x-axis limits
    scale_fill_manual(values=c('deepskyblue2','grey80')) +
    geom_vline(xintercept = 0, linetype = "twodash", color = "darkmagenta", size=1.5) + 
    theme(#panel.border = element_rect(fill=NA, colour = "white", size=1),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=30),
      #axis.text.y = element_text(colour='black',size=22),
      axis.text.y = element_blank(),
      #axis.title.x = element_text(colour='black',size=22),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(colour='black',size=22),
      axis.title.y = element_blank(),
      axis.ticks = element_line(color = 'black', size=1.5),
      axis.ticks.length=unit(0.3,"cm"),
      legend.position="none",
      legend.text=element_text(size=20),
      legend.title=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  #annotate("rect", xmin = mat_limit[i,1], xmax = mat_limit[i,2], ymin = 6.5, ymax = 12.5, alpha = 0.2, fill = "grey80")
  
  ggsave(p, file=paste0("graphics/quants_loading_cluster_",i,"_custom.png"), bg = "transparent", width=8, height=10, dpi=600)
}

for (i in 2:8){
  var.i <- quants.all.lev[quants.all.lev[,2] == i,]
  
  p <- ggplot(var.i, aes(x=loading, y=factor(var_pc,levels = rev(levels(var_pc))), fill= as.factor(group))) + 
    geom_boxplot(outlier.color = "grey") +
    #geom_violin() +
    #coord_flip() +
    scale_x_continuous(limits =  c(mat_limit[i,1],mat_limit[i,2]), oob=squish) +  # Set x-axis limits
    scale_fill_manual(values=c('deepskyblue2','indianred1','grey80')) +
    geom_vline(xintercept = 0, linetype = "twodash", color = "darkmagenta", size=1.5) + 
    theme(#panel.border = element_rect(fill=NA, colour = "white", size=1),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=30),
      #axis.text.y = element_text(colour='black',size=25),
      axis.text.y = element_blank(),
      #axis.title.x = element_text(colour='black',size=25),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(colour='black',size=25),
      axis.title.y = element_blank(),
      axis.ticks = element_line(color = 'black', size=1.5),
      axis.ticks.length=unit(0.3,"cm"),
      legend.position="none",
      legend.text=element_text(size=20),
      legend.title=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  #annotate("rect", xmin = mat_limit[i,1], xmax = mat_limit[i,2], ymin = 6.5, ymax = 12.5, alpha = 0.2, fill = "grey80")
  
  ggsave(p, file=paste0("graphics/quants_loading_cluster_",i,"_custom.png"), bg = "transparent", width=8, height=10, dpi=600)
}

