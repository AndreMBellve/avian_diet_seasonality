#Script to plot dietary space at the assemblage level
#Author: Marta Jarzyna

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
require(ggdensity)
require(tibble)


##################################################
##  Plot the entire diet space regardless of the season
##################################################
spp_all_scores <- read.csv("./output/species_diet_scores.csv") #This is a result of pca: clr_diet_pca$x
spp_month_scores <- read.csv("./output/species_month_pc_scores.csv") #This is a result of pca: clr_diet_pca$x
diet_loadings <- read.csv("./output/pca_variable_corr.csv") #this is the output of the pca: clr_diet_pca$rotation (eigenvectors)

## Draft plots to get axes labels
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC1,y = PC2), alpha=0.5, colour="grey80") +
  #Variable axes
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2),
               colour = "red") +
  # Add labels for each variable
  geom_text(data = diet_loadings, aes(x = PC1 * 2, y = PC2 * 2, label = variable),
            color = "blue", hjust = 1.2, vjust = 1.2) 
p

#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC3,y = PC4), alpha=0.5, colour="grey80") +
  #Variable axes
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2),
               colour = "red") +
  # Add labels for each variable
  geom_text(data = diet_loadings, aes(x = PC3 * 2, y = PC4 * 2, label = variable),
            color = "blue", hjust = 1.2, vjust = 1.2)
p

#plot all points for pc 5 and 6
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC5, y = PC6), alpha=0.5, colour="grey80") +
  #Variable axes
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2),
               colour = "red") +
  # Add labels for each variable
  geom_text(data = diet_loadings, aes(x = PC5 * 2, y = PC6 * 2, label = variable),
            color = "blue", hjust = 1.2, vjust = 1.2) 
p


# now add the KDEs, proper vectors for each diet axis
pal <- "#607B8B"
density_probs  <-  c(0.99, 0.9, 0.8, 0.7)

#only plot diet axis that have loading >=0.6
diet_load_pc12_sign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar")
diet_load_pc12_nonsign <- diet_loadings %>%
  filter(variable == "fruit" | variable == "scavenger" | variable == "other" | variable == "vertebrate")

#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  #geom_hdr(data = spp_all_scores, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill=pal)
  geom_hdr(data = spp_all_scores, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = pal) +
  #Adding high density values
  geom_hdr_lines(data = spp_all_scores, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = pal) +
  #Variable axes
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  geom_segment(data = diet_load_pc12_nonsign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "grey60", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  # Add labels for each variable
  #geom_text(data = diet_loadings, aes(x = PC1 * 2, y = PC2 * 2, label = variable),
  #          color = "blue", hjust = 1.2, vjust = 1.2) +
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
diet_load_pc34_sign <- diet_loadings %>%
  filter(variable == "other" | variable == "fruit")
diet_load_pc34_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar" | variable == "scavenger" | variable == "vertebrate")

p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = spp_all_scores, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = pal) +
  #Adding high density values
  geom_hdr_lines(data = spp_all_scores, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = pal) +
  #Variable axes
  geom_segment(data = diet_load_pc34_sign,
               aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  geom_segment(data = diet_load_pc34_nonsign,
               aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2), colour = "grey60", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  # Add labels for each variable
  #geom_text(data = diet_loadings, aes(x = PC1 * 2, y = PC2 * 2, label = variable),
  #          color = "blue", hjust = 1.2, vjust = 1.2) +
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 5 and 6
diet_load_pc56_sign <- diet_loadings %>%
  filter(variable == "scavenger" | variable == "vertebrate")
diet_load_pc56_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "other" | variable == "fruit" | variable == "seed" | variable == "nectar")

p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  geom_point(data = spp_all_scores, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = spp_all_scores, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = pal) +
  #Adding high density values
  geom_hdr_lines(data = spp_all_scores, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = pal) +
  #Variable axes
  geom_segment(data = diet_load_pc56_sign,
               aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  geom_segment(data = diet_load_pc56_nonsign,
               aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2), colour = "grey60", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  # Add labels for each variable
  #geom_text(data = diet_loadings, aes(x = PC1 * 2, y = PC2 * 2, label = variable),
  #          color = "blue", hjust = 1.2, vjust = 1.2) +
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6.png"), bg = "transparent", width=8, height=8, dpi=600)


##################################################
##  Plot the entire diet space for Jun-Aug and Dec-Feb separately
##################################################
winter_pal <- "#180F3EFF"
summer_pal <- "#FFC125"
density_probs  <-  c(0.9, 0.8, 0.7)

#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = dec_feb, aes(x = PC1,y = PC2), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +

  #summer
  geom_point(data = jun_aug, aes(x = PC1,y = PC2), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +

  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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

ggsave(p, file=paste0("./graphics/traitspace_kde_SEASWINSUM_pc1_2.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = dec_feb, aes(x = PC3,y = PC4), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
  #summer
  geom_point(data = jun_aug, aes(x = PC3,y = PC4), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +

  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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


ggsave(p, file=paste0("./graphics/traitspace_kde_SEASWINWUM_pc3_4.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = dec_feb, aes(x = PC5,y = PC6), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +

  #summer
  geom_point(data = jun_aug, aes(x = PC5,y = PC6), alpha=0.8, colour="grey80", size=0.8, position="jitter") +
  geom_hdr(data = jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
  geom_segment(data = diet_loadings,
               aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2), colour = "black", size=2.5,
               arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +

  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
    panel.border = element_rect(fill=NA, colour = "black", size=1.5),
    axis.line = element_line(color = 'black', size=1.5),
    plot.title = element_text(size=15, vjust=2, family="sans"),
    axis.text.x = element_text(colour='black',size=30),
    axis.text.y = element_text(colour='black',size=30),
    #axis.text.y = element_blank(),
    #axis.title.x = element_text(colour='black',size=25),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(colour='black',size=25),
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


ggsave(p, file=paste0("./graphics/traitspace_kde_SEASWINSUM_pc5_6.png"), bg = "transparent", width=8, height=8, dpi=600)


##################################################
##  Plot the  diet space for each cluster and season (Jun-Aug, Dec-Feb)
##################################################

# Extract lists of species found in each cluster
#Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
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

# Fill ms_pca_matrix_complete with values from clusters where rows match
ms_pca_matrix_complete <- matrix(NA, nrow = length(grid.names.all), ncol = ncol(clusters),
                                 dimnames = list(grid.names.all, colnames(clusters)))
ms_pca_matrix_complete[rownames(clusters), ] <- clusters

# read in random ratser and prep cluster raster
model_rast <- rast("./output/model_raster.tif")
clus_rast <- model_rast
values(clus_rast) <- ms_pca_matrix_complete[,19] #local maximum at k=8
clus_rast <- as.factor(clus_rast)
#convert to shp
clus_poly <- as.polygons(clus_rast)
clus_poly <- st_as_sf(clus_poly)
clus_poly <- fortify(clus_poly)


###############Warning: This takes a very long time
for (k in 1:12){
# select month
cat(k)
model_rast <- rast("./output/model_raster.tif")
mos <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

month <- readRDS(paste0("./output/monthly_ranges/",mos[k],".rds"))
month <- month[!names(month) %in% c("Camptorhynchus_labradorius", "Moringilla_nivalis.tif")]
species_names <- names(month)
#species_names <- species_names[!species_names %in% c("Camptorhynchus_labradorius", "Moringilla_nivalis.tif")]
#month <- month[!names(month) %in% c("Camptorhynchus_labradorius", "Moringilla_nivalis.tif")]

sp.list.mon <- list()

for (i in 1:length(species_names)){
#for (i in 1:10){
  #cat(i)
  sp_rast <- model_rast
  sp <- month[[i]]
  values(sp_rast) <- sp
  sp_pres <- list()
  
  for (j in 1:8){
    sp_pres_out <- tibble()
    sp_pres_out <- sp_pres_out %>%
      tibble::add_column(species = character(), cluster = numeric(), status = numeric())
    clus.j <- clus_poly %>% filter(model_raster == j)
    sp_pres.i <- extract(sp_rast, clus.j, na.rm = TRUE)
    sp_pres_out[1,1] <- names(month[i])
    sp_pres_out[1,2] <- j

    if (sum(sp_pres.i$model_raster)>0){
      sp_pres_out[1,3] <- 1
    } else {
      sp_pres_out[1,3] <- 0
    }
    
    sp_pres <- sp_pres %>%
      bind_rows(sp_pres_out)
  }
  sp.list.mon <- sp.list.mon %>%
    bind_rows(sp_pres)
}
    
saveRDS(sp.list.mon, paste0("./output/clusters_traitspace/sp_list_cluster_",mos[k],".rds"))  
#rm(list = ls())
}


### read in the output for each month
jan <- readRDS("./output/clusters_traitspace/sp_list_cluster_Jan.rds") 
feb <- readRDS("./output/clusters_traitspace/sp_list_cluster_Feb.rds") 
mar <- readRDS("./output/clusters_traitspace/sp_list_cluster_Mar.rds") 
apr <- readRDS("./output/clusters_traitspace/sp_list_cluster_Apr.rds") 
may <- readRDS("./output/clusters_traitspace/sp_list_cluster_May.rds") 
jun <- readRDS("./output/clusters_traitspace/sp_list_cluster_Jun.rds") 
jul <- readRDS("./output/clusters_traitspace/sp_list_cluster_Jul.rds") 
aug <- readRDS("./output/clusters_traitspace/sp_list_cluster_Aug.rds") 
sep <- readRDS("./output/clusters_traitspace/sp_list_cluster_Sep.rds") 
oct <- readRDS("./output/clusters_traitspace/sp_list_cluster_Oct.rds") 
nov <- readRDS("./output/clusters_traitspace/sp_list_cluster_Nov.rds") 
dec <- readRDS("./output/clusters_traitspace/sp_list_cluster_Dec.rds") 

all.months <- jan %>%
  bind_cols(feb[,3], mar[,3], apr[,3], may[,3],
            jun[,3], jul[,3], aug[,3], sep[,3],
            oct[,3], nov[,3], dec[,3])

colnames(all.months) <- c("species","cluster","Jan","Feb","Mar","Apr",
                          "May","Jun","Jul","Aug",
                          "Sep","Oct","Nov","Dec")

all.months.long <- pivot_longer(all.months, cols= c("Jan","Feb","Mar","Apr",
                                                   "May","Jun","Jul","Aug",
                                                   "Sep","Oct","Nov","Dec"),
                                names_to = "month", values_to = "status")


##plots for each cluster region
spp_all_scores <- read.csv("./output/species_diet_scores.csv") #This is a result of pca: clr_diet_pca$x
diet_loadings <- read.csv("./output/pca_variable_corr.csv") #this is the output of the pca: clr_diet_pca$rotation (eigenvectors)
all.months.long <- all.months.long %>%
  mutate(species_scientific_name = gsub("_", " ", species)) %>%
  mutate(id_uniq = paste(species, month, sep = "_")) %>%
  dplyr::select(species_scientific_name, month, species, id_uniq, cluster, status)

sp_all_scores_join <- spp_all_scores %>%
  inner_join(all.months.long, by = "id_uniq") %>%
  mutate(species_scientific_name = species_scientific_name.x, month = month.x, species = species.x) %>%
  dplyr::select(id_uniq, species_scientific_name, species, month, cluster, status, 
         invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
         PC1, PC2, PC3, PC4, PC5, PC6, PC7)

#select the diet axis that had a loading >=0.6
diet_load_pc12_sign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar")
diet_load_pc12_nonsign <- diet_loadings %>%
  filter(variable == "fruit" | variable == "scavenger" | variable == "other" | variable == "vertebrate")

#plot all points for pc 3 and 4
diet_load_pc34_sign <- diet_loadings %>%
  filter(variable == "other" | variable == "fruit")
diet_load_pc34_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar" | variable == "scavenger" | variable == "vertebrate")

#plot all points for pc 5 and 6
diet_load_pc56_sign <- diet_loadings %>%
  filter(variable == "scavenger" | variable == "vertebrate")
diet_load_pc56_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "other" | variable == "fruit" | variable == "seed" | variable == "nectar")


## plot for each cluster, looping across cluster 1:8
for (r in 1:8){
  #select from the scores for all species only those species that are present in a given region and time period  
  reg.i <- sp_all_scores_join[sp_all_scores_join$cluster == r,]  
  
  #now only select those species that are present, i.e., status=1, in a given season
  dec_feb <- reg.i %>%
    filter(month %in% c("Dec", "Jan", "Feb")) %>%
    filter(status == 1)
  mar_may <- reg.i %>%
    filter(month %in% c("Mar", "Apr", "May")) %>%
    filter(status == 1)
  jun_aug <- reg.i %>%
    filter(month %in% c("Jun", "Jul", "Aug")) %>%
    filter(status == 1)
  sep_nov <- reg.i %>%
    filter(month %in% c("Sep", "Oct", "Nov")) %>%
    filter(status == 1)
  

  
  ggsave(p, file=paste0("graphics/./dietspace_kde_SEASWINSUM_pc1_2_cluster",r,".png"), bg = "transparent", width=8, height=8, dpi=600)
  
  
  
  #plot all points for pc 3 and 4
  p <- ggplot() +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    #PC scores for each diet
    #winter
    geom_point(data = dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
    geom_hdr(data = dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
    geom_hdr_lines(data = dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
       #summer
    geom_point(data = jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
    geom_hdr(data = jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
    geom_hdr_lines(data = jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
    geom_segment(data = diet_load_pc34_sign,
                 aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2), colour = "black", size=2.5,
                 arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
     theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
      panel.border = element_rect(fill=NA, colour = "black", size=1.5),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=22),
      axis.text.y = element_text(colour='black',size=22),
      #axis.text.y = element_blank(),
      #axis.title.x = element_text(colour='black',size=22),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(colour='black',size=22),
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
  
  ggsave(p, file=paste0("./graphics/dietspace_kde_SEASWINSUM_pc3_4_cluster",r,".png"), bg = "transparent", width=8, height=8, dpi=600)
  
  
  p <- ggplot() +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    #PC scores for each diet
    #winter
    geom_point(data = dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
    geom_hdr(data = dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
    geom_hdr_lines(data = dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
    #summer
    geom_point(data = jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
    geom_hdr(data = jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
    geom_hdr_lines(data = jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
    geom_segment(data = diet_load_pc56_sign,
                 aes(x = 0, y = 0, xend = PC5*2, yend = PC6*2), colour = "black", size=2.5,
                 arrow = arrow(length = unit(0.75, "cm"), type = "closed")) +
     theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
      panel.border = element_rect(fill=NA, colour = "black", size=1.5),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=22),
      axis.text.y = element_text(colour='black',size=22),
      #axis.text.y = element_blank(),
      #axis.title.x = element_text(colour='black',size=22),
      axis.title.x = element_blank(),
      #axis.title.y = element_text(colour='black',size=22),
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
  
  ggsave(p, file=paste0("./graphics/dietspace_kde_SEASWINSUM_pc5_6_cluster",r,".png"), bg = "transparent", width=8, height=8, dpi=600)
  
}


  