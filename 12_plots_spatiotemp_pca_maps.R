#Author: Marta Jarzyna

#Script to plot space-time PCA results
require(tidyverse)
require(ggplot2)
require(dplyr)
require(tidyr)
library(terra)
library(tidyterra)
require(stars)

##no needed:
require(here)
require(svglite)
require(viridis)
require(lubridate)
require(scales)
require(stars)
require(raster)
require(rnaturalearth)

require(grid)
require(maptools)
require(rgdal)
require(sf)


### This is a frustrating bit about R. Two packages both want to use the command select. Make sure we are using the correct one
select <- dplyr::select

##################################################
##  Read PCA data and output in
##################################################
ms_pca_matrix <- readRDS(file = "./output/month_space_pca/ms_pca_matrix.rds")
pca_fit <- readRDS(file = "./output/month_space_pca/ts_pca.rds")
pc_importance <- summary(pca_fit)
pc_importance


##################################################
##  Scree Plot
##################################################
#write_path <- write_figures_path
#file.path(write_figures_path, "fric")
#dir.create(write_path, recursive=TRUE, showWarnings = FALSE)

### Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:12,]

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
	+ geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
	+ scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
	+ scale_y_continuous(name = "Proportion of variance explained (%)", limits=c(0,80),breaks=c(0,10,20,30,40,50,60,70,80)) %>%
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

ggsave(p, file="./graphics/spacetime-PCA_prop_var.png", bg="transparent", width=8, height=6, dpi=600)

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Cumulative variance explained (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
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

ggsave(p, file="./graphics/spacetime-PCA_cum_var.png", bg="transparent", width=8, height=6, dpi=600)



##################################################
##  Plot Spatiotemporal PC scores
##################################################

### Loop through scores
for(i in seq(1, 3)){
  
  ### Create scores
  plot_df <- data.frame(month = seq(1,12), score = pca_fit$x[,i])
  #plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
  
  ### Create plot
  p <- ggplot(plot_df, aes(x=month, y=score)) %>%
    + geom_line(size=2) %>%
    + geom_point(alpha=1,size=4, pch=16) %>%
    #+ geom_hline(yintercept = 0, linetype = "dashed") %>%
    + geom_hline(yintercept = 0, linetype = "twodash", color = "darkmagenta", size=1.5) %>%
    + scale_x_continuous(name = "Month", breaks = c(2,4,6, 8,10,12)) %>%
    + scale_y_continuous(name = "Score") %>%
    + theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
            panel.border = element_rect(fill=NA, colour = "white", size=1),
            axis.line = element_line(color = 'black', size=1.5),
            plot.title = element_text(size=15, vjust=2, family="sans"),
            axis.text.x = element_text(colour='black',size=30),
            axis.text.y = element_text(colour='black',size=30),
           # axis.title.x = element_text(colour='black',size=22),
           #axis.title.y = element_text(colour='black',size=22),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_line(color = 'black', size=1.5),
            axis.ticks.length=unit(0.3,"cm"),
            legend.position="none",
            legend.text=element_text(size=20),
            legend.title=element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank())
  p
  
  ### Save plot
  ggsave(p, file=paste0("./graphics/spacetime_pc_",i,"_score.png"), width=8, height=6, dpi=600)
  
}

##get all scores for all PCs
### Create scores
plot_df <- data.frame(month = seq(1,12), pca_fit$x)
library(tidyr)
df_long <- pivot_longer(plot_df, # Select all columns except 'row'
  cols = colnames(plot_df[,2:13]),
  names_to = "PC", # Name of the new column for former column names
  values_to = "score"  # Name of the new column for cell values
)

df_long$PC <- factor(df_long$PC, levels = unique(df_long$PC))

p <- ggplot(df_long, aes(x = month, y = score)) +
  geom_line(size = 2) +
  geom_point(alpha = 1, size = 4, pch = 16) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "darkmagenta", size = 1.5) +
  scale_x_continuous(name = "Month", breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(name = "Score") +
  facet_wrap(~PC, ncol=4) +  # Add faceting by column 'PC'
  theme(
    panel.border = element_rect(fill = NA, colour = "white", size = 1),
    axis.line = element_line(color = 'black', size = 1.5),
    plot.title = element_text(size = 15, vjust = 2, family = "sans"),
    axis.text.x = element_text(colour = 'black', size = 30),
    axis.text.y = element_text(colour = 'black', size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(color = 'black', size = 1.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.text = element_text(size = 30),
    legend.position = "none",
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
p

### Save plot
ggsave(p, file=paste0("./graphics/spacetime_pc_score_allpcs.png"), bg="transparent", width=24, height=18, dpi=600)
  
  
##################################################
##  Plot Spatiotemporal PC loadings (maps)
##################################################
#Loading in data from script 04
ts_pca <- readRDS(file="./output/month_space_pca/ts_pca.rds")
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 
pca_geography_matrix <- readRDS("./output/pca_geography_matrix.rds")

#geography pca has 80000*6 chunks The first 80000 is pc1
ms_pca_mat_t <- t(ms_pca_matrix)

#Data frame of row indices to iterate through
row_index_df <- data.frame(begin = c(1, 80001, 160001,240001, 320001, 400001),
                           end = c(80000, 160000, 240000, 320000, 400000, 480000))

# Get the continent shapefile
globe_shp <- ne_download(scale = "medium", type = "land",
                         category = "physical", returnclass = "sf")
globe_shp <- fortify(globe_shp)

ocean_shp <- ne_download(scale = "medium", type = "ocean",
                         category = "physical", returnclass = "sf")
ocean_shp <- fortify(ocean_shp)

#Iterating through to build the rasters for each of our original six PCA axes - these are split into  80k chunks
for(pc in 1:6){
  
  cur_rast_vals <- ms_pca_mat_t[seq(from = row_index_df[pc,"begin"],
                                    to = row_index_df[pc,"end"],
                                    by = 1),]
  
  #load in raster and make a model raster
  model_rast <- rast("./output/model_raster.tif")
  #Overwriting assemblage values from the model raster as some sites were not fed into the PCA
  values(model_rast) <- NA
  
  #create 2 replicate layers 
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



#Testing to see what range and breaks makes sense for each of the ts-PCA discretizations
#First three axes look like they have similar bounds!
pc1_val <- values(ts_rast[[1]]) %>%  na.omit()
hist(pc1_val, xlim = c(-0.02, 0.02), breaks = 10000)

pc2_val <- values(ts_rast[[2]]) %>%  na.omit()
hist(pc2_val, xlim = c(-0.02, 0.02), breaks = 10000)


### This is the full grid, notice everything goes out to -180, 180, -90 and 90
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


p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  #scale_fill_distiller(palette = "RdBu", type="div",
  #                     limits=c(-0.01,0.01), oob=squish,
  #                     direction = -1, 
  #                     name = "PC-Trait Loading",
  #                     na.value = "transparent") + 
  scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
                       limits = c(-0.015, 0.015), oob = squish, name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + 
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 2) + 
  coord_sf(crs = "ESRI:54030") 

ggsave("./graphics/time_space_diet_rasters_grids_015.png", bg = "transparent", width = 10, height = 16, dpi=800)

p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  #scale_fill_distiller(palette = "RdBu", type="div",
  #                     limits=c(-0.01,0.01), oob=squish,
  #                     direction = -1, 
  #                     name = "PC-Trait Loading",
  #                     na.value = "transparent") + 
  scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
                       limits = c(-0.015, 0.015), oob = squish, name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + 
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="right",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 2) + 
  coord_sf(crs = "ESRI:54030") 

ggsave("./graphics/time_space_diet_rasters_grids_015_legend.png", bg = "transparent", width = 10, height = 16, dpi=800)


##################################################
##  Plot Dominant Spatiotemporal PC (maps)
##################################################
#Loading in data from script 04
ms_pca_matrix <- readRDS("./output/month_space_pca/ms_pca_matrix.rds") 

#geography pca has 80000*6 chunks The first 80000 is pc1
ms_pca_mat_t <- t(ms_pca_matrix)

#designate which PC has a higher absolute loading 
#(i.e., presumably contributes more to explaining temporal variation in each diet axes)
#this can be done like this because loadings for both PCs are on pretty much the same scale
ms_pca_mat_t <- as.data.frame(ms_pca_mat_t)
# 1: PC1 > PC2, and is positive
# 2: PC1 > PC2, and is negative
# 3: PC2 > PC1, and is positive
# 4: PC2 > PC1, and is negative
ms_pca_mat_t$Contr <- ifelse(abs(ms_pca_mat_t$PC1) > abs(ms_pca_mat_t$PC2) & ms_pca_mat_t$PC1 > 0, 1,
                             ifelse(abs(ms_pca_mat_t$PC1) > abs(ms_pca_mat_t$PC2) & ms_pca_mat_t$PC1 < 0, 2,
                                    ifelse(abs(ms_pca_mat_t$PC2) > abs(ms_pca_mat_t$PC1) & ms_pca_mat_t$PC2 > 0, 3, 
                                           ifelse(abs(ms_pca_mat_t$PC2) > abs(ms_pca_mat_t$PC1) & ms_pca_mat_t$PC2 < 0, 4, NA))))


#Iterating through to build the rasters for each of our original six PCA axes - these are split into  80k chunks
for(pc in 1:6){
  
  cur_rast_vals <- ms_pca_mat_t[seq(from = row_index_df[pc,"begin"],
                                    to = row_index_df[pc,"end"],
                                    by = 1),]
  cur_rast_vals <- cur_rast_vals[,3]
  
  #load in raster and make a model raster
  model_rast <- rast("./output/model_raster.tif")
  #Overwriting assemblage values from the model raster as some sites were not fed into the PCA
  values(model_rast) <- cur_rast_vals
  
  #relabel names
  pc_names <- paste0(c('Contr'), "-PC", pc)
  
  #Naming the rasters for clarity
  names(model_rast) <- pc_names
  
  #Either initialising a spatraster storage object or 
  if(pc == 1){
    ts_rast <- model_rast
  }else{
    ts_rast <- c(ts_rast, model_rast)
  }
}


#plot
p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  scale_fill_gradientn(colors = c("maroon2", "pink", "#006400", "#90EE90"), na.value = "transparent", name = "Contribution") +
  #scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
  #                    limits = c(-0.015, 0.015), oob = squish, name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 1) + 
  coord_sf(crs = "ESRI:54030") 

ggsave("./graphics/time_space_diet_pc-contr.png", bg = "transparent", width = 5, height = 16, dpi=800)

p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  scale_fill_gradientn(colors = c("maroon2", "pink", "#006400", "#90EE90"), na.value = "transparent", name = "Contribution") +
  #scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
  #                    limits = c(-0.015, 0.015), oob = squish, name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 1) + 
  coord_sf(crs = "ESRI:54030") 

ggsave("./graphics/time_space_diet_pc-contr-legend.png", bg = "transparent", width = 5, height = 16, dpi=800)



################################################################
##  Sensitivity analysis
##################################################

##########################################################25% least certain species removed
##################################################
##  Read PCA data and output in
##################################################
pca_fit <- readRDS(file = "./output/certainty_analysis/month_space_pca/ts_25_pca.rds")
ms_pca_matrix <- readRDS(file = "./output/certainty_analysis/month_space_pca/ms_pca_25_matrix.rds")
pc_importance <- summary(pca_fit)
pc_importance


##################################################
##  Scree Plot
##################################################
#write_path <- write_figures_path
#file.path(write_figures_path, "fric")
#dir.create(write_path, recursive=TRUE, showWarnings = FALSE)

### Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:12,]

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
  + scale_y_continuous(name = "Proportion of variance explained (%)", limits=c(0,80),breaks=c(0,10,20,30,40,50,60,70,80)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/spacetime-PCA_prop_var_25.png", width=8, height=6, dpi=600)

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Cumulative variance explained (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/spacetime-PCA_cum_var_25.png", width=8, height=6, dpi=600)



##################################################
##  Plot Spatiotemporal PC scores
##################################################

### Loop through scores
for(i in seq(1, 3)){
  
  #pc_folder <- file.path(write_path, paste0("pc_", i))
  #dir.create(pc_folder, recursive=TRUE, showWarnings = FALSE)
  
  ### Create scores
  plot_df <- data.frame(month = seq(1,12), score = pca_fit$x[,i])
  #plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
  
  ### Create plot
  p <- ggplot(plot_df, aes(x=month, y=score)) %>%
    + geom_line(size=2) %>%
    + geom_point(alpha=1,size=4, pch=16) %>%
    #+ geom_hline(yintercept = 0, linetype = "dashed") %>%
    + geom_hline(yintercept = 0, linetype = "twodash", color = "darkmagenta", size=1.5) %>%
    + scale_x_continuous(name = "Month", breaks = c(2,4,6, 8,10,12)) %>%
    + scale_y_continuous(name = "Score") %>%
    + theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
      panel.border = element_rect(fill=NA, colour = "white", size=1),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=30),
      axis.text.y = element_text(colour='black',size=30),
      # axis.title.x = element_text(colour='black',size=22),
      #axis.title.y = element_text(colour='black',size=22),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_line(color = 'black', size=1.5),
      axis.ticks.length=unit(0.3,"cm"),
      legend.position="none",
      legend.text=element_text(size=20),
      legend.title=element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())
  p
  
  ### Save plot
  ggsave(p, file=paste0("./graphics/spacetime_pc_",i,"_score_25.png"), width=8, height=6, dpi=600)
  
}

##################################################
##  Plot Spatiotemporal PC loadings (maps)
##################################################
#Loading in data from script 04
ts_pca <- readRDS(file = "./output/certainty_analysis/month_space_pca/ts_25_pca.rds")
ms_pca_matrix <- readRDS(file = "./output/certainty_analysis/month_space_pca/ms_pca_25_matrix.rds")
pca_geography_matrix <- readRDS("./output/certainty_analysis/pca_25_geography_matrix.rds")

#geography pca has 80000*6 chunks The first 80000 is pc1
ms_pca_mat_t <- t(ms_pca_matrix)

#Data frame of row indices to iterate through
ms_pca_mat_t <- ms_pca_mat_t[,1:2]
row_index_df <- data.frame(begin = c(1, 80001, 160001,240001, 320001, 400001),
                           end = c(80000, 160000, 240000, 320000, 400000, 480000))

# Get the continent shapefile
globe_shp <- ne_download(scale = "medium", type = "land",
                         category = "physical", returnclass = "sf")
globe_shp <- fortify(globe_shp)

ocean_shp <- ne_download(scale = "medium", type = "ocean",
                         category = "physical", returnclass = "sf")
ocean_shp <- fortify(ocean_shp)

#Iterating through to build the rasters for each of our original six PCA axes - these are split into  80k chunks
for(pc in 1:6){
  
  cur_rast_vals <- ms_pca_mat_t[seq(from = row_index_df[pc,"begin"],
                                    to = row_index_df[pc,"end"],
                                    by = 1),]
  
  #load in raster and make a model raster
  model_rast <- rast("./output/model_raster.tif")
  #Overwriting assemblage values from the model raster as some sites were not fed into the PCA
  values(model_rast) <- NA
  
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



#Testing to see what range and breaks makes sense for each of the ts-PCA discretizations
#First three axes look like they have similar bounds!
pc1_val <- values(ts_rast[[1]]) %>%  na.omit()
hist(pc1_val, xlim = c(-0.015, 0.015), breaks = 10000)

pc2_val <- values(ts_rast[[2]]) %>%  na.omit()
hist(pc2_val, xlim = c(-0.015, 0.015), breaks = 10000)


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


p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  #scale_fill_distiller(palette = "RdBu", type="div",
  #                     limits=c(-0.01,0.01), oob=squish,
  #                     direction = -1, 
  #                     name = "PC-Trait Loading",
  #                     na.value = "transparent") + 
  scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
                       limits = c(-0.015, 0.015), oob = squish, 
                       name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 2) + 
  coord_sf(crs = "ESRI:54030") 

#Saving the PCA output
ggsave("./graphics/time_space_diet_rasters_grids_015_25.png", bg = "transparent", width = 10, height = 16, dpi=800)



##########################################################50% least certain species removed

##################################################
##  Read PCA data and output in
##################################################
pca_fit <- readRDS(file = "./output/certainty_analysis/month_space_pca/ts_50_pca.rds")
ms_pca_matrix <- readRDS(file = "./output/certainty_analysis/month_space_pca/ms_pca_50_matrix.rds")
pc_importance <- summary(pca_fit)
pc_importance


##################################################
##  Scree Plot
##################################################
#write_path <- write_figures_path
#file.path(write_figures_path, "fric")
#dir.create(write_path, recursive=TRUE, showWarnings = FALSE)

### Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:12,]

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
  + scale_y_continuous(name = "Proportion of variance explained (%)", limits=c(0,80),breaks=c(0,10,20,30,40,50,60,70,80)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/spacetime-PCA_prop_var_50.png", width=8, height=6, dpi=600)

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Cumulative variance explained (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/spacetime-PCA_cum_var_50.png", width=8, height=6, dpi=600)



##################################################
##  Plot Spatiotemporal PC scores
##################################################

### Loop through scores
for(i in seq(1, 3)){
  
  #pc_folder <- file.path(write_path, paste0("pc_", i))
  #dir.create(pc_folder, recursive=TRUE, showWarnings = FALSE)
  
  ### Create scores
  plot_df <- data.frame(month = seq(1,12), score = pca_fit$x[,i])
  #plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
  
  ### Create plot
  p <- ggplot(plot_df, aes(x=month, y=score)) %>%
    + geom_line(size=2) %>%
    + geom_point(alpha=1,size=4, pch=16) %>%
    #+ geom_hline(yintercept = 0, linetype = "dashed") %>%
    + geom_hline(yintercept = 0, linetype = "twodash", color = "darkmagenta", size=1.5) %>%
    + scale_x_continuous(name = "Month", breaks = c(2,4,6, 8,10,12)) %>%
    + scale_y_continuous(name = "Score") %>%
    + theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
      panel.border = element_rect(fill=NA, colour = "white", size=1),
      axis.line = element_line(color = 'black', size=1.5),
      plot.title = element_text(size=15, vjust=2, family="sans"),
      axis.text.x = element_text(colour='black',size=30),
      axis.text.y = element_text(colour='black',size=30),
      # axis.title.x = element_text(colour='black',size=22),
      #axis.title.y = element_text(colour='black',size=22),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_line(color = 'black', size=1.5),
      axis.ticks.length=unit(0.3,"cm"),
      legend.position="none",
      legend.text=element_text(size=20),
      legend.title=element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())
  p
  
  ### Save plot
  ggsave(p, file=paste0("./graphics/spacetime_pc_",i,"_score_50.png"), width=8, height=6, dpi=600)
  
}

##################################################
##  Plot Spatiotemporal PC loadings (maps)
##################################################
#Loading in data from script 04
ts_pca <- readRDS(file = "./output/certainty_analysis/month_space_pca/ts_50_pca.rds")
ms_pca_matrix <- readRDS(file = "./output/certainty_analysis/month_space_pca/ms_pca_50_matrix.rds")
pca_geography_matrix <- readRDS("./output/certainty_analysis/pca_50_geography_matrix.rds")

#geography pca has 80000*6 chunks The first 80000 is pc1
ms_pca_mat_t <- t(ms_pca_matrix)

#Data frame of row indices to iterate through
ms_pca_mat_t <- ms_pca_mat_t[,1:2]
row_index_df <- data.frame(begin = c(1, 80001, 160001,240001, 320001, 400001),
                           end = c(80000, 160000, 240000, 320000, 400000, 480000))

# Get the continent shapefile
globe_shp <- ne_download(scale = "medium", type = "land",
                         category = "physical", returnclass = "sf")
globe_shp <- fortify(globe_shp)

ocean_shp <- ne_download(scale = "medium", type = "ocean",
                         category = "physical", returnclass = "sf")
ocean_shp <- fortify(ocean_shp)

#Iterating through to build the rasters for each of our original six PCA axes - these are split into  80k chunks
for(pc in 1:6){
  
  cur_rast_vals <- ms_pca_mat_t[seq(from = row_index_df[pc,"begin"],
                                    to = row_index_df[pc,"end"],
                                    by = 1),]
  
  #load in raster and make a model raster
  model_rast <- rast("./output/model_raster.tif")
  #Overwriting assemblage values from the model raster as some sites were not fed into the PCA
  values(model_rast) <- NA
  
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



#Testing to see what range and breaks makes sense for each of the ts-PCA discretizations
#First three axes look like they have similar bounds!
pc1_val <- values(ts_rast[[1]]) %>%  na.omit()
hist(pc1_val, xlim = c(-0.015, 0.015), breaks = 10000)

pc2_val <- values(ts_rast[[2]]) %>%  na.omit()
hist(pc2_val, xlim = c(-0.015, 0.015), breaks = 10000)

# for plotting
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


p <- ggplot() + 
  geom_spatraster(data = ts_rast,
                  maxcell = 80000) +
  #scale_fill_distiller(palette = "RdBu", type="div",
  #                     limits=c(-0.01,0.01), oob=squish,
  #                     direction = -1, 
  #                     name = "PC-Trait Loading",
  #                     na.value = "transparent") + 
  scale_fill_gradient2(low = "#076ba7", mid = "white", high = "#ca2207", midpoint = 0, 
                       limits = c(-0.015, 0.015), oob = squish, 
                       name = "PC-Trait Loading", na.value = "transparent") + 
  geom_sf(data = ocean_shp, fill = "grey70", alpha = 1, size=0.2) +
  geom_sf(data = grid_full, linewidth = 0.2, colour = "grey40") + ### Swithc this out teo be grid_full
  geom_sf(data = globe_shp, fill = NA, alpha = 1, size=0.2) +
  #geom_sf(data = circle, fill = NA, color = "red", size = 1) +
  #theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",    
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA), 
        strip.text = element_blank()) +
  facet_wrap(~lyr, ncol = 2) + 
  coord_sf(crs = "ESRI:54030") 

#Saving the PCA output
ggsave("./graphics/time_space_diet_rasters_grids_015_50.png", bg = "transparent", width = 10, height = 16, dpi=800)



