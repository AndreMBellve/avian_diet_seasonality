#Script to plot dietary space across resident, short-distance, and long-distance migrants, and plot species-level diet variability

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

setwd("/Users/jarzyna.1/Documents/Dropbox/RESEARCH_TemporallyVaryingDiet")

##################################################
##  Step 1: histograms of dietary variability
##################################################
merged_df_all <- read.csv(file = "./output/species_diet_variability.csv")

res <- mean(merged_df_all$summed_wgt_var[merged_df_all$resident_only == TRUE], na.rm = TRUE)
mig <- mean(merged_df_all$summed_wgt_var[merged_df_all$resident_only == FALSE], na.rm = TRUE)
shor <- mean(merged_df_all$summed_wgt_var[merged_df_all$mig_type == 1], na.rm = TRUE) #short distance migrants
long <- mean(merged_df_all$summed_wgt_var[merged_df_all$mig_type == 2], na.rm = TRUE) #long distance migrants

#plot for all species, including the ones that showed no variability
p <- ggplot(merged_df_all, aes(x = summed_wgt_var, fill = factor(mig_type), color = factor(mig_type))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  scale_color_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  geom_vline(xintercept = res, linetype = "twodash", color = "grey60", size=1.5) +
  geom_vline(xintercept = shor, linetype = "twodash", color = "skyblue1", size=1.5) +
  geom_vline(xintercept = long, linetype = "twodash", color = "indianred1", size=1.5) +
  labs(x = "Diet variability",y = "Density") +
  theme(
    #panel.background = element_rect(fill = 'white', colour = 'grey85'),
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
ggsave(p, file="./graphics/dietvar_all.png", bg="transparent", width=8, height=6, dpi=600)

#cut of x axis at 0.2
p <- ggplot(merged_df_all, aes(x = summed_wgt_var, fill = factor(mig_type), color = factor(mig_type))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  scale_color_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  geom_vline(xintercept = res, linetype = "twodash", color = "grey60", size=1.5) +
  geom_vline(xintercept = shor, linetype = "twodash", color = "skyblue1", size=1.5) +
  geom_vline(xintercept = long, linetype = "twodash", color = "indianred1", size=1.5) +
  labs(x = "Diet variability",y = "Density") +
  scale_x_continuous(limits = c(0, 0.2)) + 
  theme(
    #panel.background = element_rect(fill = 'white', colour = 'grey85'),
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
ggsave(p, file="./graphics/dietvar_all_01_V3.png", bg="transparent", width=8, height=6, dpi=600)


#plot for species that showed variability >0
merged_df_var <- merged_df_all %>%
  filter(summed_wgt_var > 0)
res <- mean(merged_df_var$summed_wgt_var[merged_df_var$resident_only == TRUE], na.rm = TRUE)
mig <- mean(merged_df_var$summed_wgt_var[merged_df_var$resident_only == FALSE], na.rm = TRUE)
shor <- mean(merged_df_var$summed_wgt_var[merged_df_var$mig_type == 1], na.rm = TRUE) #short distance migrants
long <- mean(merged_df_var$summed_wgt_var[merged_df_var$mig_type == 2], na.rm = TRUE) #long distance migrants

p <- ggplot(merged_df_var, aes(x = summed_wgt_var, fill = factor(mig_type), color = factor(mig_type))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  scale_color_manual(values = c("0" = "grey60", "1" = "skyblue1", "2" = "indianred1")) +
  geom_vline(xintercept = res, linetype = "twodash", color = "grey60", size=1.5) +
  geom_vline(xintercept = shor, linetype = "twodash", color = "skyblue1", size=1.5) +
  geom_vline(xintercept = long, linetype = "twodash", color = "indianred1", size=1.5) +
  labs(x = "Diet variability",y = "Density") +
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
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
ggsave(p, file="./graphics/dietvar_varonly_V3.png", bg="transparent", width=8, height=6, dpi=600)


##################################################
##  Step 2: latitudinal variation
##################################################
merged_df_all <- merged_df_all %>%
  mutate(absmedianlat = abs(medianlat))

resident <- merged_df_all[merged_df_all$resident_only == TRUE,]
migrant <- merged_df_all[merged_df_all$resident_only == FALSE,]
short <- merged_df_all[merged_df_all$mig_type == 1,] #short distance migrants
long <- merged_df_all[merged_df_all$mig_type == 2,] #long distance migrants

summary(lm(migrant$summed_wgt_var~migrant$difflat))##this is what we're reporting, no apparent relationship with migration distance
#summary(lm(migrant$summed_wgt_var~migrant$difflat + I(migrant$difflat^2)))

p <- ggplot(migrant, aes(x = difflat, y = summed_wgt_var)) +
  geom_point(shape = 21, alpha=0.6, colour="black", fill="grey80", size=5) +
  #geom_smooth() +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 1)+
  labs(x = "Migration distance",y = "Diet variability") +
  #scale_y_continuous(limits=c(0,6.2))+
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/dietvar_migrationdistance_allsp.png", bg="transparent", width=8, height=6, dpi=600)


#repeat for only species that have variable diet
merged_df_var <- merged_df_all %>%
  filter(summed_wgt_var > 0)
resident <- merged_df_var[merged_df_var$resident_only == TRUE,]
migrant <- merged_df_var[merged_df_var$resident_only == FALSE,]
short <- merged_df_var[merged_df_var$mig_type == 1,] #short distance migrants
long <- merged_df_var[merged_df_var$mig_type == 2,] #long distance migrants

summary(lm(migrant$summed_wgt_var~migrant$difflat)) #this is what we're reporting
summary(lm(migrant$summed_wgt_var~migrant$difflat + I(migrant$difflat^2)))

p <- ggplot(migrant, aes(x = difflat, y = summed_wgt_var)) +
  geom_point(shape = 21, alpha=0.6, colour="black", fill="grey80", size=5) +
  #geom_smooth() +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 1)+
  labs(x = "Migration distance",y = "Diet variability") +
  #scale_y_continuous(limits=c(0,6.2))+
  theme(#panel.background = element_rect(fill = 'white', colour = 'grey85'),
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

ggsave(p, file="./graphics/dietvar_migrationdistance_varonlysp.png", bg="transparent", width=8, height=6, dpi=600)


##################################################
##  Step 3: Diet space for resident versus migrants
##################################################
### let's plot the trait space for resident and migrant (short- and long-distance) birds
spp_all_scores <- read.csv("./output/species_diet_scores.csv") #This is a result of pca: clr_diet_pca$x
diet_loadings <- read.csv("./output/pca_variable_corr.csv") #this is the output of the pca: clr_diet_pca$rotation (eigenvectors)

#pc 1 and 2
diet_load_pc12_sign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar")
diet_load_pc12_nonsign <- diet_loadings %>%
  filter(variable == "fruit" | variable == "scavenger" | variable == "other" | variable == "vertebrate")

#pc 3 and 4
diet_load_pc34_sign <- diet_loadings %>%
  filter(variable == "other" | variable == "fruit")
diet_load_pc34_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "seed" | variable == "nectar" | variable == "scavenger" | variable == "vertebrate")

#pc 5 and 6
diet_load_pc56_sign <- diet_loadings %>%
  filter(variable == "scavenger" | variable == "vertebrate")
diet_load_pc56_nonsign <- diet_loadings %>%
  filter(variable == "invertebrate" | variable == "other" | variable == "fruit" | variable == "seed" | variable == "nectar")


#prepare plots for all species, regardless if they showed diet variability
resident <- merged_df_all[merged_df_all$resident_only == TRUE,]
migrant <- merged_df_all[merged_df_all$resident_only == FALSE,]
short <- merged_df_all[merged_df_all$mig_type == 1,] #short distance migrants
long <- merged_df_all[merged_df_all$mig_type == 2,] #long distance migrants


# residents
res_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% resident[, "species_scientific_name"], ]
# all migrants
mig_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% migrant[, "species_scientific_name"], ]
# short-distance migrants
short_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% short[, "species_scientific_name"], ]
# long-distance migrants
long_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% long[, "species_scientific_name"], ]

res_dec_feb <- res_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
res_mar_may <- res_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
res_jun_aug <- res_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
res_sep_nov <- res_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

mig_dec_feb <- mig_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
mig_mar_may <- mig_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
mig_jun_aug <- mig_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
mig_sep_nov <- mig_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

short_dec_feb <- short_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
short_mar_may <- short_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
short_jun_aug <- short_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
short_sep_nov <- short_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

long_dec_feb <- long_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
long_mar_may <- long_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
long_jun_aug <- long_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
long_sep_nov <- long_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

winter_pal <- "#180F3EFF"
summer_pal <- "#FFC125"

density_probs  <-  c(0.9, 0.8, 0.7)


##plot diet space for resident birds
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_resident.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4_resident.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_resident.png"), bg = "transparent", width=8, height=8, dpi=600)


##for short-distance migrants
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_short-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4_short-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_short-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


##for long-distance migrants
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_long-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4_long-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_long-migrant.png"), bg = "transparent", width=8, height=8, dpi=600)


##################################################
##  Step 4: Diet space for resident versus migrants, but for species with variable diets only
##################################################
merged_df_var <- merged_df_all %>%
  filter(summed_wgt_var > 0)
resident <- merged_df_var[merged_df_var$resident_only == TRUE,]
migrant <- merged_df_var[merged_df_var$resident_only == FALSE,]
short <- merged_df_var[merged_df_var$mig_type == 1,] #short distance migrants
long <- merged_df_var[merged_df_var$mig_type == 2,] #long distance migrants

#for residents
res_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% resident[, "species_scientific_name"], ]
#for migrants
mig_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% migrant[, "species_scientific_name"], ]
#for short-distance migrants
short_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% short[, "species_scientific_name"], ]
#for long-distance migrants
long_var_scores <- spp_all_scores[spp_all_scores[, "species_scientific_name"] %in% long[, "species_scientific_name"], ]

res_dec_feb <- res_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
res_mar_may <- res_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
res_jun_aug <- res_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
res_sep_nov <- res_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

mig_dec_feb <- mig_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
mig_mar_may <- mig_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
mig_jun_aug <- mig_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
mig_sep_nov <- mig_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

short_dec_feb <- short_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
short_mar_may <- short_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
short_jun_aug <- short_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
short_sep_nov <- short_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

long_dec_feb <- long_var_scores %>%
  filter(month %in% c("Dec", "Jan", "Feb"))
long_mar_may <- long_var_scores %>%
  filter(month %in% c("Mar", "Apr", "May"))
long_jun_aug <- long_var_scores %>%
  filter(month %in% c("Jun", "Jul", "Aug"))
long_sep_nov <- long_var_scores %>%
  filter(month %in% c("Sep", "Oct", "Nov"))

winter_pal <- "#180F3EFF"
summer_pal <- "#FFC125"

density_probs  <-  c(0.9, 0.8, 0.7)


#for resident birds with dietary variability
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_resident_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4_resident_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = res_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = res_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = res_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = res_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = res_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_resident_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


#for short-distance miigrants with dietary variability
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_short-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc3_4_short-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = short_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = short_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = short_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = short_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = short_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_short-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


#for long-distance migrants with dietary varibaility
#plot all points for pc 1 and 2
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC1,y = PC2), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC1, y = PC2), method = "kde", probs = density_probs, colour = summer_pal) +
  #loadings
  geom_segment(data = diet_load_pc12_sign,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), colour = "black", size=2.5,
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc1_2_long-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


#plot all points for pc 3 and 4
p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC3,y = PC4), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC3, y = PC4), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("graphics/dietspace_kde_pc3_4_long-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)


p <- ggplot() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #PC scores for each diet
  #winter
  geom_point(data = long_dec_feb, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = winter_pal) +
  geom_hdr_lines(data = long_dec_feb, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = winter_pal) +
  #summer
  geom_point(data = long_jun_aug, aes(x = PC5,y = PC6), alpha=1, colour="grey70", size=1.5, position="jitter") +
  geom_hdr(data = long_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, fill = summer_pal) +
  geom_hdr_lines(data = long_jun_aug, aes(x = PC5, y = PC6), method = "kde", probs = density_probs, colour = summer_pal) +
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

ggsave(p, file=paste0("./graphics/dietspace_kde_pc5_6_long-migrant_varonly.png"), bg = "transparent", width=8, height=8, dpi=600)
