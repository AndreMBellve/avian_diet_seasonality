#This script fits a PCA to CLR transformed diet data for the unique diets which occur among Avian species throughout the year. It then joins these back onto the full species × month dataframe, so that these values can subsequently be projected onto space and time.

#Author: Andre Bellve

# Libraries ---------------------------------------------------------------

#Data manipulation
require(dplyr)
require(tidyr)
require(stringr)
require(janitor)

#Analysis packages
require(compositions)
require(zCompositions)

#Visualisation
require(ggplot2)

# Data preparation --------------------------------------------------------

#Note, these data need to be downloaded from the SAviTraits publication (DOI: doi.org/10.1111/geb.13738)
#Cluster directory
raw_trait_df <- read.csv('./data/diet_data/SAviTraits_1-0_1.csv') %>%
  #Trimming out unnecessary variables
  dplyr::select(Species_Scientific_Name, Diet_Sub_Cat,
                Jan, Feb, Mar, Apr, May, Jun,
                Jul, Aug, Sep, Oct, Nov, Dec)


#Collapses the vertebrate diet groups into a single vertebrate diet class to reduce the number of traits analysed and because most of the diet proportions lie in the "unknown vertebrate category anyway
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

#Unique diets to create diet space from these
uni_diet_df <- full_trait_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) 
nrow(uni_diet_df)
#Setting seed so that the multiplicative simple replacement is reproducible
set.seed(42)
#Collapsing to unique entries, effectively ignoring species traits - this will allow us to calculate positions in the ordination space, which we can then use to map onto functional space
clr_diet_df <- full_trait_df %>%
  
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

#Invertebrates looks good...
ggplot(corr_test_df,
       aes(y = untransformed,
           x = clr)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2)) +
  facet_wrap(~diet)

#Generally all the values are monotonically transformed

#Fitting a PCA to the centered log ratio transformed data
clr_diet_pca <- prcomp(clr_diet_df,
                       center = FALSE)
#The complete diet PCA
saveRDS(clr_diet_pca, 
        "./output/diet_pca/clr_diet_pca.RDS")

#Variance explained by each PC axis
plot(clr_diet_pca)
round(((clr_diet_pca$sdev^2)/sum(clr_diet_pca$sdev^2)), 
      digits = 2)
#First two axes explain 58% of the variation in the data - first three explain 75%

# Species PC scores -------------------------------------------------------

#Creating a df that has the unique CLR transformed diets and their PC scores to join back onto the species dataframe
diet_scores_df <- full_trait_df %>%
  distinct(invertebrate, scavenger, fruit, nectar, seed, other, vertebrate,
           .keep_all = TRUE) %>% 
  dplyr::select(-c(species_scientific_name, species, id_uniq, month)) %>% 
  as.data.frame() %>% 
  bind_cols(clr_diet_pca$x)

#Preparing the full data frame to have species scores joined on so we can map these back onto the world through space and time
spp_scores_df <- full_trait_df %>% 
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
          "./output/species_diet_scores.csv",
          row.names = FALSE)

#Creating a stripped back wide version of the species × month PC scores so we can visualise the diet space for different biogeographic regions
spp_month_pc_df <- spp_scores_df %>% 
  dplyr::select(species_scientific_name, species, month, starts_with("PC")) %>% 
  pivot_wider(id_cols = c(species_scientific_name, species),
              names_from = "month",
              values_from = starts_with("PC"),
              names_sep = "_")
#Saving output...
write.csv(spp_month_pc_df,
          "./output/species_month_pc_scores.csv",
          row.names = FALSE)

#Saving the PC loadings to retain variable correlations with the PC axes
clr_diet_pca$rotation %>% 
  as.data.frame() %>% 
  mutate_all(round, digits = 4) %>% 
  mutate(variable = row.names(.)) %>% 
  write.csv("./output/pca_variable_corr.csv",
            row.names = FALSE)

