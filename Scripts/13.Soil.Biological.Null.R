# Code to generate the biologically plausible null
# WITHIN species: Calculate geographic distance and soil distance
# Randomly sample 50 points WITHIN each species - preserves soil and geographic distance within species
# Export the distance matrices as vectors 
# This null model does NOT break species-specific relationships between geographic and microclimate distance, 
# does NOT break biogeographic realms/distributions (i.e. southern distribution versus northern distribution), 
# accounts for dispersal limitation, but is constrained by species-specific adaptations to the environment.

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)
library(gghalves)

# Fraxinus albicans only has 49 data points.

# read in complete dataset

all.data = read.csv("./Formatted.Data/gbif.final.csv")

# split all.data into species-specific dataframe
split.data = split(all.data, all.data$species)

# randomly select 50 points for each species,
# calculate microclimate and geographic distance
# turn into vectors and save in a list

# Prepare a list to hold results from all 50 iterations
all_runs_data <- list()

# Repeat the process 50 times
for(run in 1:50){
  dist_values <- list()

for(i in 1:length(split.data)){
  sp = as.data.frame(split.data[[i]])
  sp.sample = slice_sample(sp, n = 50)

  # creating spatial matrix
  spat.dat = as.matrix(sp.sample[,c("decimalLongitude", "decimalLatitude")])
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/1000 # convert to km

  # creating soil data
  soil.dat = sp.sample %>%
    dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)
  
  # calculate gower distance for scaled microclimate data
  soil.dist = gowdis(soil.dat)
  
  # convert to vectors and combined into dataframe
  df_vecs <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    SoilDist = as.vector(soil.dist)
  )
  
  df_vecs$species <- unique(sp$species)

  # adding dataframe to the list
  dist_values[[length(dist_values) + 1]] <- df_vecs
}

# Combine values for all species

extracted_data <- do.call(rbind, dist_values)
all_runs_data[[run]] <- extracted_data
}

# Combine all 50 runs into one large data frame
combined_data <- do.call(rbind, all_runs_data)

# write.csv(extracted_data, file = "./Results/Soil.Biological.Null.dist.csv")

# save the model
# saveRDS(combined_data, file = "./Results/Soil.Biological.Null.50.RDS")

combined_data = readRDS("./Results/Soil.Biological.Null.50.RDS")

cor.test(combined_data$GeoDist, combined_data$SoilDist)
# R = 0.35


#### Plotting ####

SoilBioNull.plot = ggplot(extracted_data, aes(y = SoilDist, x = GeoDist)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Soil Biological Null", x = "Geographical Distance (km)", y = "Soil Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
SoilBioNull.plot

#ggsave(SoilBioNull.plot, file = "./Plots/Soil.Biological.Null.pdf", height = 5, width = 5)

Mean.SoilBioNull.plot = ggplot(combined_data, aes(y = SoilDist, x = GeoDist)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "50 Soil Biological Nulls", x = "Geographical Distance (km)", y = "Soil Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Mean.SoilBioNull.plot

ggsave(Mean.SoilBioNull.plot, file = "./Plots/Mean.Soil.Biological.Null.pdf", height = 5, width = 5)

#### Geographic Distance Bins ####

range(combined_data$GeoDist)
# 9.297789e-05 to 3.953471e+03

# Create quantile bins (5 bins) and convert to factor
combined_data$Geographic.Distance.Quantile <- cut(
  combined_data$GeoDist,
  breaks = quantile(combined_data$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(combined_data$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 9.297789e-05  1.556511e+02  3.647269e+02  6.600376e+02  1.087065e+03 3.953471e+03 

Bio.Null.Quant.Plot = ggplot(combined_data, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–155",
                              "2" = "156–364",
                              "3" = "365–660",
                              "4" = "661–1087",
                              "5" = "1088–3953")) +
  labs(y = "Soil Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Bio.Null.Quant.Plot

#ggsave(Bio.Null.Quant.Plot, file = "./Plots/Soil.Biol.Null.Quantile.pdf", height = 5, width = 5)


# bins for just 1 run

dat = read.csv("./Results/Soil.Biological.Null.dist.csv")

cor.test(dat$GeoDist,dat$SoilDist) # r = 0.35

range(dat$GeoDist)
# 1.113195e-04 3.774878e+03
# so range of 50 iterations is similar to range of 1 null

# Create quantile bins (5 bins) and convert to factor
dat$Geographic.Distance.Quantile <- cut(
  dat$GeoDist,
  breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 1.113195e-04  1.522914e+02  3.610014e+02  6.585560e+02  1.090959e+03 3.774878e+03  

Bio.1.Null.Quant.Plot = ggplot(dat, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–152",
                              "2" = "153–361",
                              "3" = "362–658",
                              "4" = "659–1090",
                              "5" = "1091–3774")) +
  labs(y = "Soil Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Bio.1.Null.Quant.Plot

#ggsave(Bio.1.Null.Quant.Plot, file = "./Plots/Soil.Biol.1.Null.Quantile.pdf", height = 5, width = 5)


 