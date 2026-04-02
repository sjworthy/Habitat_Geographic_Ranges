# Code to generate the global null for soil
# Randomly sample 40 points from each species
# Calculate geographic distance and soil distance
# Repeat 999 times
# Run MRM
# plot distributions of intercept, slope, R2
# # Plot every 100th iteration

# 40 samples per species was chosen so that the global species would have 4,880 individuals,
# a similar value to the median number of individuals per species in the data (n = 4710).

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)
library(gghalves)

# read in complete dataset
all.data = read.csv("All.Final.Data.csv", row.names = 1)

# Initialize vector to store intercept and slope values
slopes <- numeric()
intercepts <- numeric()
R2 <- numeric()

for(i in 1:999){
# Sample 40 individuals from each species
sampled_indices <- all.data %>%
  group_by(species) %>%
  slice_sample(n = 40) %>%
  ungroup()
  
# creating spatial matrix
spat.dat = as.matrix(sampled_indices[,c("decimalLongitude", "decimalLatitude")])
  
# calculate Haversine distance for spatial data
geo.dist = distm(spat.dat, fun = distHaversine)
geo.dist.2 = as.dist(geo.dist) # convert to dist object
geo.dist.3 = geo.dist.2/10000 # convert to hectometer or 1/10th km
# 0.5 hectometers = 5000 meters and 5 km
  
# creating soil data
soil.dat = sampled_indices %>%
  dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)
  
# calculate gower distance for scaled soil data
soil.dist = gowdis(soil.dat)

# Run the MRM model with the bootstrapped matrices
  model <- MRM(soil.dist ~ geo.dist.3)

# Extract the intercept and slope values
  intercept_value <- model$coef[1,1]
  slope_value <- model$coef[2,1]
  R2_value <- model$r.squared[1]

 intercepts[i] <- intercept_value
 slopes[i] <- slope_value
 R2[i] <- R2_value

# Plot every 100th iteration
 if (i %% 100 == 0) {
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    SoilDist = as.vector(soil.dist))
    
  plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoilDist)) +
    geom_point(shape = 21, fill = "grey", color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ylim(0,1)+
    ggtitle(paste("Iteration", i)) +
    xlab("Geographic Distance (hectometer)") +
    ylab("Soil Distance") +
    theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
          axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
          axis.title= element_text(face = "bold", size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black"))
plot_obj
    
# Save the plot
ggsave(filename = paste0("./Soil.Global.Null/Soil_Global_Null_plot_", i, ".png"),
       plot = plot_obj,
       width = 5, height = 5)

write.csv(df_plot, file = paste("./Soil.Global.Null/soil.global.dist_", i, ".csv"))
 }
}

# combine intercept and slope vectors into a dataframe
null_output = as.data.frame(intercepts)
null_output$slopes = slopes
null_output$R2 = R2

write.csv(null_output, file = "./Soil.Global.Null/global.null.999.results.csv")
