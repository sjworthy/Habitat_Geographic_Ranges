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


#### Geographic Distance Bins ####

# read in data for first plot
# just picked 500 because it has the largest range
dat = read.csv("./Results/Soil.Global.Null/soil.global.dist_ 500 .csv", row.names = 1)

range(dat$GeoDist)
# 0.002836702 to 346.7406

# Create quantile bins (5 bins) and convert to factor
dat$Geographic.Distance.Quantile <- cut(
  dat$GeoDist,
  breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#       0%       20%       40%       60%       80%      100% 
# 0.002836702  52.84461  89.90785  124.9091  163.797 346.7406

Global.Null.Quant.Plot = ggplot(dat, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–53",
                              "2" = "54–90",
                              "3" = "91–124",
                              "4" = "125–163",
                              "5" = "164–346")) +
  labs(y = "Soil Distance",
       x = "Geographic Distance Quantile (hm)", fill = " ") +
  theme_classic()
Global.Null.Quant.Plot

ggsave(Global.Null.Quant.Plot, file = "./Plots/Soil.Global.Null.Quantile.png", height = 5, width = 5)


# Plot the data regular as well

plot_obj = ggplot(dat, aes(x = GeoDist, y = SoilDist)) +
  geom_point(shape = 21, fill = "grey", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ylim(0,1)+
  xlab("Geographic Distance (hm)") + # should be hectometer
  ylab("Soil Distance") +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
plot_obj

# Save the plot
ggsave(plot_obj, file = "./Plots/Soil.Global.Null.png", height = 5, width = 5)




