# Code to generate the global null
# Randomly sample 5 points from each species
# Calculate geographic distance and soil distance
# Repeat 999 times
# Run MRM
# plot distributions of intercept, slope, R2
# plot MRM every 100th interaction
# This null model breaks species-specific relationship between geographic and microclimate, 
# breaks biogeographic realms/distributions (i.e. southern distribution versus northern distribution), 
# makes the assumption of global dispersal, no ecological constraint except that a tree can grow here.

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# read in complete dataset

all.data = read.csv("./Formatted.Data/gbif.final.csv", row.names = 1)

# Initialize vector to store intercept and slope values
slopes <- numeric()
intercepts <- numeric()
R2 <- numeric()

for(i in 1:999){
# Sample 5 individuals from each species
sampled_indices <- all.data %>%
  group_by(species) %>%
  slice_sample(n = 5) %>%
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
  
# calculate gower distance for scaled microclimate data
soil.dist = gowdis(as.data.frame(soil.dat))

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
    xlab("Geographic Distance (km)") +
    ylab("Soil Distance") +
    theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
          axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
          axis.title= element_text(face = "bold", size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black"))
plot_obj
    
# Save the plot
ggsave(filename = paste0("./Results/Global.Null/Soil_Global_Null_plot_", i, ".png"),
       plot = plot_obj,
       width = 5, height = 5)

write.csv(df_plot, file = paste("./Results/Global.Null/soil.global.dist_", i, ".csv"))
 }
}

# combine intercept and slope vectors into a dataframe
null_output = as.data.frame(intercepts)
null_output$slopes = slopes
null_output$R2 = R2

mean(null_output$intercepts)
# 0.2418721
mean(null_output$slopes)
# 0.0004634458
mean(null_output$R2)
# 0.05967644

ggplot(null_output, aes(intercepts))+
  geom_density()
ggplot(null_output, aes(slopes))+
  geom_density()
ggplot(null_output, aes(R2))+
  geom_density()

write.csv(null_output, file = "./Results/global.null.999.results.csv")


#### Geographic Distance Bins ####

range(df_plot$GeoDist)
# 0 to 3538.93

# Create quantile bins (5 bins) and convert to factor
df_plot$Geographic.Distance.Quantile <- cut(
  df_plot$GeoDist,
  breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0  524.6038  908.4138  1273.9590  1646.7984 3538.9299

Global.Null.Quant.Plot = ggplot(df_plot, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–524",
                              "2" = "525–908",
                              "3" = "909–1273",
                              "4" = "1274–1646",
                              "5" = "1647–3538")) +
  labs(y = "Soil Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Global.Null.Quant.Plot

ggsave(Global.Null.Quant.Plot, file = "./Plots/Soil.Global.Null.Quantile.pdf", height = 5, width = 5)






