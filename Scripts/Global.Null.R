# Code to generate the global null
# Randomly sample 50 points from each species
# Calculate geographic distance and microclimate distance – ignoring species identity but keeping location and microclimates of points together
# Run MRM – extract slopes, intercepts, R2 values
# Repeat 1000 times
# This null model breaks species-specific relationship between geographic and microclimate, 
# breaks biogeographic realms/distributions (i.e. southern distribution versus northern distribution), 
# makes the assumption of global dispersal, no ecological constraint except that a tree can grow here.

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# read in complete dataset

all.data = read.csv("./Results/climate.data.output.1.csv", row.names = 1)

# Initialize vector to store intercept and slope values
slopes <- numeric()
intercepts <- numeric()
R2 <- numeric()

for (i in 1:1000) {
# Sample 50 individuals from each species
sampled_indices <- all.data %>%
  group_by(species) %>%
  slice_sample(n = 50, replace = TRUE)

# creating spatial matrix
spat.dat = as.matrix(sampled_indices[,c("decimalLongitude", "decimalLatitude")])

# calculate Haversine distance for spatial data
geo.dist = distm(spat.dat, fun = distHaversine)
geo.dist.2 = as.dist(geo.dist) # convert to dist object
geo.dist.3 = geo.dist.2/1000 # convert to km

# creating microclimate data
microclim.dat = sampled_indices %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for microclimate data
microclim.dist = gowdis(as.matrix(microclim.dat))

# Run the MRM model with the bootstrapped matrices
model <- MRM(microclim.dist ~ geo.dist.3)

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
    MicroclimDist = as.vector(microclim.dist)
  )
  
  plot_obj = ggplot(df_plot, aes(x = GeoDist, y = MicroclimDist)) +
    geom_point(shape = 21, fill = "grey", color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggtitle(paste("MRM Plot - Iteration", i)) +
    xlab("Geographic Distance (km)") +
    ylab("Microclimate Distance") +
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"))
  
  # Save the plot
  ggsave(filename = paste0("MRM_Global_Null_plot_", i, ".pdf"),
         plot = plot_obj,
         width = 5, height = 5)
}

}

# combine intercept and slope vectors into a dataframe
null_output = as.data.frame(intercepts)
null_output$slopes = slopes
null_output$R2 = R2

write.csv(null_output, file = "global.null.results.csv")




