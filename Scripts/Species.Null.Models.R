# Species-level Null Models

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# each species data
species_files <- list.files("./Species.Nulls/", pattern = "*.csv", full.names = TRUE)

for(i in 1:length(species_files)){
# read in the species data
sp = read.csv(species_files[i])

# remove the species the null model is built for
sp.name = unique(sp$species)

sp.2 = sp %>%
  filter(!real.species %in% sp.name)

# randomly sample 45000 points b/c that gives ~ 1 billion data points (max is 2.1, but have issues > 50,000)
if (nrow(sp.2) > 45000) {
  set.seed(13)
  sp.2 <- sp.2[sample(nrow(sp.2), 45000), ]
}

# creating spatial matrix
spat.dat = as.matrix(sp.2[,c("decimalLongitude","decimalLatitude")])

# creating soil data
soil.dat = sp.2 %>%
  dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)

# calculate gower distance for scaled microclimate data
soil.dist = gowdis(soil.dat)

# calculate Haversine distance for spatial data
geo.dist = distm(spat.dat, fun = distHaversine)
geo.dist.2 = as.dist(geo.dist) # convert to dist object
geo.dist.3 = geo.dist.2/10000 # convert from m to hectometer

# Perform MRM
MRM.soil = MRM(soil.dist ~ geo.dist.3)

# save the model
saveRDS(MRM.soil, file = paste0("./Species.Null.Results/MRM_soil_species_null_", sp.name, ".RDS"))

# Plot for each species
df_plot <- data.frame(
  GeoDist = as.vector(geo.dist.3),
  SoilDist = as.vector(soil.dist)
)

write.csv(df_plot, file = paste0("./Species.Null.Results/MRM_Soil_species_null_data_", sp.name, ".csv"))

plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoilDist)) +
  geom_point(shape = 21, fill = "grey", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ggtitle(paste("MRM Soil Species Null", sp.name)) +
  xlab("Geographic Distance (hectometer)") +
  ylab("Soil Distance") +
  ylim(0,1)+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))

# Save the plot
ggsave(filename = paste0("./Species.Null.Results/MRM_Soil_species_null_", sp.name, ".png"),
       plot = plot_obj,
       width = 5, height = 5)

}
