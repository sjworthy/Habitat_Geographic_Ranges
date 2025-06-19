# Multiple Regression on Distance Matrices
# relationship between microclimate and geogrpahic distance
# code written to run on UNL HCC SWAN

# https://github.com/csdambros/BioGeoAmazonia

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

species.list = list.files("./Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  Intercept = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # creating spatial matrix
  spat.dat = as.matrix(species[,c("decimalLongitude","decimalLatitude")])
  
  # creating microclim data, scaled
  microclim.dat = species %>%
    dplyr::select(high_temp_C,low_temp_C,moisture_mm) %>%
    scale()
  
  # calculate gower distance for scaled microclimate data
  microclim.dist = gowdis(as.matrix(microclim.dat))
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/1000 # convert from m to km
  
  # Perform MRM
  MRM.mod = MRM(microclim.dist ~ geo.dist.3)
  
  # Extract from model
  Intercept = MRM.mod$coef[1,1]
  Slope = MRM.mod$coef[2,1]
  Slope.p.value = MRM.mod$coef[2,2]
  R2 = MRM.mod$r.squared[1]
  F.value = MRM.mod$F.test[1]
  F.test.p.value = MRM.mod$F.test[2]
  
  # Create a row with all the results for the current species
  species_name <- gsub("^clim_data_|\\.csv$", "", basename(species_file)) %>%
    gsub(" ", "_", .)  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    n = number.individuals,
    Intercept = Intercept,
    Slope = Slope,
    Slope.p.value = Slope.p.value,
    R2 = R2,
    F.value = F.value,
    F.test.p.value = F.test.p.value,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
  
  # save the model
  saveRDS(MRM.mod, file = paste0("./Microclimate/MRM_microclim_model_", species_name, ".RDS"))
  
  
  # Plot for each species
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    SoiDist = as.vector(microclim.dist)
  )
  
  plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoiDist)) +
    geom_point(shape = 21, fill = "grey", color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggtitle(paste("MRM Microclimate Plot", species.name)) +
    xlab("Geographic Distance (km)") +
    ylab("Microclimate Distance") +
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"))
  
  # Save the plot
  ggsave(filename = paste0("./Microclimate/MRM_Microclim_plot_", species.name, ".pdf"),
         plot = plot_obj,
         width = 5, height = 5)
}

write.csv(results, file = "./Microclimate/microclimate.MRM.results.csv")