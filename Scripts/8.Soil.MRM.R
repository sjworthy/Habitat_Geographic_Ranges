# Multiple Regression on Distance Matrices
# relationship between soil variables and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# read in soil data
soil.data = read.csv("gbif.data.soils.0.100cm.csv")
soil.data = read.csv("./Formatted.Data/gbif.data.soils.0.100cm.csv")

# remove rows with NAs
soil.data.2 = soil.data %>%
  drop_na()

# split into species
species.list = split(soil.data.2, soil.data.2$species)

soil.results <- data.frame(
  species = character(),
  n = numeric(),
  Intercept = numeric(),
  Intercept.p.value = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  R2.p.value = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  stringsAsFactors = FALSE
)

for(species.name in names(species.list)){
  
  species.data <- species.list[[species.name]]
  number.individuals = nrow(species.data)
  
  # creating spatial matrix
  spat.dat = as.matrix(species.data[,c("longitude","latitude")])
  
  # creating soil data
  soil.dat = species.data %>%
    dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)
  
  # calculate gower distance for scaled microclimate data
  soil.dist = gowdis(soil.dat)
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/1000 # convert from m to km
  
  # Perform MRM
  MRM.soil = MRM(soil.dist ~ geo.dist.3)
  
  # Extract from model
  Intercept = MRM.soil$coef[1,1]
  Intercept.p.value = MRM.soil$coef[1,2]
  Slope = MRM.soil$coef[2,1]
  Slope.p.value = MRM.soil$coef[2,2]
  R2 = MRM.soil$r.squared[1]
  R2.p.value = MRM.soil$r.squared[2]
  F.value = MRM.soil$F.test[[1]]
  F.test.p.value = MRM.soil$F.test[[2]]

  # Create a row with all the results for the current species
  clean_species_name = unique(species.data$species)

  species.results <- data.frame(
    species = clean_species_name,
    n = number.individuals,
    Intercept = Intercept,
    Intercept.p.value = Intercept.p.value,
    Slope = Slope,
    Slope.p.value = Slope.p.value,
    R2 = R2,
    R2.p.value = R2.p.value,
    F.value = F.value,
    F.test.p.value = F.test.p.value,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  soil.results <- rbind(soil.results, species.results)
  
  # save the model
  saveRDS(MRM.soil, file = paste0("./Soil/MRM_soil_model_", clean_species_name, ".RDS"))
  
  
  # Plot for each species
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    SoilDist = as.vector(soil.dist)
  )
  
    plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoilDist)) +
      geom_point(shape = 21, fill = "grey", color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggtitle(paste("MRM Soil Plot", clean_species_name)) +
      xlab("Geographic Distance (km)") +
      ylab("Soil Distance") +
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"))
    
    # Save the plot
    ggsave(filename = paste0("./Soil/MRM_Soil_plot_", clean_species_name, ".pdf"),
           plot = plot_obj,
           width = 5, height = 5)
  
}

write.csv(soil.results, file = "./Soil/soil.MRM.results.csv")