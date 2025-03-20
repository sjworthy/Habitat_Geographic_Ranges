# Multiple Regression on Distance Matrices

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

species.list = list.files("./Results/Species.Climate/", pattern = "*.csv", full.names = TRUE)

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
  spat.dat = as.matrix(species[,c(2,3)])
  
  # creating microclim data, scaled
  microclim.dat = species %>%
    select(high_temp_C,low_temp_C,moisture_mm) %>%
    scale()
  
  # calculate gower distance for scaled microclimate data
  microclim.dist = gowdis(as.matrix(microclim.dat))
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  
  # Perform MRM
  MRM.mod = MRM(microclim.dist ~ geo.dist.2)
  
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
}



  









