# Multiple Regression on Distance Matrices
# Relationship between topography and macroclimate

# https://github.com/csdambros/BioGeoAmazonia

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

species.list = list.files("./Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  macro.corr = numeric(),
  geo.temp.corr = numeric(),
  geo.ppt.corr = numeric(),
  Intercept = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  Intercept.temp = numeric(),
  Slope.temp = numeric(),
  Slope.p.value.temp = numeric(),
  R2.temp = numeric(),
  F.value.temp = numeric(),
  F.test.p.value.temp = numeric(),
  Intercept.ppt = numeric(),
  Slope.ppt = numeric(),
  Slope.p.value.ppt = numeric(),
  R2.ppt = numeric(),
  F.value.ppt = numeric(),
  F.test.p.value.ppt = numeric(),
  Intercept.geo = numeric(),
  Slope.geo = numeric(),
  Slope.p.value.geo = numeric(),
  R2.geo = numeric(),
  F.value.geo = numeric(),
  F.test.p.value.geo = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # creating spatial matrix
  macroclim.dat = species %>%
    select(macro_bio1_mean_annual_temp_C,macro_bio12_total_annual_precip_mm,) %>%
    scale()
  
  macro.temp = macroclim.dat[,1]
  macro.ppt = macroclim.dat[,2]
  
  macro.corr = cor(macro.temp,macro.ppt)
  
  # creating microclim data, scaled
  topo.dat = species %>%
    select(slope,aspect,mTPI) %>%
    scale()
  
  # calculate gower distance for scaled microclimate data
  topo.dist = gowdis(as.matrix(topo.dat))
  
  # calculate gower distance for macroclimate data
  macroclim.dist = gowdis(as.matrix(macroclim.dat))
  macro.temp.dist = gowdis(as.matrix(macro.temp))
  macro.ppt.dist = gowdis(as.matrix(macro.ppt))
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  
  geo.temp.corr = cor(geo.dist.2,macro.temp.dist)
  geo.ppt.corr = cor(geo.dist.2,macro.ppt.dist)
  
  # Perform MRM
  MRM.mod = MRM(topo.dist ~ macroclim.dist)
  
  # Extract from model
  Intercept = MRM.mod$coef[1,1]
  Slope = MRM.mod$coef[2,1]
  Slope.p.value = MRM.mod$coef[2,2]
  R2 = MRM.mod$r.squared[1]
  F.value = MRM.mod$F.test[1]
  F.test.p.value = MRM.mod$F.test[2]
  
  # Perform MRM with temp alone
  MRM.mod.temp = MRM(topo.dist ~ macro.temp.dist)
  
  # Extract from model
  Intercept.temp = MRM.mod.temp$coef[1,1]
  Slope.temp = MRM.mod.temp$coef[2,1]
  Slope.p.value.temp = MRM.mod.temp$coef[2,2]
  R2.temp = MRM.mod.temp$r.squared[1]
  F.value.temp = MRM.mod.temp$F.test[1]
  F.test.p.value.temp = MRM.mod.temp$F.test[2]
  
  # Perform MRM with ppt alone
  MRM.mod.ppt = MRM(topo.dist ~ macro.ppt.dist)
  
  # Extract from model
  Intercept.ppt = MRM.mod.ppt$coef[1,1]
  Slope.ppt = MRM.mod.ppt$coef[2,1]
  Slope.p.value.ppt = MRM.mod.ppt$coef[2,2]
  R2.ppt = MRM.mod.ppt$r.squared[1]
  F.value.ppt = MRM.mod.ppt$F.test[1]
  F.test.p.value.ppt = MRM.mod.ppt$F.test[2]
  
  # Perform MRM with geo dist
  MRM.mod.geo = MRM(topo.dist ~ geo.dist.2)
  
  # Extract from model
  Intercept.geo = MRM.mod.geo$coef[1,1]
  Slope.geo = MRM.mod.geo$coef[2,1]
  Slope.p.value.geo = MRM.mod.geo$coef[2,2]
  R2.geo = MRM.mod.geo$r.squared[1]
  F.value.geo = MRM.mod.geo$F.test[1]
  F.test.p.value.geo = MRM.mod.geo$F.test[2]
  
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
    Intercept.temp = Intercept.temp,
    Slope.temp = Slope.temp,
    Slope.p.value.temp = Slope.p.value.temp,
    R2.temp = R2.temp,
    F.value.temp = F.value.temp,
    F.test.p.value.temp = F.test.p.value.temp,
    Intercept.ppt = Intercept.ppt,
    Slope.ppt = Slope.ppt,
    Slope.p.value.ppt = Slope.p.value.ppt,
    R2.ppt = R2.ppt,
    F.value.ppt = F.value.ppt,
    F.test.p.value.ppt = F.test.p.value.ppt,
    Intercept.geo = Intercept.geo,
    Slope.geo = Slope.geo,
    Slope.p.value.geo = Slope.p.value.geo,
    R2.geo = R2.geo,
    F.value.geo = F.value.geo,
    F.test.p.value.geo = F.test.p.value.geo,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}

write.csv(results, file = "./Results/Q3.MRM.results.csv")
