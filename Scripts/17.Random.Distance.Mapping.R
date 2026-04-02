# Generating microclimate, topography, and soil distance matrices
# Calculating mean distance from each point to all others for plotting
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# must sample points from each species. Can't compute the mean distance between a point 
# and all other points b/c need to convert to a matrix and that exhaust vector memory

mean_gower_sample = function(data, sample_size = 1000) {
  n = nrow(data)
  result = numeric(n)
  
  for (i in seq_len(n)) {
    # sample indices excluding self
    idx <- sample(setdiff(seq_len(n), i), 
                  size = min(sample_size, n - 1))
    
    # compute distances from point i to sampled points
    d <- gowdis(rbind(data[i, ], data[idx, ]))
    
    # distances from first row to others
    result[i] <- mean(as.numeric(d[1:length(idx)]), na.rm = TRUE)
  }
  
  return(result)
}

# read in data
all.data = read.csv("./Formatted.Data/All.Final.Data.csv", row.names = 1)

# split into species
species.list = split(all.data, all.data$species)

for(species.name in names(species.list)){
  
  species.data <- species.list[[species.name]]
  
  # randomly sample 48000 points (max vector length is 2.1, but have issues > 50,000)
  # this mean only Liquidambar styraciflua (n = 51647), Acer rubrum (n = 70995), 
  # and Quercus palustric (n = 90009) are randomly sampled
  
  if (nrow(species.data) > 48000) {
    set.seed(13)
    species.data <- species.data[sample(nrow(species.data), 48000), ]
  }
  
  # creating dataframe:
  microclimat.dat = species.data %>%
    dplyr::select(high_temp_C,low_temp_C,moisture_mm)
  
  topo.dat = species.data %>%
    dplyr::select(northness,eastness,mTPI,slope,elevation)
  
  soil.dat = species.data %>%
    dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)

  # compute sampled mean distances
  species.data$microclim.dist.mean = mean_gower_sample(microclimat.dat, 1000)
  species.data$topo.dist.mean = mean_gower_sample(topo.dat, 1000)
  species.data$soil.dist.mean = mean_gower_sample(soil.dat, 1000)
  
  # Create a row with all the results for the current species
  clean_species_name = unique(species.data$species)
  
  write.csv(species.data, file = paste0("./Mean_Dist/Mean.Dist.",clean_species_name,".csv"))
}
  