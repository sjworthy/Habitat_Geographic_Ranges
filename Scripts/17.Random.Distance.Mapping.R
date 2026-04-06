# Generating microclimate, topography, and soil distance matrices
# Calculating mean distance from each point to a random sampling of points
# Calculating min and max distance of each point to all other points

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

#### Min and Max Distance #####
# read in data
all.data = read.csv("./Formatted.Data/All.Final.Data.csv", row.names = 1)

# split into species
species.list = split(all.data, all.data$species)

all.results <- vector("list", length(species.list))
names(all.results) <- names(species.list)
                            
for(species.name in names(species.list)){
  
  species.data = species.list[[species.name]]
  
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

  # calculate gower distance
  microclim.dist = gowdis(microclimat.dat)
  topo.dist = gowdis(topo.dat)
  soil.dist = gowdis(soil.dat)
  
  n = nrow(species.data)
  
  min_dist_micro = numeric(n)
  max_dist_micro = numeric(n)
  min_dist_topo = numeric(n)
  max_dist_topo = numeric(n)
  min_dist_soil = numeric(n)
  max_dist_soil = numeric(n)
  
  for(i in seq_len(n)){
    microclim.di = microclim.dist[i,] # distances from species i to all others
    microclim.di.2 = microclim.di[microclim.di > 0] # drop self-distance (0)
    min_dist_micro[i] = min(microclim.di.2, na.rm = TRUE)
    max_dist_micro[i] = max(microclim.di.2, na.rm = TRUE)
    
    topo.di = topo.dist[i,] # distances from species i to all others
    topo.di.2 = topo.di[topo.di > 0] # drop self-distance (0)
    min_dist_topo[i] = min(topo.di.2, na.rm = TRUE)
    max_dist_topo[i] = max(topo.di.2, na.rm = TRUE)

    soil.di = soil.dist[i,] # distances from species i to all others
    soil.di.2 = soil.di[soil.di > 0] # drop self-distance (0)
    min_dist_soil[i] = min(soil.di.2, na.rm = TRUE)
    max_dist_soil[i] = max(soil.di.2, na.rm = TRUE)
  }
  
  all.results[[species.name]] <- data.frame(
    species = rep(species.name, n),
    min_microclim_distance = min_dist_micro,
    max_microclim_distance = max_dist_micro,
    min_topo_distance = min_dist_topo,
    max_topo_distance = max_dist_topo,
    min_soil_distance = min_dist_soil,
    max_soil_distance = max_dist_soil)
 
  safe_name <- gsub(" ", "_", species.name) 
  
  write.csv(
    all.results[[species.name]],
    file = paste0("./Distance.Results/", safe.name, "_gower_distances.csv"),
    row.names = FALSE)
  
}







### Mean Distance #####
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
  