install.packages("divvy")
library(divvy)
library(tidyverse)

gbif = read.csv("./Formatted.Data/gbif.final.csv", row.names = 1)

# count and sort by number of occurrences
gbif.sorted = gbif %>%
  count(species) %>%                 
  arrange(n) 

gbif.2 <- gbif %>%
  mutate(name = factor(species, levels = gbif.sorted$species)) %>%
  arrange(species)

# make the species list
species_list <- split(gbif.2, gbif.2$name)

for(species in names(species_list)) {

  species.sub = species_list[[species]]
  coords = species.sub[,c(3,2)]
  crs = "epsg:4267"
  
  coords <- unique(coords)
  if (!isa(coords, "data.frame")) {
  coords <- as.data.frame(coords)
  }
  
  n <- nrow(coords)
  print(n)
  sfPts <- sf::st_as_sf(coords, coords = 1:2, crs = crs)
  if (!sf::st_is_longlat(sfPts)) {
  coords <- sf::sf_project(from = crs, to = "epsg:4326", 
                           coords, keep = TRUE, warn = TRUE)
  }
  
  latDiff <- max(coords[, 2]) - min(coords[, 2])
  latRng <- abs(latDiff)
  
  ptsGrp <- sf::st_union(sfPts)
  cntr <- unlist(sf::st_centroid(ptsGrp))
  centroidX = cntr[1]
  print(centroidX)
  centroidY = cntr[2]
  print(centroidY)
  print(latRng)

  gcdists <- sf::st_distance(sfPts)
  gcdists <- units::set_units(gcdists, "km")
  gcMax <- max(gcdists)
  print(gcMax)

  du = units::drop_units(gcdists)
  mst <- vegan::spantree(du)

  agg <- sum(mst$dist)
  diag(gcdists) <- NA
  mpd <- mean(gcdists, na.rm = TRUE)
  print(mpd)
  print(agg)
}



# count and sort by number of occurrences
gbif.sorted = gbif.2 %>%
  count(species) %>%                 
  arrange(n) 

gbif.3 <- gbif.2 %>%
  mutate(name = factor(species, levels = gbif.sorted$species)) %>%
  arrange(species)

# make the species list
species_list <- split(gbif.3, gbif.3$name)

results <- list()

for(species in names(species_list)) {
  # Create a file name for each species
  species.sub = species_list[[species]]
  species.occ = species.sub[,c(3,2)]
  range.size = as.data.frame(rangeSize(species.occ, crs = "epsg:4267"))
  range.size$species = unique(species.sub$species)
  results[[species]] = range.size
  
}

species.sub = species_list[[1]]

final.df = do.call(rbind, results)

range.size = rangeSize(species.occ, crs = "epsg:4267")



species.list.2 = species.list[30:122]
