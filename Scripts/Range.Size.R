install.packages("divvy")
library(divvy)

gbif = read.csv("final.gbif.data.csv", row.names = 1)

species_list <- split(gbif, gbif$species)

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