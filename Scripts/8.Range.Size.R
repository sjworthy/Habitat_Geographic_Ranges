install.packages("divvy")
library(divvy)

gbif = read.csv("./Formatted.Data/final.gbif.data.csv", row.names = 1)
gbif.2 = gbif %>%
  dplyr::filter(species %in% c("Acer nigrum",
                         "Acer rubrum",
                         "Fagus grandifolia",
                         "Fraxinus nigra",
                         "Larix laricina",
                         "Liquidambar styraciflua",
                         "Picea glauca",
                         "Picea mariana",
                         "Pinus banksiana",
                         "Pinus resinosa",
                         "Populus balsamifera",
                         "Populus deltoides",
                         "Quercus palustris",
                         "Quercus prinus",
                         "Taxodium ascendens",
                         "Thuja occidentalis"))

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