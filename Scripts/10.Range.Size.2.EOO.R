library(divvy)
library(tidyverse)

gbif = read.csv("gbif.final.csv", row.names = 1)

# count and sort by number of occurrences
gbif.sorted = gbif %>%
  count(species) %>%                 
  arrange(n) 

gbif.2 <- gbif %>%
  mutate(name = factor(species, levels = gbif.sorted$species)) %>%
  arrange(species)

# make the species list
species_list <- split(gbif.2, gbif.2$name)

results <- list()

for(species in names(species_list)) {
  # Create a file name for each species
  species.sub = species_list[[species]]
  species.occ = species.sub[,c(3,2)]
  range.size = as.data.frame(rangeSize(species.occ, crs = "epsg:4267"))
  range.size$species = unique(species.sub$species)
  results[[species]] = range.size
}

final.df = do.call(rbind, results)

write.csv(final.df, file = "multiple.range.size.csv")
