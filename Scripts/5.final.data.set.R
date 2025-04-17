# Removing occurrence data for species that have more distribution in Canada than US

library(tidyverse)

# Read in current data
gbif = read.csv("./Formatted.Data/final.gbif.data.csv", row.names = 1)

# Get list of species with bad ranking

species = read.csv("./Raw.Data/species.csv")

# filter for BAD species - more of distribution in Canada than US

bad.species = species %>%
  filter(US.extent == "BAD")

# Remove occurence data of these species from gbif

gbif.2 = gbif %>%
  filter(!species %in% c(bad.species$Species))

sort(unique(gbif$species)) # 141 species
sort(unique(gbif.2$species)) # 123 species

write.csv(gbif.2, file = "./Formatted.Data/good.gbif.data.csv")
