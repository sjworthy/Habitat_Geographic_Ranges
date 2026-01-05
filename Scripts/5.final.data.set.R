# Removing occurrence points for a variety of reasons

library(tidyverse)

#### Removing occurrence points for species that have more distribution in Canada than US ####

# Read in current data
gbif = read.csv("./Formatted.Data/final.gbif.data.csv", row.names = 1)

# Get list of species with bad ranking

species = read.csv("./Raw.Data/species.csv")

# filter for BAD species - more of distribution in Canada than US

bad.species = species %>%
  filter(US.extent == "BAD")

# Remove occurrence data of these species from gbif

gbif.2 = gbif %>%
  filter(!species %in% c(bad.species$Species))

sort(unique(gbif$species)) # 141 species
sort(unique(gbif.2$species)) # 123 species

#### remove Aesculus flava with only 8 individuals ####

sort(table(gbif.2$species)) 

gbif.3 = gbif.2 %>%
  filter(species != "Aesculus flava")

sort(unique(gbif.3$species)) # 122 species

# write.csv(gbif.3, file = "./Formatted.Data/good.gbif.data.csv")

#### Remove occurrence points without soil data ####
# read in previously filtered occurrence data
gbif = read.csv("./Formatted.Data/good.gbif.data.csv", row.names = 1)

# read in soil data
soil = read.csv("./Formatted.Data/gbif.data.soils.0.100cm.csv")
# 67889 are NA: ph_d0_100
# 66180 are NA: clay_d0_100
# 66773 are NA: sand_d0_100
# 67405 are NA: silt_d0_100
# 66547 are NA: db_d0_100
# 76944 are NA: ec_d0_100
# 67417 are NA: texture_d0_100

# drop the rows in soil with NA in any column
soil.final = soil %>%
  drop_na()
# total NA 82436

# looking at how many points for each species are missing
counts = as.data.frame(table(gbif$species))
colnames(counts)[1] = "species"
colnames(counts)[2] = "gbif"
counts$soil = table(soil.final$species)
counts$difference = counts$gbif - counts$soil
# write.csv(counts, "./Raw.Data/missing.soil.csv")

# filter gbif data for occurrence points with soil
soil.final = soil.final %>%
  rename(decimalLatitude = latitude, decimalLongitude = longitude)

gbif.soil = inner_join(gbif, soil.final)
# 1254429 occurrence points
# 122 species

#### Merge and filter for missing elevation ####

elev.final = read.csv("./Formatted.Data/elev.final.csv", row.names = 1)

gbif.elev = inner_join(gbif.soil,elev.final)
# 1254426, lose 3 individuals without elevation

# write.csv(gbif.elev, file = "Formatted.Data/gbif.final.csv")

#### Filter out final 45 data points without microclimate data ####

final.data.01.05.26 = read.csv("./Formatted.Data/final.data.01.05.26.csv", row.names = 1)
colSums(is.na(final.data.01.05.26))

all.final.data = na.omit(final.data.01.05.26)
colSums(is.na(all.final.data))

#write.csv(all.final.data, file = "./Formatted.Data/All.Final.Data.csv")