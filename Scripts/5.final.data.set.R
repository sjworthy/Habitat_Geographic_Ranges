# Removing occurrence points for a variety of reasons

library(tidyverse)
library(rgbif)

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

gbif.2 = gbif %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

# read in soil data
soil = read.csv("./Formatted.Data/gbif.data.soils.0.100cm.csv")

# filter out introduced species
soil.2 = soil %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

# 62183 are NA: ph_d0_100
# 60521 are NA: clay_d0_100
# 61087 are NA: sand_d0_100
# 61690 are NA: silt_d0_100
# 60873 are NA: db_d0_100
# 71454 are NA: ec_d0_100
# 61700 are NA: texture_d0_100

# drop the rows in soil with NA in any column
soil.final = soil.2 %>%
  drop_na()
# total NA 51714

# looking at how many points for each species are missing
counts = as.data.frame(table(gbif.2$species))
colnames(counts)[1] = "species"
colnames(counts)[2] = "gbif"
counts$soil = table(soil.final$species)
counts$difference = counts$gbif - counts$soil
# write.csv(counts, "./Raw.Data/missing.soil.csv")

# filter gbif data for occurrence points with soil
soil.final = soil.final %>%
  rename(decimalLatitude = latitude, decimalLongitude = longitude)

gbif.soil = inner_join(gbif.2, soil.final)
# 1208687 occurrence points
# 119 species

#### Merge and filter for missing elevation ####

elev.final = read.csv("./Formatted.Data/elev.final.csv", row.names = 1)

gbif.elev = inner_join(gbif.soil,elev.final)
# 1208684, lose 3 individuals without elevation

# write.csv(gbif.elev, file = "Formatted.Data/gbif.final.csv")

#### Filter out final 45 data points without microclimate data ####

final.data.01.05.26 = read.csv("./Formatted.Data/final.data.01.05.26.csv", row.names = 1)

# filter introduced species
final.data.01.05.26.2 = final.data.01.05.26 %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

colSums(is.na(final.data.01.05.26.2))

all.final.data = na.omit(final.data.01.05.26.2)
colSums(is.na(all.final.data))

#write.csv(all.final.data, file = "./Formatted.Data/All.Final.Data.csv")

#### Filter out 3 introducted species ####

all.final.data = read.csv("./Formatted.Data/All.Final.Data.csv")
all.final.data.2 = all.final.data %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

write.csv(all.final.data.2, file = "./Formatted.Data/All.Final.Data.csv")

#### Summary Stats ####

# trying to get the original gbif columns to match with final data set
all.final.data = read.csv("./Formatted.Data/All.Final.Data.csv")
d = occ_download_get('0010272-241126133413365', path = "/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges/Raw.Data") %>%
  occ_download_import()

# merge the two 
merge.df = left_join(all.final.data,d, by = join_by("species","decimalLatitude","decimalLongitude"))
# final data is 1254381 rows
# merge data is 1373679 rows
# difference is 119298 rows

write.csv(merge.df, file = "./Formatted.Data/Final.Data.GBIF.csv")

# can't distinguish between the duplicates based on Lat./Long.

# average coordinate uncertainy
mean(merge.df$coordinateUncertaintyInMeters, na.rm = TRUE)

# mean number of individuals per species
table(all.final.data$species)
mean(table(all.final.data$species)) # 10,156.73
median(table(all.final.data$species)) # 4544

### Removing introduced species from gbif, post hoc ####
gbif.4 = gbif.3 %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))



