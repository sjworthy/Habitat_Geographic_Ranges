# Code to generate the species-specific nulls
# Subset all GBIF occurrence data for points that fall inside the range of the species
# Remove occurrence points of the particular species we are working with
# # Calculate geographic distance and microclimate distance – ignoring species identity but keeping location and microclimates of points together
# Run MRM – extract slopes, intercepts, R2 values
# Export null distance matrices for later plotting on top of observed matrices
# This model allows us to address three issues:
# 1) Similar microclimate doesn't exist at certain geographic distances, when we see
# steep slopes between microclimate and geographic distance (C.tex example)
# 2) Similar microclimate exists, but species isn't finding it, doesn't care it exists,
# or it's occupied by another specie. Less overlap between available microclimate and what
# microclimate the species actually occupies
# 3) Population differences in physiological tolerances to climatic variation. 
# Local adaptation, more overlap between available microclimate and what microclimate species
# actually occupy.

library(sf)
library(maps)
library(tigris)
options(tigris_use_cache = TRUE)
library(tidyverse)
library(CoordinateCleaner)
library(BIEN)

# read in complete dataset, using soil data to get final list
all.data = read.csv("./Formatted.Data/gbif.final.csv", row.names = 1)
# make another column with the real species names
all.data$real.species = all.data$species

### Abies fraseri ####
# read in the range map of the first species
ABFR.range = st_read("../USTreeAtlas/shp/abiefras/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ABFR.range)  <- 4267 

# get US base map, only done for the first species
us_states <- states(cb = TRUE)
continental_states <- us_states %>%
  filter(!NAME %in% (c("Alaska","American Samoa","Guam","Commonwealth of the Northern Mariana Islands","Hawaii","United States Virgin Islands",
                       "Puerto Rico")))
states.map = continental_states %>%
  st_as_sf %>%
  st_transform(st_crs(ABFR.range))

# clip the species range based on US map
ABFR_clipped = st_intersection(ABFR.range, states.map)

# Add the species name back to the map info, also make all points same species so
# cc_iucn will work and think all points are same species for the range
all.data$species = "Abies fraseri"
ABFR_clipped$species = "Abies fraseri"

# clip all points outside of range with 50000 m buffer
ABFR_flag = cc_iucn(x = all.data, range = ABFR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ABFR_occ_final = all.data[ABFR_flag, ]

# plot to make sure things are working correctly
ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = ABFR_clipped, col = "red", linewidth = 2)+
  geom_point(data = ABFR_occ_final, aes(x = decimalLongitude, y = decimalLatitude), color = "black", alpha = 0.7)+
  theme_classic()

# export the values
write.csv(ABFR_occ_final, file = "./Formatted.Data/Species.Nulls/Abies.fraseri.null.csv")

#### Acer barbatum ####
# read in the range map of the first species
ACBA3.range = st_read("../USTreeAtlas/shp/acerbarb/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACBA3.range)  <- 4267 

# clip the species range based on US map
ACBA3_clipped = st_intersection(ACBA3.range, states.map)

# Add the species name back to the map info
all.data$species = "Acer barbatum"
ACBA3_clipped$species = "Acer barbatum"

# clip all points outside of range with 50000 m buffer
ACBA3_flag = cc_iucn(x = all.data, range = ACBA3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ACBA3_occ_final = all.data[ACBA3_flag, ]

# export the values
write.csv(ACBA3_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.barbatum.null.csv")

#### Acer leucoderme ####
# read in the range map of the first species
ACLE.range = st_read("../USTreeAtlas/shp/acerleuc/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACLE.range)  <- 4267 

# clip the species range based on US map
ACLE_clipped = st_intersection(ACLE.range, states.map)

# Add the species name back to the map info
ACLE_clipped$species = "Acer leucoderme"
all.data$species = "Acer leucoderme"

# clip all points outside of range with 50000 m buffer
ACLE_flag = cc_iucn(x = all.data, range = ACLE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ACLE_occ_final = all.data[ACLE_flag, ]

# export the values
write.csv(ACLE_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.leucoderme.null.csv")

#### Acer negundo ####
# read in the range map of the first species
ACNE3.range = st_read("../USTreeAtlas/shp/acernegu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACNE3.range)  <- 4267 

# clip the species range based on US map
ACNE3_clipped = st_intersection(ACNE3.range, states.map)

# Add the species name back to the map info
ACNE3_clipped$species = "Acer negundo"
all.data$species = "Acer negundo"

# clip all points outside of range with 50000 m buffer
ACNE3_flag = cc_iucn(x = all.data, range = ACNE3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ACNE3_occ_final = all.data[ACNE3_flag, ]

# export the values
write.csv(ACNE3_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.negundo.null.csv")

#### Acer nigrum ####
# read in the range map of the first species
ACNI5.range = st_read("../USTreeAtlas/shp/acernigr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACNI5.range)  <- 4267 

# clip the species range based on US map
ACNI5_clipped = st_intersection(ACNI5.range, states.map)

# Add the species name back to the map info
ACNI5_clipped$species = "Acer nigrum"
all.data$species = "Acer nigrum"

# clip all points outside of range with 50000 m buffer
ACNI5_flag = cc_iucn(x = all.data, range = ACNI5_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ACNI5_occ_final = all.data[ACNI5_flag, ]

# export the values
write.csv(ACNI5_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.nigrum.null.csv")

#### Acer pensylvanicum ####
# read in the range map of the first species
ACPE.range = st_read("../USTreeAtlas/shp/acerpens/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACPE.range)  <- 4267 

# clip the species range based on US map
ACPE_clipped = st_intersection(ACPE.range, states.map)

# Add the species name back to the map info
ACPE_clipped$species = "Acer pensylvanicum"
all.data$species = "Acer pensylvanicum"

# clip all points outside of range with 50000 m buffer
ACPE_flag = cc_iucn(x = all.data, range = ACPE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ACPE_occ_final = all.data[ACPE_flag, ]

# export the values
write.csv(ACPE_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.pensylvanicum.null.csv")

#### Acer rubrum ####
# read in the range map of the first species
ACRU.range = st_read("../USTreeAtlas/shp/acerrubr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACRU.range)  <- 4267 

# clip the species range based on US map
ACRU_clipped = st_intersection(ACRU.range, states.map)

# Add the species name back to the map info
ACRU_clipped$species = "Acer rubrum"
all.data$species = "Acer rubrum"

# clip all points outside of range with 50000 m buffer
ACRU_flag = cc_iucn(x = all.data, range = ACRU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ACRU_occ_final = all.data[ACRU_flag, ]

# export the values
write.csv(ACRU_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.rubrum.null.csv")

#### Acer saccharinum ####
# read in the range map of the first species
ACSA2.range = st_read("../USTreeAtlas/shp/acersacc/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACSA2.range)  <- 4267 

# clip the species range based on US map
ACSA2_clipped = st_intersection(ACSA2.range, states.map)

# Add the species name back to the map info
ACSA2_clipped$species = "Acer saccharinum"
all.data$species = "Acer saccharinum"

# clip all points outside of range with 50000 m buffer
ACSA2_flag = cc_iucn(x = all.data, range = ACSA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ACSA2_occ_final = all.data[ACSA2_flag, ]

# export the values
write.csv(ACSA2_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.Saccharinum.null.csv")

#### Acer saccharum ####
# read in the range map of the first species
ACSA3.range = st_read("../USTreeAtlas/shp/acersacr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ACSA3.range)  <- 4267 

# clip the species range based on US map
ACSA3_clipped = st_intersection(ACSA3.range, states.map)

# Add the species name back to the map info
ACSA3_clipped$species = "Acer saccharum"
all.data$species = "Acer saccharum"

# clip all points outside of range with 50000 m buffer
ACSA3_flag = cc_iucn(x = all.data, range = ACSA3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ACSA3_occ_final = all.data[ACSA3_flag, ]

# export the values
write.csv(ACSA3_occ_final, file = "./Formatted.Data/Species.Nulls/Acer.saccharum.null.csv")

#### Aesculus glabra ####
# read in the range map of the first species
AEGL.range = st_read("../USTreeAtlas/shp/aescglab/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(AEGL.range)  <- 4267 

# clip the species range based on US map
AEGL_clipped = st_intersection(AEGL.range, states.map)

# Add the species name back to the map info
AEGL_clipped$species = "Aesculus glabra"
all.data$species = "Aesculus glabra"

# clip all points outside of range with 50000 m buffer
AEGL_flag = cc_iucn(x = all.data, range = AEGL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
AEGL_occ_final = all.data[AEGL_flag, ]

# export the values
write.csv(AEGL_occ_final, file = "./Formatted.Data/Species.Nulls/Aesculus.glabra.null.csv")

#### Ailanthus altissima ####
# clip the species range based on US map, from BIEN
(AIAL.range.sf <- BIEN_ranges_load_species('Ailanthus altissima'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
AIAL.range.2 = AIAL.range.sf %>%
  st_transform(st_crs(states.map))

AIAL.range.3 = terra::vect(AIAL.range.2)
states.map.2 = terra::vect(states.map)
AIAL.range.4 = terra::intersect(AIAL.range.3, states.map.2)
AIAL.range.5 = st_as_sf(AIAL.range.4)

# Add the species name back to the map info
AIAL.range.5$species = "Ailanthus altissima"
all.data$species = "Ailanthus altissima"

# clip all points outside of range with 50000 m buffer
AIAL_flag = cc_iucn(x = all.data, range = AIAL.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
AIAL_occ_final = all.data[AIAL_flag, ]

# export the values
write.csv(AIAL_occ_final, file = "./Formatted.Data/Species.Nulls/Ailanthus.altissima.null.csv")

#### Annona glabra ####
# clip the species range based on US map, from BIEN
(ANGL4.range.sf <- BIEN_ranges_load_species('Annona glabra'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
ANGL4.range.2 = ANGL4.range.sf %>%
  st_transform(st_crs(states.map))

ANGL4.range.3 = terra::vect(ANGL4.range.2)
states.map.2 = terra::vect(states.map)
ANGL4.range.4 = terra::intersect(ANGL4.range.3, states.map.2)
ANGL4.range.5 = st_as_sf(ANGL4.range.4)

# Add the species name back to the map info
ANGL4.range.5$species = "Annona glabra"
all.data$species = "Annona glabra"

# clip all points outside of range with 50000 m buffer
ANGL4_flag = cc_iucn(x = all.data, range = ANGL4.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
ANGL4_occ_final = all.data[ANGL4_flag, ]

# export the values
write.csv(ANGL4_occ_final, file = "./Formatted.Data/Species.Nulls/Annona.glabra.null.csv")

#### Asimina triloba ####
# read in the range map of the first species
ASTR.range = st_read("../USTreeAtlas/shp/asimtril/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ASTR.range)  <- 4267 

# clip the species range based on US map
ASTR_clipped = st_intersection(ASTR.range, states.map)

# Add the species name back to the map info
ASTR_clipped$species = "Asimina triloba"
all.data$species = "Asimina triloba"

# clip all points outside of range with 50000 m buffer
ASTR_flag = cc_iucn(x = all.data, range = ASTR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ASTR_occ_final = all.data[ASTR_flag, ]

# export the values
write.csv(ASTR_occ_final, file = "./Formatted.Data/Species.Nulls/Asimina.triloba.null.csv")

#### Betula alleghaniensis ####
# read in the range map of the first species
BEAL2.range = st_read("../USTreeAtlas/shp/betualle/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(BEAL2.range)  <- 4267 

# clip the species range based on US map
BEAL2_clipped = st_intersection(BEAL2.range, states.map)

# Add the species name back to the map info
BEAL2_clipped$species = "Betula alleghaniensis"
all.data$species = "Betula alleghaniensis"

# clip all points outside of range with 50000 m buffer
BEAL2_flag = cc_iucn(x = all.data, range = BEAL2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
BEAL2_occ_final = all.data[BEAL2_flag, ]

# export the values
write.csv(BEAL2_occ_final, file = "./Formatted.Data/Species.Nulls/Betula.alleghaniensis.null.csv")

#### Betula lenta ####
# read in the range map of the first species
BELE.range = st_read("../USTreeAtlas/shp/betulent/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(BELE.range)  <- 4267 

# clip the species range based on US map
BELE_clipped = st_intersection(BELE.range, states.map)

# Add the species name back to the map info
BELE_clipped$species = "Betula lenta"
all.data$species = "Betula lenta"

# clip all points outside of range with 50000 m buffer
BELE_flag = cc_iucn(x = all.data, range = BELE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
BELE_occ_final = all.data[BELE_flag, ]

# export the values
write.csv(BELE_occ_final, file = "./Formatted.Data/Species.Nulls/Betula.lenta.null.csv")

#### Betula nigra ####
# read in the range map of the first species
BENI.range = st_read("../USTreeAtlas/shp/betunigr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(BENI.range)  <- 4267 

# clip the species range based on US map
BENI_clipped = st_intersection(BENI.range, states.map)

# Add the species name back to the map info
BENI_clipped$species = "Betula nigra"
all.data$species = "Betula nigra"

# clip all points outside of range with 50000 m buffer
BENI_flag = cc_iucn(x = all.data, range = BENI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
BENI_occ_final = all.data[BENI_flag, ]

# export the values
write.csv(BENI_occ_final, file = "./Formatted.Data/Species.Nulls/Betula.nigra.null.csv")

#### Betula populifolia ####
# read in the range map of the first species
BEPO.range = st_read("../USTreeAtlas/shp/betupopu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(BEPO.range)  <- 4267 

# clip the species range based on US map
BEPO_clipped = st_intersection(BEPO.range, states.map)

# Add the species name back to the map info
BEPO_clipped$species = "Betula populifolia"
all.data$species = "Betula populifolia"

# clip all points outside of range with 50000 m buffer
BEPO_flag = cc_iucn(x = all.data, range = BEPO_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
BEPO_occ_final = all.data[BEPO_flag, ]

# export the values
write.csv(BEPO_occ_final, file = "./Formatted.Data/Species.Nulls/Betula.populifolia.null.csv")

#### Carpinus caroliniana ####
# read in the range map of the first species
CACA18.range = st_read("../USTreeAtlas/shp/carpcaro/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CACA18.range)  <- 4267 

# clip the species range based on US map
CACA18_clipped = st_intersection(CACA18.range, states.map)

# Add the species name back to the map info
CACA18_clipped$species = "Carpinus caroliniana"
all.data$species = "Carpinus caroliniana"

# clip all points outside of range with 50000 m buffer
CACA18_flag = cc_iucn(x = all.data, range = CACA18_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                      value = "flagged", buffer = 50000)
# remove the points
CACA18_occ_final = all.data[CACA18_flag, ]

# export the values
write.csv(CACA18_occ_final, file = "./Formatted.Data/Species.Nulls/Carpinus.caroliniana.null.csv")

#### Carya alba ####
# clip the species range based on US map, from BIEN
(CAAL27.range.sf <- BIEN_ranges_load_species('Carya alba'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
CAAL27.range.2 = CAAL27.range.sf %>%
  st_transform(st_crs(states.map))

CAAL27.range.3 = terra::vect(CAAL27.range.2)
states.map.2 = terra::vect(states.map)
CAAL27.range.4 = terra::intersect(CAAL27.range.3, states.map.2)
CAAL27.range.5 = st_as_sf(CAAL27.range.4)

# Add the species name back to the map info
CAAL27.range.5$species = "Carya alba"
all.data$species = "Carya alba"

# clip all points outside of range with 50000 m buffer
CAAL27_flag = cc_iucn(x = all.data, range = CAAL27.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                      value = "flagged", buffer = 50000)
# remove the points
CAAL27_occ_final = all.data[CAAL27_flag, ]

# export the values
write.csv(CAAL27_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.alba.null.csv")

#### Carya aquatica ####
# read in the range map of the first species
CAAQ2.range = st_read("../USTreeAtlas/shp/caryaqua/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CAAQ2.range)  <- 4267 

# clip the species range based on US map
CAAQ2_clipped = st_intersection(CAAQ2.range, states.map)

# Add the species name back to the map info
CAAQ2_clipped$species = "Carya aquatica"
all.data$species = "Carya aquatica"

# clip all points outside of range with 50000 m buffer
CAAQ2_flag = cc_iucn(x = all.data, range = CAAQ2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CAAQ2_occ_final = all.data[CAAQ2_flag, ]

# export the values
write.csv(CAAQ2_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.aquatica.null.csv")

#### Carya cordiformis ####
# read in the range map of the first species
CACO15.range = st_read("../USTreeAtlas/shp/carycord/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CACO15.range)  <- 4267 

# clip the species range based on US map
CACO15_clipped = st_intersection(CACO15.range, states.map)

# Add the species name back to the map info
CACO15_clipped$species = "Carya cordiformis"
all.data$species = "Carya cordiformis"

# clip all points outside of range with 50000 m buffer
CACO15_flag = cc_iucn(x = all.data, range = CACO15_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                      value = "flagged", buffer = 50000)
# remove the points
CACO15_occ_final = all.data[CACO15_flag, ]

# export the values
write.csv(CACO15_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.cordiformis.null.csv")

#### Carya glabra ####
# read in the range map of the first species
CAGL8.range = st_read("../USTreeAtlas/shp/caryglab/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CAGL8.range)  <- 4267 

# clip the species range based on US map
CAGL8_clipped = st_intersection(CAGL8.range, states.map)

# Add the species name back to the map info
CAGL8_clipped$species = "Carya glabra"
all.data$species = "Carya glabra"

# clip all points outside of range with 50000 m buffer
CAGL8_flag = cc_iucn(x = all.data, range = CAGL8_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CAGL8_occ_final = all.data[CAGL8_flag, ]

# export the values
write.csv(CAGL8_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.glabra.null.csv")

#### Carya illinoinensis ####
# read in the range map of the first species
CAIL2.range = st_read("../USTreeAtlas/shp/caryilli/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CAIL2.range)  <- 4267 

# clip the species range based on US map
CAIL2_clipped = st_intersection(CAIL2.range, states.map)

# Add the species name back to the map info
CAIL2_clipped$species = "Carya illinoinensis"
all.data$species = "Carya illinoinensis"

# clip all points outside of range with 50000 m buffer
CAIL2_flag = cc_iucn(x = all.data, range = CAIL2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CAIL2_occ_final = all.data[CAIL2_flag, ]

# export the values
write.csv(CAIL2_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.illinoinensis.null.csv")

#### Carya ovata ####
# read in the range map of the first species
CAOV2.range = st_read("../USTreeAtlas/shp/caryovat/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CAOV2.range)  <- 4267 

# clip the species range based on US map
CAOV2_clipped = st_intersection(CAOV2.range, states.map)

# Add the species name back to the map info
CAOV2_clipped$species = "Carya ovata"
all.data$species = "Carya ovata"

# clip all points outside of range with 50000 m buffer
CAOV2_flag = cc_iucn(x = all.data, range = CAOV2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CAOV2_occ_final = all.data[CAOV2_flag, ]

# export the values
write.csv(CAOV2_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.ovata.null.csv")

#### Carya texana ####
# read in the range map of the first species
CATE9.range = st_read("../USTreeAtlas/shp/carytexa/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CATE9.range)  <- 4267 

# clip the species range based on US map
CATE9_clipped = st_intersection(CATE9.range, states.map)

# Add the species name back to the map info
CATE9_clipped$species = "Carya texana"
all.data$species = "Carya texana"

# clip all points outside of range with 50000 m buffer
CATE9_flag = cc_iucn(x = all.data, range = CATE9_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CATE9_occ_final = all.data[CATE9_flag, ]

# export the values
write.csv(CATE9_occ_final, file = "./Formatted.Data/Species.Nulls/Carya.texana.null.csv")

#### Castanea dentata ####
# read in the range map of the first species
CADE12.range = st_read("../USTreeAtlas/shp/castdent/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CADE12.range)  <- 4267 

# clip the species range based on US map
CADE12_clipped = st_intersection(CADE12.range, states.map)

# Add the species name back to the map info
CADE12_clipped$species = "Castanea dentata"
all.data$species = "Castanea dentata"

# clip all points outside of range with 50000 m buffer
CADE12_flag = cc_iucn(x = all.data, range = CADE12_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                      value = "flagged", buffer = 50000)
# remove the points
CADE12_occ_final = all.data[CADE12_flag, ]

# export the values
write.csv(CADE12_occ_final, file = "./Formatted.Data/Species.Nulls/Castanea.dentata.null.csv")

#### Celtis laevigata ####
# read in the range map of the first species
CELA.range = st_read("../USTreeAtlas/shp/celtlaev/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CELA.range)  <- 4267 

# clip the species range based on US map
CELA_clipped = st_intersection(CELA.range, states.map)

# Add the species name back to the map info
CELA_clipped$species = "Celtis laevigata"
all.data$species = "Celtis laevigata"

# clip all points outside of range with 50000 m buffer
CELA_flag = cc_iucn(x = all.data, range = CELA_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
CELA_occ_final = all.data[CELA_flag, ]

# export the values
write.csv(CELA_occ_final, file = "./Formatted.Data/Species.Nulls/Celtis.laevigata.null.csv")

#### Celtis occidentalis ####
# read in the range map of the first species
CEOC.range = st_read("../USTreeAtlas/shp/celtocci/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CEOC.range)  <- 4267 

# clip the species range based on US map
CEOC_clipped = st_intersection(CEOC.range, states.map)

# Add the species name back to the map info
CEOC_clipped$species = "Celtis occidentalis"
all.data$species = "Celtis occidentalis"

# clip all points outside of range with 50000 m buffer
CEOC_flag = cc_iucn(x = all.data, range = CEOC_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
CEOC_occ_final = all.data[CEOC_flag, ]

# export the values
write.csv(CEOC_occ_final, file = "./Formatted.Data/Species.Nulls/Celtis.occidentalis.null.csv")

#### Cercis canadensis ####
# clip the species range based on US map, from BIEN
(CECA4.range.sf <- BIEN_ranges_load_species('Cercis canadensis'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
CECA4.range.2 = CECA4.range.sf %>%
  st_transform(st_crs(states.map))

CECA4.range.3 = terra::vect(CECA4.range.2)
states.map.2 = terra::vect(states.map)
CECA4.range.4 = terra::intersect(CECA4.range.3, states.map.2)
CECA4.range.5 = st_as_sf(CECA4.range.4)

# Add the species name back to the map info
CECA4.range.5$species = "Cercis canadensis"
all.data$species = "Cercis canadensis"

# clip all points outside of range with 50000 m buffer
CECA4_flag = cc_iucn(x = all.data, range = CECA4.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CECA4_occ_final = all.data[CECA4_flag, ]

# export the values
write.csv(CECA4_occ_final, file = "./Formatted.Data/Species.Nulls/Cercis.canadensis.null.csv")

#### Chamaecyparis thyoides ####
# read in the range map of the first species
CHTH2.range = st_read("../USTreeAtlas/shp/chamthyo/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CHTH2.range)  <- 4267 

# clip the species range based on US map
CHTH2_clipped = st_intersection(CHTH2.range, states.map)

# Add the species name back to the map info
CHTH2_clipped$species = "Chamaecyparis thyoides"
all.data$species = "Chamaecyparis thyoides"

# clip all points outside of range with 50000 m buffer
CHTH2_flag = cc_iucn(x = all.data, range = CHTH2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CHTH2_occ_final = all.data[CHTH2_flag, ]

# export the values
write.csv(CHTH2_occ_final, file = "./Formatted.Data/Species.Nulls/Chamaecyparis.thyoides.null.csv")

#### Cornus florida ####
# read in the range map of the first species
COFL2.range = st_read("../USTreeAtlas/shp/cornflor/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(COFL2.range)  <- 4267 

# clip the species range based on US map
COFL2_clipped = st_intersection(COFL2.range, states.map)

# Add the species name back to the map info
COFL2_clipped$species = "Cornus florida"
all.data$species = "Cornus florida"

# clip all points outside of range with 50000 m buffer
COFL2_flag = cc_iucn(x = all.data, range = COFL2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
COFL2_occ_final = all.data[COFL2_flag, ]

# export the values
write.csv(COFL2_occ_final, file = "./Formatted.Data/Species.Nulls/Cornus.florida.null.csv")

#### Diospyros virginiana ####
# read in the range map of the first species
DIVI5.range = st_read("../USTreeAtlas/shp/diosvirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(DIVI5.range)  <- 4267 

# clip the species range based on US map
DIVI5_clipped = st_intersection(DIVI5.range, states.map)

# Add the species name back to the map info
DIVI5_clipped$species = "Diospyros virginiana"
all.data$species = "Diospyros virginiana"

# clip all points outside of range with 50000 m buffer
DIVI5_flag = cc_iucn(x = all.data, range = DIVI5_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
DIVI5_occ_final = all.data[DIVI5_flag, ]

# export the values
write.csv(DIVI5_occ_final, file = "./Formatted.Data/Species.Nulls/Diospyros.virginiana.null.csv")

#### Fagus grandifolia ####
# read in the range map of the first species
FAGR.range = st_read("../USTreeAtlas/shp/fagugran/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FAGR.range)  <- 4267 

# clip the species range based on US map
FAGR_clipped = st_intersection(FAGR.range, states.map)

# Add the species name back to the map info
FAGR_clipped$species = "Fagus grandifolia"
all.data$species = "Fagus grandifolia"

# clip all points outside of range with 50000 m buffer
FAGR_flag = cc_iucn(x = all.data, range = FAGR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
FAGR_occ_final = all.data[FAGR_flag, ]

# export the values
write.csv(FAGR_occ_final, file = "./Formatted.Data/Species.Nulls/Fagus.grandifolia.null.csv")

#### Fraxinus americana ####
# read in the range map of the first species
FRAM2.range = st_read("../USTreeAtlas/shp/fraxamer/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRAM2.range)  <- 4267 

# clip the species range based on US map
FRAM2_clipped = st_intersection(FRAM2.range, states.map)

# Add the species name back to the map info
FRAM2_clipped$species = "Fraxinus americana"
all.data$species = "Fraxinus americana"

# clip all points outside of range with 50000 m buffer
FRAM2_flag = cc_iucn(x = all.data, range = FRAM2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
FRAM2_occ_final = all.data[FRAM2_flag, ]

# export the values
write.csv(FRAM2_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.americana.null.csv")

#### Fraxinus caroliniana ####
# read in the range map of the first species
FRCA3.range = st_read("../USTreeAtlas/shp/fraxcaro/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRCA3.range)  <- 4267 

# clip the species range based on US map
FRCA3_clipped = st_intersection(FRCA3.range, states.map)

# Add the species name back to the map info
FRCA3_clipped$species = "Fraxinus caroliniana"
all.data$species = "Fraxinus caroliniana"

# clip all points outside of range with 50000 m buffer
FRCA3_flag = cc_iucn(x = all.data, range = FRCA3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
FRCA3_occ_final = all.data[FRCA3_flag, ]

# export the values
write.csv(FRCA3_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.caroliniana.null.csv")

#### Fraxinus pennsylvanica ####
# read in the range map of the first species
FRPE.range = st_read("../USTreeAtlas/shp/fraxpenn/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRPE.range)  <- 4267 

# clip the species range based on US map
FRPE_clipped = st_intersection(FRPE.range, states.map)

# Add the species name back to the map info
FRPE_clipped$species = "Fraxinus pennsylvanica"
all.data$species = "Fraxinus pennsylvanica"

# clip all points outside of range with 50000 m buffer
FRPE_flag = cc_iucn(x = all.data, range = FRPE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
FRPE_occ_final = all.data[FRPE_flag, ]

# export the values
write.csv(FRPE_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.pennsylvanica.null.csv")

#### Fraxinus profunda ####
# read in the range map of the first species
FRPR.range = st_read("../USTreeAtlas/shp/fraxprof/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRPR.range)  <- 4267 

# clip the species range based on US map
FRPR_clipped = st_intersection(FRPR.range, states.map)

# Add the species name back to the map info
FRPR_clipped$species = "Fraxinus profunda"
all.data$species = "Fraxinus profunda"

# clip all points outside of range with 50000 m buffer
FRPR_flag = cc_iucn(x = all.data, range = FRPR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
FRPR_occ_final = all.data[FRPR_flag, ]

# export the values
write.csv(FRPR_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.profunda.null.csv")

#### Fraxinus quadrangulata ####
# read in the range map of the first species
FRQU.range = st_read("../USTreeAtlas/shp/fraxquad/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRQU.range)  <- 4267 

# clip the species range based on US map
FRQU_clipped = st_intersection(FRQU.range, states.map)

# Add the species name back to the map info
FRQU_clipped$species = "Fraxinus quadrangulata"
all.data$species = "Fraxinus quadrangulata"

# clip all points outside of range with 50000 m buffer
FRQU_flag = cc_iucn(x = all.data, range = FRQU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
FRQU_occ_final = all.data[FRQU_flag, ]

# export the values
write.csv(FRQU_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.quadrangulata.null.csv")

#### Fraxinus albicans ####
# Synonym with Fraxinus texensis

# read in the range map of the first species
FRTE.range = st_read("../USTreeAtlas/shp/fraxtexe/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(FRTE.range)  <- 4267 

# clip the species range based on US map
FRTE_clipped = st_intersection(FRTE.range, states.map)

# Add the species name back to the map info
FRTE_clipped$species = "Fraxinus albicans"
all.data$species = "Fraxinus albicans"

# clip all points outside of range with 50000 m buffer
FRTE_flag = cc_iucn(x = all.data, range = FRTE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
FRTE_occ_final = all.data[FRTE_flag, ]

# export the values
write.csv(FRTE_occ_final, file = "./Formatted.Data/Species.Nulls/Fraxinus.albicans.null.csv")

#### Gleditsia aquatica ####
# read in the range map of the first species
GLAQ.range = st_read("../USTreeAtlas/shp/gledaqua/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(GLAQ.range)  <- 4267 

# clip the species range based on US map
GLAQ_clipped = st_intersection(GLAQ.range, states.map)

# Add the species name back to the map info
GLAQ_clipped$species = "Gleditsia aquatica"
all.data$species = "Gleditsia aquatica"

# clip all points outside of range with 50000 m buffer
GLAQ_flag = cc_iucn(x = all.data, range = GLAQ_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
GLAQ_occ_final = all.data[GLAQ_flag, ]

# export the values
write.csv(GLAQ_occ_final, file = "./Formatted.Data/Species.Nulls/Gleditsia.aquatica.null.csv")

#### Gleditsia triacanthos ####
# read in the range map of the first species
GLTR.range = st_read("../USTreeAtlas/shp/gledtria/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(GLTR.range)  <- 4267 

# clip the species range based on US map
GLTR_clipped = st_intersection(GLTR.range, states.map)

# Add the species name back to the map info
GLTR_clipped$species = "Gleditsia triacanthos"
all.data$species = "Gleditsia triacanthos"

# clip all points outside of range with 50000 m buffer
GLTR_flag = cc_iucn(x = all.data, range = GLTR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
GLTR_occ_final = all.data[GLTR_flag, ]

# export the values
write.csv(GLTR_occ_final, file = "./Formatted.Data/Species.Nulls/Gleditsia.triacanthos.null.csv")

#### Gordonia lasianthus ####
# read in the range map of the first species
GOLA.range = st_read("../USTreeAtlas/shp/gordlasi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(GOLA.range)  <- 4267 

# clip the species range based on US map
GOLA_clipped = st_intersection(GOLA.range, states.map)

# Add the species name back to the map info
GOLA_clipped$species = "Gordonia lasianthus"
all.data$species = "Gordonia lasianthus"

# clip all points outside of range with 50000 m buffer
GOLA_flag = cc_iucn(x = all.data, range = GOLA_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
GOLA_occ_final = all.data[GOLA_flag, ]

# export the values
write.csv(GOLA_occ_final, file = "./Formatted.Data/Species.Nulls/Gordonia.lasianthus.null.csv")

#### Gymnocladus dioicus ####
# read in the range map of the first species
GYDI.range = st_read("../USTreeAtlas/shp/gymndioi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(GYDI.range)  <- 4267 

# clip the species range based on US map
GYDI_clipped = st_intersection(GYDI.range, states.map)

# Add the species name back to the map info
GYDI_clipped$species = "Gymnocladus dioicus"
all.data$species = "Gymnocladus dioicus"

# clip all points outside of range with 50000 m buffer
GYDI_flag = cc_iucn(x = all.data, range = GYDI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
GYDI_occ_final = all.data[GYDI_flag, ]

# export the values
write.csv(GYDI_occ_final, file = "./Formatted.Data/Species.Nulls/Gymnocladus.dioicus.null.csv")

#### Ilex opaca ####
# read in the range map of the first species
ILOP.range = st_read("../USTreeAtlas/shp/ilexopac/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ILOP.range)  <- 4267 

# clip the species range based on US map
ILOP_clipped = st_intersection(ILOP.range, states.map)

# Add the species name back to the map info
ILOP_clipped$species = "Ilex opaca"
all.data$species = "Ilex opaca"

# clip all points outside of range with 50000 m buffer
ILOP_flag = cc_iucn(x = all.data, range = ILOP_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ILOP_occ_final = all.data[ILOP_flag, ]

# export the values
write.csv(ILOP_occ_final, file = "./Formatted.Data/Species.Nulls/Ilex.opaca.null.csv")

#### Juglans cinerea ####
# read in the range map of the first species
JUCI.range = st_read("../USTreeAtlas/shp/juglcine/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(JUCI.range)  <- 4267 

# clip the species range based on US map
JUCI_clipped = st_intersection(JUCI.range, states.map)

# Add the species name back to the map info
JUCI_clipped$species = "Juglans cinerea"
all.data$species = "Juglans cinerea"

# clip all points outside of range with 50000 m buffer
JUCI_flag = cc_iucn(x = all.data, range = JUCI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
JUCI_occ_final = all.data[JUCI_flag, ]

# export the values
write.csv(JUCI_occ_final, file = "./Formatted.Data/Species.Nulls/Juglans.cinerea.null.csv")

#### Juglans nigra ####
# read in the range map of the first species
JUNI.range = st_read("../USTreeAtlas/shp/juglnigr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(JUNI.range)  <- 4267 

# clip the species range based on US map
JUNI_clipped = st_intersection(JUNI.range, states.map)

# Add the species name back to the map info
JUNI_clipped$species = "Juglans nigra"
all.data$species = "Juglans nigra"

# clip all points outside of range with 50000 m buffer
JUNI_flag = cc_iucn(x = all.data, range = JUNI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
JUNI_occ_final = all.data[JUNI_flag, ]

# export the values
write.csv(JUNI_occ_final, file = "./Formatted.Data/Species.Nulls/Juglans.nigra.null.csv")

#### Juniperus virginiana ####
# read in the range map of the first species
JUVI.range = st_read("../USTreeAtlas/shp/junivirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(JUVI.range)  <- 4267 

# clip the species range based on US map
JUVI_clipped = st_intersection(JUVI.range, states.map)

# Add the species name back to the map info
JUVI_clipped$species = "Juniperus virginiana"
all.data$species = "Juniperus virginiana"

# clip all points outside of range with 50000 m buffer
JUVI_flag = cc_iucn(x = all.data, range = JUVI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
JUVI_occ_final = JUVI_clipped[JUVI_flag, ]

# export the values
write.csv(JUVI_occ_final, file = "./Formatted.Data/Species.Nulls/Juniperus.virginiana.null.csv")

#### Liquidambar styraciflua ####
# read in the range map of the first species
LIST2.range = st_read("../USTreeAtlas/shp/liqustyr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(LIST2.range)  <- 4267 

# clip the species range based on US map
LIST2_clipped = st_intersection(LIST2.range, states.map)

# Add the species name back to the map info
LIST2_clipped$species = "Liquidambar styraciflua"
all.data$species = "Liquidambar styraciflua"

# clip all points outside of range with 50000 m buffer
LIST2_flag = cc_iucn(x = all.data, range = LIST2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
LIST2_occ_final = all.data[LIST2_flag, ]

# export the values
write.csv(LIST2_occ_final, file = "./Formatted.Data/Species.Nulls/Liquidambar.styraciflua.null.csv")

#### Liriodendron tulipifera ####
# read in the range map of the first species
LITU.range = st_read("../USTreeAtlas/shp/lirituli/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(LITU.range)  <- 4267 

# clip the species range based on US map
LITU_clipped = st_intersection(LITU.range, states.map)

# Add the species name back to the map info
LITU_clipped$species = "Liriodendron tulipifera"
all.data$species = "Liriodendron tulipifera"

# clip all points outside of range with 50000 m buffer
LITU_flag = cc_iucn(x = all.data, range = LITU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
LITU_occ_final = all.data[LITU_flag, ]

# export the values
write.csv(LITU_occ_final, file = "./Formatted.Data/Species.Nulls/Liriodendron.tulipifera.null.csv")

#### Maclura pomifera ####
# read in the range map of the first species
MAPO.range = st_read("../USTreeAtlas/shp/maclpomi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAPO.range)  <- 4267 

# clip the species range based on US map
MAPO_clipped = st_intersection(MAPO.range, states.map)

# Add the species name back to the map info
MAPO_clipped$species = "Maclura pomifera"
all.data$species = "Maclura pomifera"

# clip all points outside of range with 50000 m buffer
MAPO_flag = cc_iucn(x = all.data, range = MAPO_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
MAPO_occ_final = all.data[MAPO_flag, ]

# export the values
write.csv(MAPO_occ_final, file = "./Formatted.Data/Species.Nulls/Maclura.pomifera.null.csv")

#### Magnolia acuminata ####
# read in the range map of the first species
MAAC.range = st_read("../USTreeAtlas/shp/magnacum/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAAC.range)  <- 4267 

# clip the species range based on US map
MAAC_clipped = st_intersection(MAAC.range, states.map)

# Add the species name back to the map info
MAAC_clipped$species = "Magnolia acuminata"
all.data$species = "Magnolia acuminata"

# clip all points outside of range with 50000 m buffer
MAAC_flag = cc_iucn(x = all.data, range = MAAC_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
MAAC_occ_final = all.data[MAAC_flag, ]

# export the values
write.csv(MAAC_occ_final, file = "./Formatted.Data/Species.Nulls/Magnolia.acuminata.null.csv")

#### Magnolia fraseri ####
# read in the range map of the first species
MAFR.range = st_read("../USTreeAtlas/shp/magnfras/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAFR.range)  <- 4267 

# clip the species range based on US map
MAFR_clipped = st_intersection(MAFR.range, states.map)

# Add the species name back to the map info
MAFR_clipped$species = "Magnolia fraseri"
all.data$species = "Magnolia fraseri"

# clip all points outside of range with 50000 m buffer
MAFR_flag = cc_iucn(x = all.data, range = MAFR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
MAFR_occ_final = all.data[MAFR_flag, ]

# export the values
write.csv(MAFR_occ_final, file = "./Formatted.Data/Species.Nulls/Magnolia.fraseri.null.csv")

#### Magnolia grandiflora ####
# read in the range map of the first species
MAGR4.range = st_read("../USTreeAtlas/shp/magngran/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAGR4.range)  <- 4267 

# clip the species range based on US map
MAGR4_clipped = st_intersection(MAGR4.range, states.map)

# Add the species name back to the map info
MAGR4_clipped$species = "Magnolia grandiflora"
all.data$species = "Magnolia grandiflora"

# clip all points outside of range with 50000 m buffer
MAGR4_flag = cc_iucn(x = all.data, range = MAGR4_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
MAGR4_occ_final = all.data[MAGR4_flag, ]

# export the values
write.csv(MAGR4_occ_final, file = "./Formatted.Data/Species.Nulls/Magnolia.grandiflora.null.csv")

#### Magnolia macrophylla ####
# read in the range map of the first species
MAMA2.range = st_read("../USTreeAtlas/shp/magnmacr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAMA2.range)  <- 4267 

# clip the species range based on US map
MAMA2_clipped = st_intersection(MAMA2.range, states.map)

# Add the species name back to the map info
MAMA2_clipped$species = "Magnolia macrophylla"
all.data$species = "Magnolia macrophylla"

# clip all points outside of range with 50000 m buffer
MAMA2_flag = cc_iucn(x = all.data, range = MAMA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
MAMA2_occ_final = all.data[MAMA2_flag, ]

# export the values
write.csv(MAMA2_occ_final, file = "./Formatted.Data/Species.Nulls/Magnolia.macrophylla.null.csv")

#### Magnolia virginiana ####
# read in the range map of the first species
MAVI2.range = st_read("../USTreeAtlas/shp/magnvirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MAVI2.range)  <- 4267 

# clip the species range based on US map
MAVI2_clipped = st_intersection(MAVI2.range, states.map)

# Add the species name back to the map info
MAVI2_clipped$species = "Magnolia virginiana"
all.data$species = "Magnolia virginiana"

# clip all points outside of range with 50000 m buffer
MAVI2_flag = cc_iucn(x = all.data, range = MAVI2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
MAVI2_occ_final = all.data[MAVI2_flag, ]

# export the values
write.csv(MAVI2_occ_final, file = "./Formatted.Data/Species.Nulls/Magnolia.virginiana.null.csv")

#### Morus rubra ####
# read in the range map of the first species
MORU2.range = st_read("../USTreeAtlas/shp/morurubr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(MORU2.range)  <- 4267 

# clip the species range based on US map
MORU2_clipped = st_intersection(MORU2.range, states.map)

# Add the species name back to the map info
MORU2_clipped$species = "Morus rubra"
all.data$species = "Morus rubra"

# clip all points outside of range with 50000 m buffer
MORU2_flag = cc_iucn(x = all.data, range = MORU2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
MORU2_occ_final = all.data[MORU2_flag, ]

# export the values
write.csv(MORU2_occ_final, file = "./Formatted.Data/Species.Nulls/Morus.rubra.null.csv")

#### Nyssaa aquatica ####
# read in the range map of the first species
NYAQ2.range = st_read("../USTreeAtlas/shp/nyssaqua/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(NYAQ2.range)  <- 4267 

# clip the species range based on US map
NYAQ2_clipped = st_intersection(NYAQ2.range, states.map)

# Add the species name back to the map info
NYAQ2_clipped$species = "Nyssaa aquatica"
all.data$species = "Nyssaa aquatica"

# clip all points outside of range with 50000 m buffer
NYAQ2_flag = cc_iucn(x = all.data, range = NYAQ2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
NYAQ2_occ_final = all.data[NYAQ2_flag, ]

# export the values
write.csv(NYAQ2_occ_final, file = "./Formatted.Data/Species.Nulls/Nyssa.aquatica.null.csv")

#### Nyssa ogeche ####
# read in the range map of the first species
NYOG.range = st_read("../USTreeAtlas/shp/nyssogec/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(NYOG.range)  <- 4267 

# clip the species range based on US map
NYOG_clipped = st_intersection(NYOG.range, states.map)

# Add the species name back to the map info
NYOG_clipped$species = "Nyssa ogeche"
all.data$species = "Nyssa ogeche"

# clip all points outside of range with 50000 m buffer
NYOG_flag = cc_iucn(x = all.data, range = NYOG_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
NYOG_occ_final = all.data[NYOG_flag, ]

# export the values
write.csv(NYOG_occ_final, file = "./Formatted.Data/Species.Nulls/Nyssa.ogeche.null.csv")

#### Nyssa sylvatica ####
# read in the range map of the first species
NYSY.range = st_read("../USTreeAtlas/shp/nysssylv/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(NYSY.range)  <- 4267 

# clip the species range based on US map
NYSY_clipped = st_intersection(NYSY.range, states.map)

# Add the species name back to the map info
NYSY_clipped$species = "Nyssa sylvatica"
all.data$species = "Nyssa sylvatica"

# clip all points outside of range with 50000 m buffer
NYSY_flag = cc_iucn(x = all.data, range = NYSY_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
NYSY_occ_final = all.data[NYSY_flag, ]

# export the values
write.csv(NYSY_occ_final, file = "./Formatted.Data/Species.Nulls/Nyssa.sylvatica.null.csv")

#### Ostrya virginiana ####
# read in the range map of the first species
OSVI.range = st_read("../USTreeAtlas/shp/ostrvirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(OSVI.range)  <- 4267 

# clip the species range based on US map
OSVI_clipped = st_intersection(OSVI.range, states.map)

# Add the species name back to the map info
OSVI_clipped$species = "Ostrya virginiana"
all.data$species = "Ostrya virginiana"

# clip all points outside of range with 50000 m buffer
OSVI_flag = cc_iucn(x = all.data, range = OSVI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
OSVI_occ_final = all.data[OSVI_flag, ]

# export the values
write.csv(OSVI_occ_final, file = "./Formatted.Data/Species.Nulls/Ostrya.virginiana.null.csv")

#### Oxydendrum arboreum ####
# read in the range map of the first species
OXAR.range = st_read("../USTreeAtlas/shp/oxydarbo/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(OXAR.range)  <- 4267 

# clip the species range based on US map
OXAR_clipped = st_intersection(OXAR.range, states.map)

# Add the species name back to the map info
OXAR_clipped$species = "Oxydendrum arboreum"
all.data$species = "Oxydendrum arboreum"

# clip all points outside of range with 50000 m buffer
OXAR_flag = cc_iucn(x = all.data, range = OXAR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
OXAR_occ_final = all.data[OXAR_flag, ]

# export the values
write.csv(OXAR_occ_final, file = "./Formatted.Data/Species.Nulls/Oxydendrum.arboreum.null.csv")

#### Paulownia tomentosa ####
# clip the species range based on US map, from BIEN
(PATO2.range.sf <- BIEN_ranges_load_species('Paulownia tomentosa'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
PATO2.range.2 = PATO2.range.sf %>%
  st_transform(st_crs(states.map))

PATO2.range.3 = terra::vect(PATO2.range.2)
states.map.2 = terra::vect(states.map)
PATO2.range.4 = terra::intersect(PATO2.range.3, states.map.2)
PATO2.range.5 = st_as_sf(PATO2.range.4)

# Add the species name back to the map info
PATO2.range.5$species = "Paulownia tomentosa"
all.data$species = "Paulownia tomentosa"

# clip all points outside of range with 50000 m buffer
PATO2_flag = cc_iucn(x = all.data, range = PATO2.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PATO2_occ_final = all.data[PATO2_flag, ]

# export the values
write.csv(PATO2_occ_final, file = "./Formatted.Data/Species.Nulls/Paulownia.tomentosa.null.csv")

#### Persea borbonia ####
# read in the range map of the first species
PEBO.range = st_read("../USTreeAtlas/shp/persborb/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PEBO.range)  <- 4267 

# clip the species range based on US map
PEBO_clipped = st_intersection(PEBO.range, states.map)

# Add the species name back to the map info
PEBO_clipped$species = "Persea borbonia"
all.data$species = "Persea borbonia"

# clip all points outside of range with 50000 m buffer
PEBO_flag = cc_iucn(x = all.data, range = PEBO_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PEBO_occ_final = all.data[PEBO_flag, ]

# export the values
write.csv(PEBO_occ_final, file = "./Formatted.Data/Species.Nulls/Persea.borbonia.null.csv")

#### Pinus clausa ####
# read in the range map of the first species
PICL.range = st_read("../USTreeAtlas/shp/pinuclau/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PICL.range)  <- 4267 

# clip the species range based on US map
PICL_clipped = st_intersection(PICL.range, states.map)

# Add the species name back to the map info
PICL_clipped$species = "Pinus clausa"
all.data$species = "Pinus clausa"

# clip all points outside of range with 50000 m buffer
PICL_flag = cc_iucn(x = all.data, range = PICL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PICL_occ_final = all.data[PICL_flag, ]

# export the values
write.csv(PICL_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.clausa.null.csv")

#### Pinus echinata ####
# read in the range map of the first species
PIEC2.range = st_read("../USTreeAtlas/shp/pinuechi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIEC2.range)  <- 4267 

# clip the species range based on US map
PIEC2_clipped = st_intersection(PIEC2.range, states.map)

# Add the species name back to the map info
PIEC2_clipped$species = "Pinus echinata"
all.data$species = "Pinus echinata"

# clip all points outside of range with 50000 m buffer
PIEC2_flag = cc_iucn(x = all.data, range = PIEC2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PIEC2_occ_final = all.data[PIEC2_flag, ]

# export the values
write.csv(PIEC2_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.echinata.null.csv")

#### Pinus elliottii ####
# read in the range map of the first species
PIEL.range = st_read("../USTreeAtlas/shp/pinuelli/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIEL.range)  <- 4267 

# clip the species range based on US map
PIEL_clipped = st_intersection(PIEL.range, states.map)

# Add the species name back to the map info
PIEL_clipped$species = "Pinus elliottii"
all.data$species = "Pinus elliottii"

# clip all points outside of range with 50000 m buffer
PIEL_flag = cc_iucn(x = all.data, range = PIEL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PIEL_occ_final = all.data[PIEL_flag, ]

# export the values
write.csv(PIEL_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.elliottii.null.csv")

#### Pinus glabra ####
# read in the range map of the first species
PIGL2.range = st_read("../USTreeAtlas/shp/pinuglab/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIGL2.range)  <- 4267 

# clip the species range based on US map
PIGL2_clipped = st_intersection(PIGL2.range, states.map)

# Add the species name back to the map info
PIGL2_clipped$species = "Pinus glabra"
all.data$species = "Pinus glabra"

# clip all points outside of range with 50000 m buffer
PIGL2_flag = cc_iucn(x = all.data, range = PIGL2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PIGL2_occ_final = all.data[PIGL2_flag, ]

# export the values
write.csv(PIGL2_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.glabra.null.csv")

#### Pinus palustris ####
# read in the range map of the first species
PIPA2.range = st_read("../USTreeAtlas/shp/pinupalu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIPA2.range)  <- 4267 

# clip the species range based on US map
PIPA2_clipped = st_intersection(PIPA2.range, states.map)

# Add the species name back to the map info
PIPA2_clipped$species = "Pinus palustris"
all.data$species = "Pinus palustris"

# clip all points outside of range with 50000 m buffer
PIPA2_flag = cc_iucn(x = all.data, range = PIPA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PIPA2_occ_final = all.data[PIPA2_flag, ]

# export the values
write.csv(PIPA2_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.palustris.null.csv")

#### Pinus pungens ####
# read in the range map of the first species
PIPUS.range = st_read("../USTreeAtlas/shp/pinupung/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIPUS.range)  <- 4267 

# clip the species range based on US map
PIPUS_clipped = st_intersection(PIPUS.range, states.map)

# Add the species name back to the map info
PIPUS_clipped$species = "Pinus pungens"
all.data$species = "Pinus pungens"

# clip all points outside of range with 50000 m buffer
PIPUS_flag = cc_iucn(x = all.data, range = PIPUS_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PIPUS_occ_final = all.data[PIPUS_flag, ]

# export the values
write.csv(PIPUS_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.pungens.null.csv")

#### Pinus rigida ####
# read in the range map of the first species
PIRI.range = st_read("../USTreeAtlas/shp/pinurigi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIRI.range)  <- 4267 

# clip the species range based on US map
PIRI_clipped = st_intersection(PIRI.range, states.map)

# Add the species name back to the map info
PIRI_clipped$species = "Pinus rigida"
all.data$species = "Pinus rigida"

# clip all points outside of range with 50000 m buffer
PIRI_flag = cc_iucn(x = all.data, range = PIRI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PIRI_occ_final = all.data[PIRI_flag, ]

# export the values
write.csv(PIRI_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.rigida.null.csv")

#### Pinus serotina ####
# read in the range map of the first species
PISE.range = st_read("../USTreeAtlas/shp/pinusero/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PISE.range)  <- 4267 

# clip the species range based on US map
PISE_clipped = st_intersection(PISE.range, states.map)

# Add the species name back to the map info
PISE_clipped$species = "Pinus serotina"
all.data$species = "Pinus serotina"

# clip all points outside of range with 50000 m buffer
PISE_flag = cc_iucn(x = all.data, range = PISE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PISE_occ_final = all.data[PISE_flag, ]

# export the values
write.csv(PISE_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.serotina.null.csv")

#### Pinus taeda ####
# read in the range map of the first species
PITA.range = st_read("../USTreeAtlas/shp/pinutaed/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PITA.range)  <- 4267 

# clip the species range based on US map
PITA_clipped = st_intersection(PITA.range, states.map)

# Add the species name back to the map info
PITA_clipped$species = "Pinus taeda"
all.data$species = "Pinus taeda"

# clip all points outside of range with 50000 m buffer
PITA_flag = cc_iucn(x = all.data, range = PITA_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PITA_occ_final = all.data[PITA_flag, ]

# export the values
write.csv(PITA_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.taeda.null.csv")

#### Pinus virginiana ####
# read in the range map of the first species
PIVI2.range = st_read("../USTreeAtlas/shp/pinuvirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PIVI2.range)  <- 4267 

# clip the species range based on US map
PIVI2_clipped = st_intersection(PIVI2.range, states.map)

# Add the species name back to the map info
PIVI2_clipped$species = "Pinus virginiana"
all.data$species = "Pinus virginiana"

# clip all points outside of range with 50000 m buffer
PIVI2_flag = cc_iucn(x = all.data, range = PIVI2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PIVI2_occ_final = all.data[PIVI2_flag, ]

# export the values
write.csv(PIVI2_occ_final, file = "./Formatted.Data/Species.Nulls/Pinus.virginiana.null.csv")

#### Planera aquatica ####
# read in the range map of the first species
PLAQ.range = st_read("../USTreeAtlas/shp/planaqua/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PLAQ.range)  <- 4267 

# clip the species range based on US map
PLAQ_clipped = st_intersection(PLAQ.range, states.map)

# Add the species name back to the map info
PLAQ_clipped$species = "Planera aquatica"
all.data$species = "Planera aquatica"

# clip all points outside of range with 50000 m buffer
PLAQ_flag = cc_iucn(x = all.data, range = PLAQ_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PLAQ_occ_final = all.data[PLAQ_flag, ]

# export the values
write.csv(PLAQ_occ_final, file = "./Formatted.Data/Species.Nulls/Planera.aquatica.null.csv")

#### Platanus occidentalis ####
# read in the range map of the first species
PLOC.range = st_read("../USTreeAtlas/shp/platocci/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PLOC.range)  <- 4267 

# clip the species range based on US map
PLOC_clipped = st_intersection(PLOC.range, states.map)

# Add the species name back to the map info
PLOC_clipped$species = "Platanus occidentalis"
all.data$species = "Platanus occidentalis"

# clip all points outside of range with 50000 m buffer
PLOC_flag = cc_iucn(x = all.data, range = PLOC_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PLOC_occ_final = all.data[PLOC_flag, ]

# export the values
write.csv(PLOC_occ_final, file = "./Formatted.Data/Species.Nulls/Platanus.occidentalis.null.csv")

#### Populus deltoides ####
# read in the range map of the first species
PODE3.range = st_read("../USTreeAtlas/shp/popudelt/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PODE3.range)  <- 4267 

# clip the species range based on US map
PODE3_clipped = st_intersection(PODE3.range, states.map)

# Add the species name back to the map info
PODE3_clipped$species = "Populus deltoides"
all.data$species = "Populus deltoides"

# clip all points outside of range with 50000 m buffer
PODE3_flag = cc_iucn(x = all.data, range = PODE3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
PODE3_occ_final = all.data[PODE3_flag, ]

# export the values
write.csv(PODE3_occ_final, file = "./Formatted.Data/Species.Nulls/Populus.deltoides.null.csv")

#### Populus heterophylla ####
# read in the range map of the first species
POHE4.range = st_read("../USTreeAtlas/shp/popuhete/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(POHE4.range)  <- 4267 

# clip the species range based on US map
POHE4_clipped = st_intersection(POHE4.range, states.map)

# Add the species name back to the map info
POHE4_clipped$species = "Populus heterophylla"
all.data$species = "Populus heterophylla"

# clip all points outside of range with 50000 m buffer
POHE4_flag = cc_iucn(x = all.data, range = POHE4_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
POHE4_occ_final = all.data[POHE4_flag, ]

# export the values
write.csv(POHE4_occ_final, file = "./Formatted.Data/Species.Nulls/Populus.heterophylla.null.csv")

#### Prunus americana ####
# read in the range map of the first species
PRAM.range = st_read("../USTreeAtlas/shp/prunamer/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PRAM.range)  <- 4267 

# clip the species range based on US map
PRAM_clipped = st_intersection(PRAM.range, states.map)

# Add the species name back to the map info
PRAM_clipped$species = "Prunus americana"
all.data$species = "Prunus americana"

# clip all points outside of range with 50000 m buffer
PRAM_flag = cc_iucn(x = all.data, range = PRAM_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
PRAM_occ_final = all.data[PRAM_flag, ]

# export the values
write.csv(PRAM_occ_final, file = "./Formatted.Data/Species.Nulls/Prunus.americana.null.csv")

#### Prunus serotina ####
# read in the range map of the first species
PRSE2.range = st_read("../USTreeAtlas/shp/prunsero/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(PRSE2.range)  <- 4267 

# clip the species range based on US map
PRSE2_clipped = st_intersection(PRSE2.range, states.map)

# Add the species name back to the map info
PRSE2_clipped$species = "Prunus serotina"
all.data$species = "Prunus serotina"

# clip all points outside of range with 50000 m buffer
PRSE2_flag = cc_iucn(x = all.data, range = PRSE2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                      value = "flagged", buffer = 50000)
# remove the points
PRSE2_occ_final = all.data[PRSE2_flag, ]

# export the values
write.csv(PRSE2_occ_final, file = "./Formatted.Data/Species.Nulls/Prunus.serotina.null.csv")

#### Quercus alba ####
# read in the range map of the first species
QUAL.range = st_read("../USTreeAtlas/shp/queralba/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUAL.range)  <- 4267 

# clip the species range based on US map
QUAL_clipped = st_intersection(QUAL.range, states.map)

# Add the species name back to the map info
QUAL_clipped$species = "Quercus alba"
all.data$species = "Quercus alba"

# clip all points outside of range with 50000 m buffer
QUAL_flag = cc_iucn(x = all.data, range = QUAL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUAL_occ_final = all.data[QUAL_flag, ]

# export the values
write.csv(QUAL_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.alba.null.csv")

#### Quercus bicolor ####
# read in the range map of the first species
QUBI.range = st_read("../USTreeAtlas/shp/querbico/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUBI.range)  <- 4267 

# clip the species range based on US map
QUBI_clipped = st_intersection(QUBI.range, states.map)

# Add the species name back to the map info
QUBI_clipped$species = "Quercus bicolor"
all.data$species = "Quercus bicolor"

# clip all points outside of range with 50000 m buffer
QUBI_flag = cc_iucn(x = all.data, range = QUBI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUBI_occ_final = all.data[QUBI_flag, ]

# export the values
write.csv(QUBI_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.bicolor.null.csv")

#### Quercus coccinea ####
# read in the range map of the first species
QUCO2.range = st_read("../USTreeAtlas/shp/quercocc/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUCO2.range)  <- 4267 

# clip the species range based on US map
QUCO2_clipped = st_intersection(QUCO2.range, states.map)

# Add the species name back to the map info
QUCO2_clipped$species = "Quercus coccinea"
all.data$species = "Quercus coccinea"

# clip all points outside of range with 50000 m buffer
QUCO2_flag = cc_iucn(x = all.data, range = QUCO2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUCO2_occ_final = all.data[QUCO2_flag, ]

# export the values
write.csv(QUCO2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.coccinea.null.csv")

#### Quercus ellipsoidalis ####
# read in the range map of the first species
QUEL.range = st_read("../USTreeAtlas/shp/querelli/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUEL.range)  <- 4267 

# clip the species range based on US map
QUEL_clipped = st_intersection(QUEL.range, states.map)

# Add the species name back to the map info
QUEL_clipped$species = "Quercus ellipsoidalis"
all.data$species = "Quercus ellipsoidalis"

# clip all points outside of range with 50000 m buffer
QUEL_flag = cc_iucn(x = all.data, range = QUEL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUEL_occ_final = all.data[QUEL_flag, ]

# export the values
write.csv(QUEL_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.ellipsoidalis.null.csv")

#### Quercus falcata ####
# read in the range map of the first species
QUFA.range = st_read("../USTreeAtlas/shp/querfalc/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUFA.range)  <- 4267 

# clip the species range based on US map
QUFA_clipped = st_intersection(QUFA.range, states.map)

# Add the species name back to the map info
QUFA_clipped$species = "Quercus falcata"
all.data$species = "Quercus falcata"

# clip all points outside of range with 50000 m buffer
QUFA_flag = cc_iucn(x = all.data, range = QUFA_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUFA_occ_final = all.data[QUFA_flag, ]

# export the values
write.csv(QUFA_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.falcata.null.csv")

#### Quercus ilicifolia ####
# read in the range map of the first species
QUIL.range = st_read("../USTreeAtlas/shp/querilic/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUIL.range)  <- 4267 

# clip the species range based on US map
QUIL_clipped = st_intersection(QUIL.range, states.map)

# Add the species name back to the map info
QUIL_clipped$species = "Quercus ilicifolia"
all.data$species = "Quercus ilicifolia"

# clip all points outside of range with 50000 m buffer
QUIL_flag = cc_iucn(x = all.data, range = QUIL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUIL_occ_final = all.data[QUIL_flag, ]

# export the values
write.csv(QUIL_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.ilicifolia.null.csv")

#### Quercus imbricaria ####
# read in the range map of the first species
QUIM.range = st_read("../USTreeAtlas/shp/querimbr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUIM.range)  <- 4267 

# clip the species range based on US map
QUIM_clipped = st_intersection(QUIM.range, states.map)

# Add the species name back to the map info
QUIM_clipped$species = "Quercus imbricaria"
all.data$species = "Quercus imbricaria"

# clip all points outside of range with 50000 m buffer
QUIM_flag = cc_iucn(x = all.data, range = QUIM_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUIM_occ_final = all.data[QUIM_flag, ]

# export the values
write.csv(QUIM_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.imbricaria.null.csv")

#### Quercus incana ####
# read in the range map of the first species
QUIN.range = st_read("../USTreeAtlas/shp/querinca/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUIN.range)  <- 4267 

# clip the species range based on US map
QUIN_clipped = st_intersection(QUIN.range, states.map)

# Add the species name back to the map info
QUIN_clipped$species = "Quercus incana"
all.data$species = "Quercus incana"

# clip all points outside of range with 50000 m buffer
QUIN_flag = cc_iucn(x = all.data, range = QUIN_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUIN_occ_final = all.data[QUIN_flag, ]

# export the values
write.csv(QUIN_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.incana.null.csv")

#### Quercus laevis ####
# read in the range map of the first species
QULA2.range = st_read("../USTreeAtlas/shp/querlaev/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QULA2.range)  <- 4267 

# clip the species range based on US map
QULA2_clipped = st_intersection(QULA2.range, states.map)

# Add the species name back to the map info
QULA2_clipped$species = "Quercus laevis"
all.data$species = "Quercus laevis"

# clip all points outside of range with 50000 m buffer
QULA2_flag = cc_iucn(x = all.data, range = QULA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QULA2_occ_final = all.data[QULA2_flag, ]

# export the values
write.csv(QULA2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.laevis.null.csv")

#### Quercus laurifolia ####
# read in the range map of the first species
QULA3.range = st_read("../USTreeAtlas/shp/querlaur/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QULA3.range)  <- 4267 

# clip the species range based on US map
QULA3_clipped = st_intersection(QULA3.range, states.map)

# Add the species name back to the map info
QULA3_clipped$species = "Quercus laurifolia"
all.data$species = "Quercus laurifolia"

# clip all points outside of range with 50000 m buffer
QULA3_flag = cc_iucn(x = all.data, range = QULA3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QULA3_occ_final = all.data[QULA3_flag, ]

# export the values
write.csv(QULA3_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.laurifolia.null.csv")

#### Quercus lyrata ####
# read in the range map of the first species
QULY.range = st_read("../USTreeAtlas/shp/querlyra/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QULY.range)  <- 4267 

# clip the species range based on US map
QULY_clipped = st_intersection(QULY.range, states.map)

# Add the species name back to the map info
QULY_clipped$species = "Quercus lyrata"
all.data$species = "Quercus lyrata"

# clip all points outside of range with 50000 m buffer
QULY_flag = cc_iucn(x = all.data, range = QULY_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QULY_occ_final = all.data[QULY_flag, ]

# export the values
write.csv(QULY_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.lyrata.null.csv")

#### Quercus macrocarpa ####
# read in the range map of the first species
QUMA2.range = st_read("../USTreeAtlas/shp/quermacr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUMA2.range)  <- 4267 

# clip the species range based on US map
QUMA2_clipped = st_intersection(QUMA2.range, states.map)

# Add the species name back to the map info
QUMA2_clipped$species = "Quercus macrocarpa"
all.data$species = "Quercus macrocarpa"

# clip all points outside of range with 50000 m buffer
QUMA2_flag = cc_iucn(x = all.data, range = QUMA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUMA2_occ_final = all.data[QUMA2_flag, ]

# export the values
write.csv(QUMA2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.macrocarpa.null.csv")

#### Quercus margaretta ####
# clip the species range based on US map, from BIEN
(QUMA6.range.sf <- BIEN_ranges_load_species('Quercus margarettae'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
QUMA6.range.2 = QUMA6.range.sf %>%
  st_transform(st_crs(states.map))

QUMA6.range.3 = terra::vect(QUMA6.range.2)
states.map.2 = terra::vect(states.map)
QUMA6.range.4 = terra::intersect(QUMA6.range.3, states.map.2)
QUMA6.range.5 = st_as_sf(QUMA6.range.4)

# Add the species name back to the map info
QUMA6.range.5$species = "Quercus margaretta"
all.data$species = "Quercus margaretta"

# clip all points outside of range with 50000 m buffer
QUMA6_flag = cc_iucn(x = all.data, range = QUMA6.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUMA6_occ_final = all.data[QUMA6_flag, ]

# export the values
write.csv(QUMA6_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.margaretta.null.csv")

#### Quercus marilandica ####
# read in the range map of the first species
QUMA3.range = st_read("../USTreeAtlas/shp/quermari/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUMA3.range)  <- 4267 

# clip the species range based on US map
QUMA3_clipped = st_intersection(QUMA3.range, states.map)

# Add the species name back to the map info
QUMA3_clipped$species = "Quercus marilandica"
all.data$species = "Quercus marilandica"

# clip all points outside of range with 50000 m buffer
QUMA3_flag = cc_iucn(x = all.data, range = QUMA3_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUMA3_occ_final = all.data[QUMA3_flag, ]

# export the values
write.csv(QUMA3_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.marilandica.null.csv")

#### Quercus michauxii ####
# includes range map of Q. michauxii and Q. prinus
# read in the range map of the first species
QUMI.range = st_read("../USTreeAtlas/shp/quermich/")
QUPR2.range = st_read("../USTreeAtlas/shp/querprin/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUMI.range)  <- 4267 
st_crs(QUPR2.range)  <- 4267

QUMI_QUPR2.range = bind_rows(QUMI.range, QUPR2.range)

# clip the species range based on US map
QUMI_QUPR2_clipped = st_intersection(QUMI_QUPR2.range, states.map)

# Add the species name back to the map info
QUMI_QUPR2_clipped$species = "Quercus michauxii"
all.data$species = "Quercus michauxii"

# clip all points outside of range with 50000 m buffer
QUMI_QUPR2_flag = cc_iucn(x = all.data, range = QUMI_QUPR2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                          value = "flagged", buffer = 50000)
# remove the points
QUMI_QUPR2_occ_final = all.data[QUMI_QUPR2_flag, ]

# export the values
write.csv(QUMI_QUPR2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.michauxii.null.csv")

#### Quercus minima ####
# clip the species range based on US map, from BIEN
(QUMI2.range.sf <- BIEN_ranges_load_species('Quercus minima'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
QUMI2.range.2 = QUMI2.range.sf %>%
  st_transform(st_crs(states.map))

QUMI2.range.3 = terra::vect(QUMI2.range.2)
states.map.2 = terra::vect(states.map)
QUMI2.range.4 = terra::intersect(QUMI2.range.3, states.map.2)
QUMI2.range.5 = st_as_sf(QUMI2.range.4)

# Add the species name back to the map info
QUMI2.range.5$species = "Quercus minima"
all.data$species = "Quercus minima"

# clip all points outside of range with 50000 m buffer
QUMI2_flag = cc_iucn(x = all.data, range = QUMI2.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUMI2_occ_final = all.data[QUMI2_flag, ]

# export the values
write.csv(QUMI2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.minima.null.csv")

#### Quercus muehlenbergii ####
# read in the range map of the first species
QUMU.range = st_read("../USTreeAtlas/shp/quermueh/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUMU.range)  <- 4267 

# clip the species range based on US map
QUMU_clipped = st_intersection(QUMU.range, states.map)

# Add the species name back to the map info
QUMU_clipped$species = "Quercus muehlenbergii"
all.data$species = "Quercus muehlenbergii"

# clip all points outside of range with 50000 m buffer
QUMU_flag = cc_iucn(x = all.data, range = QUMU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUMU_occ_final = all.data[QUMU_flag, ]

# export the values
write.csv(QUMU_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.muehlenbergii.null.csv")

#### Quercus nigra ####
# read in the range map of the first species
QUNI.range = st_read("../USTreeAtlas/shp/quernigr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUNI.range)  <- 4267 

# clip the species range based on US map
QUNI_clipped = st_intersection(QUNI.range, states.map)

# Add the species name back to the map info
QUNI_clipped$species = "Quercus nigra"
all.data$species = "Quercus nigra"

# clip all points outside of range with 50000 m buffer
QUNI_flag = cc_iucn(x = all.data, range = QUNI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUNI_occ_final = all.data[QUNI_flag, ]

# export the values
write.csv(QUNI_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.nigra.null.csv")

#### Quercus pagoda ####
# clip the species range based on US map, from BIEN
(QUPA5.range.sf <- BIEN_ranges_load_species('Quercus pagoda'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
QUPA5.range.2 = QUPA5.range.sf %>%
  st_transform(st_crs(states.map))

QUPA5.range.3 = terra::vect(QUPA5.range.2)
states.map.2 = terra::vect(states.map)
QUPA5.range.4 = terra::intersect(QUPA5.range.3, states.map.2)
QUPA5.range.5 = st_as_sf(QUPA5.range.4)

# Add the species name back to the map info
QUPA5.range.5$species = "Quercus pagoda"
all.data$species = "Quercus pagoda"

# clip all points outside of range with 50000 m buffer
QUPA5_flag = cc_iucn(x = all.data, range = QUPA5.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUPA5_occ_final = all.data[QUPA5_flag, ]

# export the values
write.csv(QUPA5_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.pagoda.null.csv")

#### Quercus palustris ####
# read in the range map of the first species
QUPA2.range = st_read("../USTreeAtlas/shp/querpalu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUPA2.range)  <- 4267 

# clip the species range based on US map
QUPA2_clipped = st_intersection(QUPA2.range, states.map)

# Add the species name back to the map info
QUPA2_clipped$species = "Quercus palustris"
all.data$species = "Quercus palustris"

# clip all points outside of range with 50000 m buffer
QUPA2_flag = cc_iucn(x = all.data, range = QUPA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUPA2_occ_final = all.data[QUPA2_flag, ]

# export the values
write.csv(QUPA2_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.palustris.null.csv")

#### Quercus phellos ####
# read in the range map of the first species
QUPH.range = st_read("../USTreeAtlas/shp/querphel/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUPH.range)  <- 4267 

# clip the species range based on US map
QUPH_clipped = st_intersection(QUPH.range, states.map)

# Add the species name back to the map info
QUPH_clipped$species = "Quercus phellos"
all.data$species = "Quercus phellos"

# clip all points outside of range with 50000 m buffer
QUPH_flag = cc_iucn(x = all.data, range = QUPH_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUPH_occ_final = all.data[QUPH_flag, ]

# export the values
write.csv(QUPH_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.phellos.null.csv")

#### Quercus rubra ####
# read in the range map of the first species
QURU.range = st_read("../USTreeAtlas/shp/querrubr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QURU.range)  <- 4267 

# clip the species range based on US map
QURU_clipped = st_intersection(QURU.range, states.map)

# Add the species name back to the map info
QURU_clipped$species = "Quercus rubra"
all.data$species = "Quercus rubra"

# clip all points outside of range with 50000 m buffer
QURU_flag = cc_iucn(x = all.data, range = QURU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QURU_occ_final = all.data[QURU_flag, ]

# export the values
write.csv(QURU_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.rubra.null.csv")

#### Quercus shumardii ####
# read in the range map of the first species
QUSH.range = st_read("../USTreeAtlas/shp/quershum/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUSH.range)  <- 4267 

# clip the species range based on US map
QUSH_clipped = st_intersection(QUSH.range, states.map)

# Add the species name back to the map info
QUSH_clipped$species = "Quercus shumardii"
all.data$species = "Quercus shumardii"

# clip all points outside of range with 50000 m buffer
QUSH_flag = cc_iucn(x = all.data, range = QUSH_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUSH_occ_final = all.data[QUSH_flag, ]

# export the values
write.csv(QUSH_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.shumardii.null.csv")

#### Quercus sinuata ####
# clip the species range based on US map, from BIEN
(QUSIS.range.sf <- BIEN_ranges_load_species('Quercus sinuata'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
QUSIS.range.2 = QUSIS.range.sf %>%
  st_transform(st_crs(states.map))

QUSIS.range.3 = terra::vect(QUSIS.range.2)
states.map.2 = terra::vect(states.map)
QUSIS.range.4 = terra::intersect(QUSIS.range.3, states.map.2)
QUSIS.range.5 = st_as_sf(QUSIS.range.4)

# Add the species name back to the map info
QUSIS.range.5$species = "Quercus sinuata"
all.data$species = "Quercus sinuata"

# clip all points outside of range with 50000 m buffer
QUSIS_flag = cc_iucn(x = all.data, range = QUSIS.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
QUSIS_occ_final = all.data[QUSIS_flag, ]

# export the values
write.csv(QUSIS_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.sinuata.null.csv")

#### Quercus stellata ####
# read in the range map of the first species
QUST.range = st_read("../USTreeAtlas/shp/querstel/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUST.range)  <- 4267 

# clip the species range based on US map
QUST_clipped = st_intersection(QUST.range, states.map)

# Add the species name back to the map info
QUST_clipped$species = "Quercus stellata"
all.data$species = "Quercus stellata"

# clip all points outside of range with 50000 m buffer
QUST_flag = cc_iucn(x = all.data, range = QUST_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUST_occ_final = all.data[QUST_flag, ]

# export the values
write.csv(QUST_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.stellata.null.csv")

#### Quercus texana ####
# clip the species range based on US map, from BIEN
(QUTE.range.sf <- BIEN_ranges_load_species('Quercus texana'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
QUTE.range.2 = QUTE.range.sf %>%
  st_transform(st_crs(states.map))

QUTE.range.3 = terra::vect(QUTE.range.2)
states.map.2 = terra::vect(states.map)
QUTE.range.4 = terra::intersect(QUTE.range.3, states.map.2)
QUTE.range.5 = st_as_sf(QUTE.range.4)

# Add the species name back to the map info
QUTE.range.5$species = "Quercus texana"
all.data$species = "Quercus texana"

# clip all points outside of range with 50000 m buffer
QUTE_flag = cc_iucn(x = all.data, range = QUTE.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUTE_occ_final = all.data[QUTE_flag, ]

# export the values
write.csv(QUTE_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.texana.null.csv")

#### Quercus velutina ####
# read in the range map of the first species
QUVE.range = st_read("../USTreeAtlas/shp/quervelu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUVE.range)  <- 4267 

# clip the species range based on US map
QUVE_clipped = st_intersection(QUVE.range, states.map)

# Add the species name back to the map info
QUVE_clipped$species = "Quercus velutina"
all.data$species = "Quercus velutina"

# clip all points outside of range with 50000 m buffer
QUVE_flag = cc_iucn(x = all.data, range = QUVE_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUVE_occ_final = all.data[QUVE_flag, ]

# export the values
write.csv(QUVE_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.velutina.null.csv")

#### Quercus virginiana ####
# read in the range map of the first species
QUVI.range = st_read("../USTreeAtlas/shp/quervirg/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(QUVI.range)  <- 4267 

# clip the species range based on US map
QUVI_clipped = st_intersection(QUVI.range, states.map)

# Add the species name back to the map info
QUVI_clipped$species = "Quercus virginiana"
all.data$species = "Quercus virginiana"

# clip all points outside of range with 50000 m buffer
QUVI_flag = cc_iucn(x = all.data, range = QUVI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
QUVI_occ_final = all.data[QUVI_flag, ]

# export the values
write.csv(QUVI_occ_final, file = "./Formatted.Data/Species.Nulls/Quercus.virginiana.null.csv")

#### Robinia pseudoacacia ####
# read in the range map of the first species
ROPS.range = st_read("../USTreeAtlas/shp/robipseu/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ROPS.range)  <- 4267 

# clip the species range based on US map
ROPS_clipped = st_intersection(ROPS.range, states.map)

# Add the species name back to the map info
ROPS_clipped$species = "Robinia pseudoacacia"
all.data$species = "Robinia pseudoacacia"

# clip all points outside of range with 50000 m buffer
ROPS_flag = cc_iucn(x = all.data, range = ROPS_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ROPS_occ_final = all.data[ROPS_flag, ]

# export the values
write.csv(ROPS_occ_final, file = "./Formatted.Data/Species.Nulls/Robinia.pseudoacacia.null.csv")

#### Salix caroliniana ####
# read in the range map of the first species
CACA5.range = st_read("../USTreeAtlas/shp/salicaro/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(CACA5.range)  <- 4267 

# clip the species range based on US map
CACA5_clipped = st_intersection(CACA5.range, states.map)

# Add the species name back to the map info
CACA5_clipped$species = "Salix caroliniana"
all.data$species = "Salix caroliniana"

# clip all points outside of range with 50000 m buffer
CACA5_flag = cc_iucn(x = all.data, range = CACA5_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
CACA5_occ_final = all.data[CACA5_flag, ]

# export the values
write.csv(CACA5_occ_final, file = "./Formatted.Data/Species.Nulls/Salix.caroliniana.null.csv")

#### Salix nigra ####
# read in the range map of the first species
SANI.range = st_read("../USTreeAtlas/shp/salinigr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(SANI.range)  <- 4267 

# clip the species range based on US map
SANI_clipped = st_intersection(SANI.range, states.map)

# Add the species name back to the map info
SANI_clipped$species = "Salix nigra"
all.data$species = "Salix nigra"

# clip all points outside of range with 50000 m buffer
SANI_flag = cc_iucn(x = all.data, range = SANI_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
SANI_occ_final = all.data[SANI_flag, ]

# export the values
write.csv(SANI_occ_final, file = "./Formatted.Data/Species.Nulls/Salix.nigra.null.csv")

#### Sassafras albidum ####
# read in the range map of the first species
SAAL5.range = st_read("../USTreeAtlas/shp/sassalbi/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(SAAL5.range)  <- 4267 

# clip the species range based on US map
SAAL5_clipped = st_intersection(SAAL5.range, states.map)

# Add the species name back to the map info
SAAL5_clipped$species = "Sassafras albidum"
all.data$species = "Sassafras albidum"

# clip all points outside of range with 50000 m buffer
SAAL5_flag = cc_iucn(x = all.data, range = SAAL5_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
SAAL5_occ_final = all.data[SAAL5_flag, ]

# export the values
write.csv(SAAL5_occ_final, file = "./Formatted.Data/Species.Nulls/Sassafras.albidum.null.csv")

#### Taxodium distichum ####
# read in the range map of the first species
TADI2.range = st_read("../USTreeAtlas/shp/taxodist/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(TADI2.range)  <- 4267 

# clip the species range based on US map
TADI2_clipped = st_intersection(TADI2.range, states.map)

# Add the species name back to the map info
TADI2_clipped$species = "Taxodium distichum"
all.data$species = "Taxodium distichum"

# clip all points outside of range with 50000 m buffer
TADI2_flag = cc_iucn(x = all.data, range = TADI2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
TADI2_occ_final = all.data[TADI2_flag, ]

# export the values
write.csv(TADI2_occ_final, file = "./Formatted.Data/Species.Nulls/Taxodium.distichum.null.csv")

#### Tilia americana ####
# read in the range map of the first species
TIAM.range = st_read("../USTreeAtlas/shp/tiliamer/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(TIAM.range)  <- 4267 

# clip the species range based on US map
TIAM_clipped = st_intersection(TIAM.range, states.map)

# Add the species name back to the map info
TIAM_clipped$species = "Tilia americana"
all.data$species = "Tilia americana"

# clip all points outside of range with 50000 m buffer
TIAM_flag = cc_iucn(x = all.data, range = TIAM_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
TIAM_occ_final = all.data[TIAM_flag, ]

# export the values
write.csv(TIAM_occ_final, file = "./Formatted.Data/Species.Nulls/Tilia.americana.null.csv")

#### Triadica sebifera ####
# clip the species range based on US map, from BIEN
(TRSE6.range.sf <- BIEN_ranges_load_species('Triadica sebifera'))

# The projection is geographic with latitude and longitude. 
# The datum is NAD27 which is EPSG:4267. The BIEN ranges are WG84.
TRSE6.range.2 = TRSE6.range.sf %>%
  st_transform(st_crs(states.map))

TRSE6.range.3 = terra::vect(TRSE6.range.2)
states.map.2 = terra::vect(states.map)
TRSE6.range.4 = terra::intersect(TRSE6.range.3, states.map.2)
TRSE6.range.5 = st_as_sf(TRSE6.range.4)

# Add the species name back to the map info
TRSE6.range.5$species = "Triadica sebifera"
all.data$species = "Triadica sebifera"

# clip all points outside of range with 50000 m buffer
TRSE6_flag = cc_iucn(x = all.data, range = TRSE6.range.5, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
TRSE6_occ_final = all.data[TRSE6_flag, ]

# export the values
write.csv(TRSE6_occ_final, file = "./Formatted.Data/Species.Nulls/Triadica.sebifera.null.csv")

#### Tsuga canadensis ####
# read in the range map of the first species
TSCA.range = st_read("../USTreeAtlas/shp/tsugcana/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(TSCA.range)  <- 4267 

# clip the species range based on US map
TSCA_clipped = st_intersection(TSCA.range, states.map)

# Add the species name back to the map info
TSCA_clipped$species = "Tsuga canadensis"
all.data$species = "Tsuga canadensis"

# clip all points outside of range with 50000 m buffer
TSCA_flag = cc_iucn(x = all.data, range = TSCA_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
TSCA_occ_final = all.data[TSCA_flag, ]

# export the values
write.csv(TSCA_occ_final, file = "./Formatted.Data/Species.Nulls/Tsuga.canadensis.null.csv")

#### Tsuga caroliniana ####
# read in the range map of the first species
TSCA2.range = st_read("../USTreeAtlas/shp/tsugcaro/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(TSCA2.range)  <- 4267 

# clip the species range based on US map
TSCA2_clipped = st_intersection(TSCA2.range, states.map)

# Add the species name back to the map info
TSCA2_clipped$species = "Tsuga caroliniana"
all.data$species = "Tsuga caroliniana"

# clip all points outside of range with 50000 m buffer
TSCA2_flag = cc_iucn(x = all.data, range = TSCA2_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                     value = "flagged", buffer = 50000)
# remove the points
TSCA2_occ_final = all.data[TSCA2_flag, ]

# export the values
write.csv(TSCA2_occ_final, file = "./Formatted.Data/Species.Nulls/Tsuga.caroliniana.null.csv")

#### Ulmus alata ####
# read in the range map of the first species
ULAL.range = st_read("../USTreeAtlas/shp/ulmualat/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ULAL.range)  <- 4267 

# clip the species range based on US map
ULAL_clipped = st_intersection(ULAL.range, states.map)

# Add the species name back to the map info
ULAL_clipped$species = "Ulmus alata"
all.data$species = "Ulmus alata"

# clip all points outside of range with 50000 m buffer
ULAL_flag = cc_iucn(x = all.data, range = ULAL_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ULAL_occ_final = all.data[ULAL_flag, ]

# export the values
write.csv(ULAL_occ_final, file = "./Formatted.Data/Species.Nulls/Ulmus.alata.null.csv")

#### Ulmus americana ####
# read in the range map of the first species
ULAM.range = st_read("../USTreeAtlas/shp/ulmuamer/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ULAM.range)  <- 4267 

# clip the species range based on US map
ULAM_clipped = st_intersection(ULAM.range, states.map)

# Add the species name back to the map info
ULAM_clipped$species = "Ulmus americana"
all.data$species = "Ulmus americana"

# clip all points outside of range with 50000 m buffer
ULAM_flag = cc_iucn(x = all.data, range = ULAM_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ULAM_occ_final = all.data[ULAM_flag, ]

# export the values
write.csv(ULAM_occ_final, file = "./Formatted.Data/Species.Nulls/Ulmus.americana.null.csv")

#### Ulmus crassifolia ####
# read in the range map of the first species
ULCR.range = st_read("../USTreeAtlas/shp/ulmucras/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ULCR.range)  <- 4267 

# clip the species range based on US map
ULCR_clipped = st_intersection(ULCR.range, states.map)

# Add the species name back to the map info
ULCR_clipped$species = "Ulmus crassifolia"
all.data$species = "Ulmus crassifolia"

# clip all points outside of range with 50000 m buffer
ULCR_flag = cc_iucn(x = all.data, range = ULCR_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ULCR_occ_final = all.data[ULCR_flag, ]

# export the values
write.csv(ULCR_occ_final, file = "./Formatted.Data/Species.Nulls/Ulmus.crassifolia.null.csv")

#### Ulmus rubra ####
# read in the range map of the first species
ULRU.range = st_read("../USTreeAtlas/shp/ulmurubr/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ULRU.range)  <- 4267 

# clip the species range based on US map
ULRU_clipped = st_intersection(ULRU.range, states.map)

# Add the species name back to the map info
ULRU_clipped$species = "Ulmus rubra"
all.data$species = "Ulmus rubra"

# clip all points outside of range with 50000 m buffer
ULRU_flag = cc_iucn(x = all.data, range = ULRU_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ULRU_occ_final = all.data[ULRU_flag, ]

# export the values
write.csv(ULRU_occ_final, file = "./Formatted.Data/Species.Nulls/Ulmus.rubra.null.csv")

#### Ulmus thomasii ####
# read in the range map of the first species
ULTH.range = st_read("../USTreeAtlas/shp/ulmuthom/")

# The projection is geographic with latitude and longitude. The datum is NAD27 which is EPSG:4267. 
# The ellipsoid is Clarke 1866, and the transformation parameters.
st_crs(ULTH.range)  <- 4267 

# clip the species range based on US map
ULTH_clipped = st_intersection(ULTH.range, states.map)

# Add the species name back to the map info
ULTH_clipped$species = "Ulmus thomasii"
all.data$species = "Ulmus thomasii"

# clip all points outside of range with 50000 m buffer
ULTH_flag = cc_iucn(x = all.data, range = ULTH_clipped, lon = "decimalLongitude", lat = "decimalLatitude", 
                    value = "flagged", buffer = 50000)
# remove the points
ULTH_occ_final = all.data[ULTH_flag, ]

# export the values
write.csv(ULTH_occ_final, file = "./Formatted.Data/Species.Nulls/Ulmus.thomasii.null.csv")
