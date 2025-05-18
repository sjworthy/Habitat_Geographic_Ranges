# Code to extract variables from downscaled DEMs for each GBIF data point

library(tidyverse)
library(topoclimate.pred)
library(sf)
library(data.table)

# Read in gbif file
setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
gbif = read.csv("./Formatted.Data/good.gbif.data.csv", row.names = 1)

# subset gbif to just lat/long columns for extraction

gbif.lat.long = gbif %>%
  dplyr::select(species,decimalLongitude,decimalLatitude)

# Convert to sf object
points_sf <- sf::st_as_sf(gbif.lat.long, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

setwd("/Volumes/Backup Plus/microclimate_DEMs")

# micro clim DEMs are automatically read in as WGS84. Going to convert to NAD 83 since this is what they should be.
# did a test if it matters and get the same extracted values with WGS84 and NAD83 crs.
#example.dem = brick("./microclim_USGS_13_n30w100_20211103_tile2.tif", crs = "+proj=longlat +datum=NAD83 +no_defs")
example.dem = raster("./USGS_13_n26w081_20231221_tile7.tif")

points_nad83 <- st_transform(points_sf, crs(example.dem))

gbif_sp <- as(points_nad83, "Spatial")

# Read in list of microclim DEMs

raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)

# Prepare an empty list to store results
extracted_values <- list()

# Loop through rasters
for(raster_file in raster_files){
  
  # Load the DEM
  micro_DEM <- brick(raster_file, crs = "+proj=longlat +datum=NAD83 +no_defs")
  
  gbif.lat.long.2 = gbif.lat.long

  # Extract raster values for the points within this raster's extent
  extracted_vals <- as.data.frame(raster::extract(micro_DEM, gbif_sp))
  
  all.dat = cbind(gbif.lat.long.2,extracted_vals)
  
  # get rid of rows with NA
  clean_df <- na.omit(all.dat)
  
  # Change column names to match
  colnames(clean_df) <- c("species","decimalLongitude","decimalLatitude","high_temp_C", "low_temp_C", 
                          "moisture_mm", "northness", "eastness", "windward_exposure", 
                    "mTPI", "slope", "aspect", "macro_bio1_mean_annual_temp_C", 
                    "macro_bio12_total_annual_precip_mm", "macro_bio5_max_temp_warm_month_C", 
                    "macro_bi06_min_temp_cold_month_C")
  
# Combine the point data and the extracted values
  extracted_values[[length(extracted_values) + 1]] <- clean_df
  }


# Combine all extracted data into a single dataframe. Need row numbers to help with the merge
extracted_data <- do.call(rbind, extracted_values)

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
# had to run the loop a few times because vector memory exhaust after ~1500 rasters
write.csv(extracted_data, file = "./Formatted.Data/microclim.output.10.csv")

# use when need to read each item in the list out separately
for (i in seq_along(extracted_vals.2)) {
  filename <- paste0("./Formatted.Data/microclim.output_", i, ".csv")
  write.csv(extracted_vals.2[[i]], file = filename, row.names = FALSE)
}

# reading back in and combinding

# List all the CSV files you just saved
file_list <- list.files(
  path = "Formatted.Data",
  pattern = "^microclim\\.output_\\d+\\.csv$",
  full.names = TRUE
)

# Read them all into a list of data frames
list_of_dfs <- lapply(file_list, read.csv)

# Combine them using data.table (fast and memory-friendly)
combined_df <- rbindlist(list_of_dfs, use.names = TRUE, fill = TRUE)
#write.csv(combined_df, file = "./Formatted.Data/microclim.output.6.csv")

#### Reading back in all outputs and merging together ####

file_list <- list.files(
  path = "Formatted.Data",
  pattern = "^microclim\\.output.\\d+\\.csv$",
  full.names = TRUE
)

# Read them all into a list of data frames
list_of_dfs <- lapply(file_list, read.csv)

combined_df <- rbindlist(list_of_dfs, use.names = TRUE, fill = TRUE)

combined.df.2 = combined_df %>%
  distinct(across(2:4), .keep_all = TRUE)

# write.csv(combined_df, file = "./Results/climate.data.output.1.csv")

#### Comparing combined_df with full data ####

# combined = 1247810
# full = 1336865
# 89055 missing

temp = as.data.frame(table(sort(gbif$species)))
temp$microclim = table(sort(combined.df.2$species))

write.csv(temp, file = "species.compare.csv")

### Split by species ####

species_list <- split(combined.df.2, combined.df.2$species)

# Loop through the list and save each dataframe
for(species in names(species_list)) {
  # Create a file name for each species
  file_name <- paste0("./Results/Species.Climate/", "clim_data_", species, ".csv")
  
  # Write the dataframe to a CSV file
  write.csv(species_list[[species]], file_name, row.names = FALSE)
}


species_list_2 <- split(gbif, gbif$species)

# Loop through the list and save each dataframe
for(species in names(species_list_2)) {
  # Create a file name for each species
  file_name <- paste0("./Results/Species.Climate/", "gbif_data_", species, ".csv")
  
  # Write the dataframe to a CSV file
  write.csv(species_list_2[[species]], file_name, row.names = FALSE)
}

#### Trying to get microclim data for missing occurrences ####

setwd("/Volumes/My Passport for Mac/entire.DEMs/")
raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)
# Load the DEM
micro_DEM <- raster(raster_files[532], crs = "+proj=longlat +datum=NAD83 +no_defs")
micro_DEM

ext = extent(-82.99, -82.96, 42.41, 42.65)
et <- crop(micro_DEM, ext)

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges/Scripts/topo_data")
# get microclim

terr <- terrain(et, c("slope", "aspect"))
projection(terr) = projection(et)

ne <- northeast(terr)
wind <- windex(ne, terr)
tpi <- mtpi(et)
macro <- macroclimate(et)
d <- stack(setNames(tpi, "tpi"), ne, setNames(wind, "wind"), macro) %>%
  rasterToPoints() %>% as.data.frame() %>% as_tibble() %>% mutate(id = 1:nrow(.))
md <- readRDS("model_metadata.rds")[2,]
deltas <- get_deltas(md, d, "samples_full.csv")
topo <- microclimate(md, d, deltas, macro)
topo <- stack(topo, ne, wind, tpi, terr,
              setNames(macro, paste0("macro_", names(macro))))


# Extract raster values for the points within this raster's extent
gbif.lat.long.2 = gbif.lat.long

# Extract raster values for the points within this raster's extent
extracted_vals <- as.data.frame(raster::extract(topo, gbif_sp))

all.dat = cbind(gbif.lat.long.2,extracted_vals)

# get rid of rows with NA
clean_df <- na.omit(all.dat)

# Change column names to match
colnames(clean_df) <- c("species","decimalLongitude","decimalLatitude","high_temp_C", "low_temp_C", 
                        "moisture_mm", "northness", "eastness", "windward_exposure", 
                        "mTPI", "slope", "aspect", "macro_bio1_mean_annual_temp_C", 
                        "macro_bio12_total_annual_precip_mm", "macro_bio5_max_temp_warm_month_C", 
                        "macro_bi06_min_temp_cold_month_C")
View(clean_df)

# Combine the point data and the extracted values
extracted_values[[length(extracted_values) + 1]] <- clean_df

# Combine all extracted data into a single dataframe. Need row numbers to help with the merge
extracted_data <- do.call(rbind, extracted_values)

extracted_data.2 = extracted_data %>%
  distinct(across(1:3), .keep_all = TRUE)

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
write.csv(extracted_data.2, file = "./Formatted.Data/microclim.output.21.csv")

setwd("/Volumes/My Passport for Mac")

output.name = basename(micro_DEM@file@name)

# Define output file name based on raster file name
output_file <- paste0("microclim_8_", output.name)

# Write out the data for this raster
writeRaster(topo, filename = output_file, options="INTERLEAVE=BAND", overwrite=TRUE)





### Plotting ####

clim.dat = read.csv("./Results/climate.data.output.1.csv")

pc.clim = prcomp(clim.dat[,c(5:7)], scale  = TRUE, center = TRUE)
summary(pc.clim)
pc.clim$rotation
biplot(pc.clim)

clim.dat$PC1 = pc.clim$x[,1]
clim.dat$PC2 = pc.clim$x[,2]

ggplot(clim.dat, aes(x = PC1, y = PC2, color = species))+
  geom_point(size = 0.2, show.legend = FALSE)
  

cor.test(clim.dat$PC1,clim.dat$decimalLatitude)
cor.test(clim.dat$PC2,clim.dat$decimalLatitude)
cor.test(clim.dat$PC1,clim.dat$decimalLongitude)
cor.test(clim.dat$PC2,clim.dat$decimalLongitude)

ggplot(clim.dat, aes(species, PC1))+
  geom_boxplot()

library(ggbiplot)

