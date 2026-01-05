# Elevation Extraction

library(tidyverse)
library(raster)
library(sf)
library(data.table)

# extracting elevation from original DEMs

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")

# Read in gbif data
gbif = read.csv("./Formatted.Data/good.gbif.data.csv")

# subset gbif to just lat/long columns for extraction
gbif.lat.long = gbif %>%
  dplyr::select(species,decimalLongitude,decimalLatitude)

gbif.test = unique(gbif.lat.long)

# Convert to sf object
points_sf <- sf::st_as_sf(gbif.lat.long, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Transform to NAD83
points_nad83 <- st_transform(points_sf, crs = 4269)

# Read in list of microclim DEMs

setwd("/Volumes/Backup Plus/entire.DEMs/")
raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)

# Prepare an empty list to store results
extracted_values <- list()

# Loop through rasters
for(raster_file in raster_files){
  
  # Load the DEM
  DEM <- raster(raster_file)
  
  # Extract raster values for the points within this raster's extent
  extracted_vals <- as.data.frame(raster::extract(DEM, points_nad83))
  
  # Remove rows where any of the extracted values are NA (across all points)
  valid_rows <- !apply(extracted_vals, 1, function(x) any(is.na(x)))
  
  # Subset the valid rows
  valid_points <- points_nad83[valid_rows, ]
  valid_vals <- extracted_vals[valid_rows, ]
  
  all.dat = cbind(valid_points,valid_vals)
  
  # Combine the point data and the extracted values
  extracted_values[[length(extracted_values) + 1]] <- all.dat
}

# Combine all extracted data into a single dataframe. Need row numbers to help with the merge
extracted_data <- do.call(rbind, extracted_values)
colnames(extracted_data)[2] = "elevation"

# Convert NAD83 back to WGS84 (EPSG:4326)
extracted_data$points_wgs84 <- st_transform(extracted_data$geometry, crs = 4326)

coords <- st_coordinates(extracted_data$points_wgs84)
extracted_data$decimalLongitude = coords[,1]
extracted_data$decimalLatitude = coords[,2]

extracted_data.2 = extracted_data %>%
  dplyr::select(species, elevation, decimalLongitude, decimalLatitude) %>%
  sf::st_drop_geometry(.)

write.csv(extracted_data.2, file = "./Formatted.Data/gbif.elevation.csv")

### Merge has too many data points. They are spread across many species

all.dat = merge(gbif.lat.long,extracted_data.2)
all.data.2 = unique(all.dat)

gbif.sp = as.data.frame(table(gbif.lat.long$species))
all.sp = as.data.frame(table(all.data.2$species))

gbif_list = split(gbif.lat.long, gbif.lat.long$species)
all_list = split(all.data.2, all.data.2$species)

# Loop through the list and save each dataframe
for(species in names(gbif_list)) {
  # Create a file name for each species
  file_name <- paste0("./Formatted.Data/elev/", "gbif_data_", species, ".csv")
  
  # Write the dataframe to a CSV file
  write.csv(gbif_list[[species]], file_name, row.names = FALSE)
}

# Loop through the list and save each dataframe
for(species in names(all_list)) {
  # Create a file name for each species
  file_name <- paste0("./Formatted.Data/elev/", "all_data_", species, ".csv")
  
  # Write the dataframe to a CSV file
  write.csv(all_list[[species]], file_name, row.names = FALSE)
}

# No Elevation data for 
# Acer saccharum	-73.073796	45.015647
# Carya cordiformis	-73.073308	45.015681
# Carya cordiformis	-73.072726	45.015722

### Get all elevation data back together ####

# Read in and combine data across all the data frames

file_list <- list.files(
  path = "Formatted.Data/elev",
  pattern = "^all_data_.*\\.csv$",
  full.names = TRUE
)

# Read them all into a list of data frames
list_of_dfs <- lapply(file_list, read.csv)

combined_df <- rbindlist(list_of_dfs, use.names = TRUE, fill = TRUE)

table(is.na(combined_df$elevation)) # 3 are NA

elev.final = combined_df %>%
  drop_na()

#write.csv(elev.final, file = "./Formatted.Data/elev.final.csv")
