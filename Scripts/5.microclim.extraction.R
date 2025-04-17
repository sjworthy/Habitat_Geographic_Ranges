# Code to extract variables from downscaled DEMs for each GBIF data point

library(tidyverse)
library(topoclimate.pred)

# Read in gbif file
setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
gbif = read.csv("./Formatted.Data/good.gbif.data.csv", row.names = 1)

# subset gbif to just lat/long columns for extraction

gbif.lat.long = gbif %>%
  dplyr::select(species,decimalLongitude,decimalLatitude)

# Read in list of microclim DEMs

setwd("/Volumes/Backup Plus/microclimate_DEMs")
raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)

# Prepare an empty list to store results
extracted_values <- list()

# Loop through rasters
for(raster_file in raster_files){
  
  # Load the DEM
  micro_DEM <- brick(raster_file)

  # Extract raster values for the points within this raster's extent
  extracted_vals <- as.data.frame(raster::extract(micro_DEM, gbif.lat.long[c("decimalLongitude", "decimalLatitude")]))
  
  # Remove rows where any of the extracted values are NA (across all points)
  valid_rows <- !apply(extracted_vals, 1, function(x) any(is.na(x)))
  
  # Subset the valid rows
  valid_points <- gbif.lat.long[valid_rows, ]
  valid_vals <- extracted_vals[valid_rows, ]
  
  all.dat = cbind(valid_points,valid_vals)
  
# Combine the point data and the extracted values
  extracted_values[[length(extracted_values) + 1]] <- all.dat
  }


# Change column names to match
new_colnames <- c("species","decimalLongitude","decimalLatitude","high_temp_C", "low_temp_C", "moisture_mm", "northness", "eastness", "windward_exposure", 
                  "mTPI", "slope", "aspect", "macro_bio1_mean_annual_temp_C", 
                  "macro_bio12_total_annual_precip_mm", "macro_bio5_max_temp_warm_month_C", 
                  "macro_bi06_min_temp_cold_month_C")

# Apply the new column names to each data frame in the extracted_values list
extracted_values <- lapply(extracted_values, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

# Combine all extracted data into a single dataframe. Need row numbers to help with the merge
extracted_data <- do.call(rbind, extracted_values)


extracted_data.2 = extracted_data
extracted_data.3 = extracted_data
extracted_data.4 = extracted_data

final.merge = rbind(extracted_data,extracted_data.2,extracted_data.3,extracted_data.4)

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
write.csv(final.merge, file = "./Results/climate.data.output.1.csv")

### Split by species ####

species_list <- split(final.merge, final.merge$species)

# Loop through the list and save each dataframe
for(species in names(species_list)) {
  # Create a file name for each species
  file_name <- paste0("./Results/Species.Climate/", "clim_data_", species, ".csv")
  
  # Write the dataframe to a CSV file
  write.csv(species_list[[species]], file_name, row.names = FALSE)
}

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

