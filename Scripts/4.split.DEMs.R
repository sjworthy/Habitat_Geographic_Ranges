## Script to split the DEMs in to 6 parts
# They are too large as a whole for downscaling

library(topoclimate.pred)
library(SpaDES)

setwd("/Volumes/My Passport for Mac/DEMS")

# Read in the raster files as a list
raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)

for(raster_file in raster_files){
  
  # Load the DEM
  elev <- raster(raster_file)
  names(elev) = str_remove(basename(elev@file@name), "\\.tif$")
  splitRaster(elev, nx = 3, ny = 3, path = "./split.rasters/")
}



