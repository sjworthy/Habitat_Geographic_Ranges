
# load libraries
install.packages("tidyverse")
library(tidyverse)
install.packages("devtools")
library(devtools)
install_github("matthewkling/topoclimate.pred")
library(topoclimate.pred)
install.packages("conflicted")
library(conflicted)

conflicts_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

source("functions.R")


# Load GBIF data
gbif = read.csv("final.gbif.data.csv", row.names = 1)

# Read in the raster files as a list
raster_files <- list.files(".", pattern = "*.tif", full.names = TRUE)

for(raster_file in raster_files){
  
  # Load the DEM
  elev <- raster(raster_file)
  names(elev) <- "elevation"
  
  # crop GBIF data to only points within extents of DEM
  gbif.crop = gbif %>%
    filter(decimalLatitude >= elev@extent@ymin & decimalLatitude <= elev@extent@ymax) %>%
    filter(decimalLongitude >= elev@extent@xmin & decimalLongitude <= elev@extent@xmax)
  
  # crop the DEM based on these data to make things go faster
  ymin = as.numeric(format(round(min(gbif.crop$decimalLatitude), digits = 2), nsmall = 2)) -0.01
  ymax = as.numeric(format(round(max(gbif.crop$decimalLatitude), digits = 2), nsmall = 2)) + 0.01
  
  xmin = as.numeric(format(round(min(gbif.crop$decimalLongitude), digits = 2), nsmall = 2)) -0.01
  xmax = as.numeric(format(round(max(gbif.crop$decimalLongitude), digits = 2), nsmall = 2)) + 0.01
  
  ext <- extent(xmin, xmax, ymin, ymax)
  et <- crop(elev, ext)
  
  # Get microclimate data
  clim <- bioclimate_2(et, include_inputs = TRUE)
  
  # extract data for trees in this DEM
  trees.lat.long = gbif.crop[,c(2:3)]
  tree.clim = raster::extract(clim, trees.lat.long)
  tree.elev = raster::extract(et, trees.lat.long)
  
  # change column names
  colnames(tree.clim) = c("high_temp_C","low_temp_C","moisture_mm","northness","eastness","windward_exposure","mTPI","slope","aspect",
                          "macro_bio1_mean_annual_temp_C","macro_bio12_total_annual_precip_mm","macro_bio5_max_temp_warm_month_C",
                          "macro_bi06_min_temp_cold_month_C")
  

  # combine with original data
  tree.clim.2 = cbind(gbif.crop,tree.clim,tree.elev)
  
  # Define output file name based on raster file name
  output_file <- paste0("micro_clim_", basename(raster_file), ".csv")
  
  # Write out the data for this raster
  write.csv(tree.clim.2, file = output_file)
  
  # clear objects and freeing up memory
  rm(elev, gbif.crop, et, clim, tree.clim)
  gc()
  
}


