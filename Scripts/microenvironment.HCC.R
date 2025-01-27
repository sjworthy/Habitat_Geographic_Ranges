
# load libraries
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(foreach)) install.packages("foreach")
if(!require(doParallel)) install.packages("doParallel")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("matthewkling/topoclimate.pred")

# setup parallel backend
num_cores = 8
cluster = makeCluster(num_cores)
registerDoParallel(cluster)

# Load GBIF data
gbif = read.csv("./Formatted.Data/final.gbif.data.csv", row.names = 1)

# Read in the raster files as a list
raster_files <- list.files("./Raw.Data/", pattern = "*.tif", full.names = TRUE)

# Use foreach to parallelize the for-loop
foreach(raster_file = raster_files, .packages = c("raster","tidyverse","bioclim")) %dopar% {

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
  clim <- bioclimate(et, include_inputs = TRUE)
  
  # extract data for trees in this DEM
  trees.lat.long = gbif.crop[,c(2:3)]
  tree.clim = raster::extract(clim, trees.lat.long)
  
  # change column names
  colnames(tree.clim) = c("high_temp_C","low_temp_C","moisture_mm","northness","eastness","windward_exposure","mTPI","macro_bio1_mean_annual_temp_C","macro_bio12_total_annual_precip_mm","macro_bio5_max_temp_warm_month_C","macro_bi06_min_temp_cold_month_C")
  
  # combine with original data
  tree.clim.2 = cbind(gbif.crop,tree.clim)
  
  # Define output file name based on raster file name
  output_file <- paste0("micro_clim_", basename(raster_file), ".csv")
  
  # Write out the data for this raster
  write.csv(tree.clim.2, file = output_file)
  
  # clear objects and freeing up memory
  rm(elev, gbif.crop, et, clim, tree.clim)
  gc()
  
}

stopCluster(cl = cluster)









