# load libraries

library(topoclimate.pred)

# Load the DEM
elev <- raster("USGS_13_n43w101_20220726.tif")
names(elev) <- "elevation"

# Get microclimate data
clim <- bioclimate(elev, include_inputs = TRUE)

# change raster layer names
names(clim) = c("high_temp_C","low_temp_C","moisture_mm","northness","eastness","windward_exposure","mTPI",
                "macro_bio1_mean_annual_temp_C","macro_bio12_total_annual_precip_mm","macro_bio5_max_temp_warm_month_C",
                "macro_bi06_min_temp_cold_month_C")

# Define output file name based on raster file name
output_file <- paste0("microclim_", basename(elev@file@name))

# Write out the data for this raster
writeRaster(clim, filename = output_file, options="INTERLEAVE=BAND", overwrite=TRUE)





