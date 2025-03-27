#' Topoclimate estimates
#'
#' This function produces estimates of biologically effective microclimate given a user-suppled digital elevation model (DEM). It is based on a model that leverages North American tree species occurrences as microclimate indicators. Estimates provided are for the posterior mode (penalized maximum likelihood) for model parameters.
#'
#' @param dem RasterLayer representing elevation, in meters. A spatial resolution of 10-30 m is recommended. The dataset must be in the continental US or Canada. DEMs covering areas larger than a small landscape may encounter computational challenges.
#' @param include_inputs Whether to return intermediate inputs in addition to topoclimate result (default FALSE)
#'
#' @return A raster stack of topoclimate variables, including "high_temperature" (in deg C), "low_temperature" (in deg C), and "moisture" (in mm).
#' @references Kling, Baer, & Ackerly (2023), in review.
#' @export
bioclimate_2 <- function(dem, include_inputs = FALSE){
  
  message("Calculating biogially effective topoclimate\n",
          "(Note: you can safely ignore any 'attempt to apply non-function' error messages;\n",
          "see https://github.com/rspatial/terra/issues/30 for more info on these gremlins.)")
  
  # enforce lat-lon crs
  if(!isLonLat(dem)){
    message("... projecting DEM to lon-lat for compatibility ...")
    dem <- projectRaster(dem, crs = "+proj=longlat +datum=NAD83 +no_defs")
  }
  
  # prep predictor data
  terr <- terrain(dem, c("slope", "aspect"))
  ne <- northeast(terr)
  wind <- windex(ne, terr)
  tpi <- mtpi(dem)
  macro <- macroclimate(dem)
  d <- stack(setNames(tpi, "tpi"), ne, setNames(wind, "wind"), macro) %>%
    rasterToPoints() %>% as.data.frame() %>% as_tibble() %>%
    mutate(id = 1:nrow(.))
  
  # calculate topoclimate estimates
  md <- readRDS(topo_data("model_metadata.rds"))[2,]
  deltas <- get_deltas(md, d, topo_data("samples_full.csv"))
  topo <- microclimate(md, d, deltas, macro)
  
  if(include_inputs) topo <- stack(topo, ne, wind, tpi, terr,
                                   setNames(macro, paste0("macro_", names(macro))))
  return(topo)
}