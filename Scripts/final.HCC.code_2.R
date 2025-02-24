# load libraries
library(topoclimate.pred)

# Functions
northeast <- function(terr){
  message("... computing northness & eastness ...")
  slope_layer = terr$slope
  aspect_layer = terr$aspect
  # Assign the CRS (projection) to each layer manually
  projection(slope_layer) <- projection(elev)
  projection(aspect_layer) <- projection(elev)
  northness <- sin(slope_layer) * cos(aspect_layer)
  eastness <- sin(slope_layer) * sin(aspect_layer)
  ne <- stack(northness, eastness)
  ne <- setNames(ne, c("northness", "eastness"))
  ne
}

windex <- function(ne, terr){
  message("... computing windward exposure ...")
  
  wind <- stack("./Scripts/topo_data/wind.tif")
  wind <- setNames(wind, c("wspeed", "wdir", "waniso", "wu", "wv"))
  wind <- projectRaster(wind, terr[[1]])
  
  anglediff <- function(x, y){
    z <- abs(x - y)
    z[z > pi] <- abs(2*pi - z[z > pi])
    z
  }
  windiff <- pi - anglediff(raster::atan2(wind$wv, wind$wu), # downwind direction,
                            raster::atan2(ne$northness, ne$eastness)) # downslope direction
  windward <- cos(windiff) * sin(terr$slope) * wind$waniso / 0.14896 # 0.14896 is mean FIA anisotropy
  setNames(windward, "windward_exposure")
}

mtpi <- function(dem){
  message("... computing elevational position ...")
  fun <- function(x){
    x <- na.omit(x)
    sum(x) / length(x)
  }
  fw <- function(radius){
    w <- focalWeight(dem, radius / 111.1, fillNA = T) # 111.1 km/deg latitude
    w / w
  }
  tpi <- (dem * 3 - focal(dem, fw(.5), fun, pad = T) -
            focal(dem, fw(.225), fun, pad = T) -
            focal(dem, fw(.1), fun, pad = T)) / 3
  setNames(tpi, "mTPI")
}

macroclimate <- function(dem, interpolation = "bilinear"){
  message("... preparing macroclimate ...")
  
  # load CHELSA rasters
  f <- list.files("./Scripts/topo_data/", pattern = "CHELSA", full.names = T)
  clim <- stack(f[!grepl("xml", f)])
  names(clim) <- paste0("bio", c(1, 12, 5, 6))
  
  # resample (over an expanded area to avoid edge effects)
  clim <- crop(clim, raster::extent(dem) * 3)
  climate <- crop(resample(clim, dem, method = interpolation), dem)
  climate
}

tensor_splines <- function(x, y, xbounds = NULL, ybounds = NULL,
                           knots = 2, degree = 3){
  k = seq(0, 1, length.out = knots)
  k = k[2:(length(k)-1)]
  if(knots == 2) k <- NULL
  
  if(is.null(xbounds)) xbounds <- range(x)
  if(is.null(ybounds)) xbounds <- range(y)
  B1 <- splines::bs(x, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
  B2 <- splines::bs(y, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
  B <- rep(1, nrow(B1))
  for(i in 1:ncol(B1)) for(j in 1:ncol(B2)) B <- cbind(B, B1[,i] * B2[,j])
  B <- B[,2:ncol(B)]
  B
}

get_deltas <- function(md, # model metadata
                       g, # predictor data
                       samples){ # path to cmdstanr output csv
  
  log10inc <- function(x){ # change in log10 precip
    y <- rnorm(1)
    z <- rnorm(1)
    ((y*10^(z+x)) - (y*10^z)) / (y*10^z)
  }
  
  vars <- md$vars[[1]]
  mod_vars <- md$mod_vars[[1]]
  topo_vars <- md$topo_vars[[1]]
  subplots <- md$subplots[[1]]
  d <- readRDS("./Scripts/topo_data/model_data_full.rds")
  
  dmd <- d$md
  dmv <- d$mv
  spid <- d$spid
  dmd$sp_id <- as.integer(factor(dmd$species))
  
  g <- g %>%
    dplyr::select(id, bio1, bio12) %>%
    mutate(bio1 = (bio1 - dmv$bio1_mean) / dmv$bio1_sd,
           bio12 = (log10(bio12+1) - dmv$bio12_mean) / dmv$bio12_sd)
  
  d <- read.csv(samples, comment.char = "#") %>%
    mutate(i = 1:nrow(.)) %>%
    gather(param, value, -i)
  
  topo_mods1 <- c("mean", mod_vars)
  params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                   topo = rep(topo_vars, each = length(topo_mods1)))
  
  deltas <- d %>%
    filter(str_detect(param, "delta")) %>%
    mutate(param = str_remove(param, "delta\\."),
           param = str_replace_all(param, "\\.", "_")) %>%
    separate(param, c("param", "var"))
  
  n_bases <- (md$s_knots + md$s_degree - 1) ^ 2
  
  deltas <- deltas %>%
    rename(clim_var = var, basis = param) %>%
    mutate(basis = as.integer(basis),
           topo_var = topo_vars[ceiling(basis / n_bases)],
           basis = basis %% n_bases,
           basis = ifelse(basis == 0, n_bases, basis),
           clim_var = vars[as.integer(clim_var)]) %>%
    dplyr::select(-i)
  
  message("... computing deltas ...")
  delt <- function(gs){
    
    # basis functions
    z <- tensor_splines(ecdf(dmd[[mod_vars[1]]])(gs[[mod_vars[1]]]),
                        ecdf(dmd[[mod_vars[2]]])(gs[[mod_vars[2]]]),
                        knots = md$s_knots, degree = md$s_degree,
                        xbounds = 0:1, ybounds = 0:1) %>%
      as.data.frame() %>% as_tibble() %>%
      mutate(id = gs$id) %>%
      gather(basis, basis_value, -id) %>%
      mutate(basis = as.integer(str_remove(basis, "V")))
    
    # deltas
    gg <- unique(deltas$basis) %>%
      map_dfr(function(b){
        expand_grid(filter(deltas, basis == b) %>% select(-basis),
                    filter(z, basis == b) %>% select(-basis)) %>%
          mutate(delta = value * basis_value,
                 basis = b)}) %>%
      group_by(clim_var, topo_var, id) %>% #, i) %>%
      summarize(delta = sum(delta), .groups = "drop") %>%
      left_join(gs, by = "id")
    
    # rescale to native units
    gg$delta_rscl <- NA
    for(v in topo_vars) gg$delta_rscl[gg$topo_var == v] <-
      gg$delta[gg$topo_var == v] / dmv[[paste0(v, "_sd")]]
    for(v in vars) gg$delta_rscl[gg$clim_var == v] <-
      gg$delta_rscl[gg$clim_var == v] * dmv[[paste0(v, "_sd")]]
    gg <- gg %>%
      mutate(delta_rscl = ifelse(clim_var == "bio12", log10inc(delta_rscl), delta_rscl),
             delta_rscl = ifelse(clim_var == "bio12", delta_rscl * 100, delta_rscl),
             delta_rscl = ifelse(topo_var == "tpi", delta_rscl * 100, delta_rscl),
             bio1_rscl = bio1 * dmv$bio1_sd + dmv$bio1_mean,
             bio12_rscl = bio12 * dmv$bio12_sd + dmv$bio12_mean,
             bio12_rscl = 10^bio12_rscl - 1)
    pb$tick()$print()
    return(gg)
  }
  
  nslice <- 100 # chunk dataset to stay within memory
  pb <- progress_estimated(nslice)
  gg <- g %>%
    mutate(slice = rep(1:nslice, each = ceiling(nrow(.)/nslice))[1:nrow(.)]) %>%
    split(.$slice) %>%
    map_dfr(delt)
  return(gg)
}

microclimate <- function(md, d, deltas, macro){
  message("... compiling bioclimate estimates ...")
  ddd <- readRDS("./Scripts/topo_data/model_data_full.rds")
  mv <- ddd$mv
  
  ddd <- d %>%
    select(id, x, y, northness, eastness, windward = wind, tpi, bio5, bio6, bio12) %>%
    
    # standardize
    mutate(bio5 = (bio5 - mv$bio5_mean) / mv$bio5_sd,
           bio6 = (bio6 - mv$bio6_mean) / mv$bio6_sd,
           bio12 = (log10(bio12 + 1) - mv$bio12_mean) / mv$bio12_sd,
           northness = (northness - mv$northness_mean) / mv$northness_sd,
           eastness = (eastness - mv$eastness_mean) / mv$eastness_sd,
           windward = (windward - mv$windward_mean) / mv$windward_sd,
           tpi = (tpi - mv$tpi_mean) / mv$tpi_sd) %>%
    
    # join to delta data
    gather(topo_var, topo_value, northness:tpi) %>%
    gather(clim_var, clim_value, bio5:bio12) %>%
    filter(is.finite(topo_value)) %>%
    left_join(deltas %>% select(id, clim_var, topo_var, delta),
              by = c("id", "topo_var", "clim_var"))
  
  dd <- ddd %>%
    
    # calculate microclimate
    group_by(id, x, y, clim_var) %>%
    summarize(microclim = clim_value[1] + sum(topo_value * delta), .groups = "drop") %>%
    spread(clim_var, microclim) %>%
    
    # destandardize
    mutate(h = bio5 * mv$bio5_sd + mv$bio5_mean,
           c = bio6 * mv$bio6_sd + mv$bio6_mean,
           m = 10 ^ (bio12 * mv$bio12_sd + mv$bio12_mean) - 1)
  
  d <- left_join(d, select(dd, id, h:m), by = "id") %>%
    select(h:m) %>%
    as.matrix()
  topo <- macro[[1:3]]
  topo[] <- d
  names(topo) <- c("high_temp", "low_temp", "moisture")
  topo
}

md <- readRDS("./Scripts/topo_data/model_metadata.rds")[2,]

# Read in the raster files as a list

setwd("/Volumes/My Passport for Mac/DEMS/split.rasters")

raster_files <- list.files("./job3/", pattern = "*.tif", full.names = TRUE) # direct toward particular folder

# Loop through rasters
for(raster_file in raster_files){
  
  # Load the DEM
  elev <- raster(raster_file)
  names(elev) <- "elevation"
  output.name = basename(elev@file@name)
  
  setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
  
  # prep predictor data
  terr <- terrain(elev, c("slope", "aspect"))
  projection(terr) = projection(elev)
  gc()
  ne <- northeast(terr)
  gc()
  wind <- windex(ne, terr)
  gc()
  tpi <- mtpi(elev)
  gc()
  macro <- macroclimate(elev)
  gc()
  d <- stack(setNames(tpi, "tpi"), ne, setNames(wind, "wind"), macro) %>%
    rasterToPoints() %>% as.data.frame() %>% as_tibble() %>%
    mutate(id = 1:nrow(.))
  
  # calculate topoclimate estimates
  deltas <- get_deltas(md, d, "./Scripts/topo_data/samples_full.csv")
  gc()
  topo <- microclimate(md, d, deltas, macro)
  rm(d)
  raster::removeTmpFiles(h = 0)
  gc()
  
  topo <- stack(topo, ne, wind, tpi, terr,
                setNames(macro, paste0("macro_", names(macro))))
  
  
  # Define output file name based on raster file name
  output_file <- paste0("./job3.output/microclim_", output.name)
  
  setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
  
  # Write out the data for this raster
  writeRaster(topo, filename = output_file, options="INTERLEAVE=BAND", overwrite=TRUE)
  
  # clear objects and freeing up memory
  rm(elev,terr,ne,wind,tpi,macro,deltas,topo)
  raster::removeTmpFiles(h = 0)
  gc()
  
}
