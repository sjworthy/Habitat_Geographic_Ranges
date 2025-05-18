
### Read in packages
library(soilDB)
library(dplyr)
library(sf)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

### Load Lat/Lon CSV 
gbif = read.csv("good.gbif.data.csv", row.names = 1)

### Rearrange data and rename columns
gbif = gbif[,c(1,3,2)]
colnames(gbif) = c("species","longitude","latitude")

### Create subset of 10 species
gbif_10 <- subset(gbif, species == "Abies fraseri" | species == "Acer barbatum" | species ==
                    "Acer leucoderme"  | species ==  "Acer negundo" | species ==
                    "Acer nigrum" | species == "Acer pensylvanicum" | species == 
                  "Acer rubrum" | species == "Acer saccharinum"| species ==
                    "Acer saccharum"  | species ==  "Aesculus glabra")


# Create output directory 
#dir.create("processed_species", showWarnings = FALSE)

# === Function to classify USDA texture ===
usda_texture <- function(clay, silt, sand) {
  if (any(is.na(c(clay, silt, sand)))) return(NA_character_)
  
  total <- clay + silt + sand
  if (!is.na(total) && total != 100) {
    clay <- clay / total * 100
    silt <- silt / total * 100
    sand <- sand / total * 100
  }
  
  if (clay >= 40 && silt <= 40 && sand <= 45) return("Clay")
  if (clay >= 35 && silt >= 40) return("Silty Clay")
  if (clay >= 35 && sand >= 45) return("Sandy Clay")
  if (clay >= 27 && clay < 40 && sand <= 20 && silt > 40) return("Silty Clay Loam")
  if (clay >= 27 && clay < 40 && sand > 20 && sand <= 45) return("Clay Loam")
  if (clay >= 27 && clay < 40 && sand > 45) return("Sandy Clay Loam")
  if (clay >= 20 && clay < 27 && silt >= 28 && silt < 50 && sand <= 52) return("Loam")
  if (clay >= 20 && clay < 27 && silt >= 50) return("Silt Loam")
  if (clay >= 20 && clay < 27 && sand > 52) return("Sandy Loam")
  if (clay < 20 && silt >= 80) return("Silt")
  if (clay < 20 && silt >= 50 && sand <= 50) return("Silt Loam")
  if (clay < 20 && silt < 50 && sand <= 52) return("Loam")
  if (clay < 20 && silt < 30 && sand > 52 && sand < 85) return("Sandy Loam")
  if (clay < 15 && sand >= 85) return("Loamy Sand")
  if (clay < 10 && silt < 10 && sand >= 90) return("Sand")
  return("Other")
}

########### Function to get soil variables for 4 depth intervals + texture
get_soil_properties_for_point <- function(point) {
  tryCatch({
    s <- suppressMessages(SDA_spatialQuery(point, what = "mukey"))
    if (is.null(s) || nrow(s) == 0) stop("No mukey found.")
    
    mukey <- s$mukey[1]
    
    query <- sprintf("
      SELECT hz.hzdept_r, hz.hzdepb_r, hz.ph1to1h2o_r, hz.claytotal_r, hz.sandtotal_r,
             hz.silttotal_r, hz.dbthirdbar_r, hz.ec_r
      FROM mapunit mu
      INNER JOIN component co ON mu.mukey = co.mukey
      INNER JOIN chorizon hz ON co.cokey = hz.cokey
      WHERE mu.mukey = '%s'
    ", mukey)
    
    soil_data <- suppressMessages(SDA_query(query))
    
    if (nrow(soil_data) == 0) stop("No soil data for mukey.")
    
    depth_intervals <- list(d0_20 = c(0, 20),d0_30 = c(0, 30),
      d0_50 = c(0, 50), d0_100 = c(0, 100))
    
    get_weighted_mean <- function(df, var, top, bottom) {
      df <- df %>%
        filter(hzdept_r < bottom, hzdepb_r > top) %>%
        mutate(adj_top = pmax(hzdept_r, top),
          adj_bottom = pmin(hzdepb_r, bottom),
          thickness = adj_bottom - adj_top) %>%
        filter(thickness > 0, !is.na(.data[[var]]))
      if (nrow(df) == 0) return(NA_real_)
      sum(df[[var]] * df$thickness, na.rm = TRUE) / sum(df$thickness, na.rm = TRUE)
    }
    
    results <- list()
    
    for (int in names(depth_intervals)) {
      bounds <- depth_intervals[[int]]
      top <- bounds[1]; bottom <- bounds[2]
      
      ph   <- get_weighted_mean(soil_data, "ph1to1h2o_r", top, bottom)
      clay <- get_weighted_mean(soil_data, "claytotal_r", top, bottom)
      sand <- get_weighted_mean(soil_data, "sandtotal_r", top, bottom)
      silt <- get_weighted_mean(soil_data, "silttotal_r", top, bottom)
      db   <- get_weighted_mean(soil_data, "dbthirdbar_r", top, bottom)
      ec   <- get_weighted_mean(soil_data, "ec_r", top, bottom)
      
      tex <- usda_texture(clay, silt, sand)
      
      results[[paste0("ph_", int)]]      <- ph
      results[[paste0("clay_", int)]]    <- clay
      results[[paste0("sand_", int)]]    <- sand
      results[[paste0("silt_", int)]]    <- silt
      results[[paste0("db_", int)]]      <- db
      results[[paste0("ec_", int)]]      <- ec
      results[[paste0("texture_", int)]] <- tex
    }
    
    return(as.list(results))
    
  }, error = function(e) {
    # Return named NA list with all expected names
    interval_names <- c("d0_20", "d0_30", "d0_50", "d0_100")
    vars <- c("ph", "clay", "sand", "silt", "db", "ec", "texture")
    all_names <- as.vector(sapply(interval_names, function(d) paste0(vars, "_", d)))
    na_list <- as.list(rep(NA, length(all_names)))
    names(na_list) <- all_names
    return(na_list)
  })
}

########### Wrapper function to process one species
process_species_soil <- function(species_name, df) {
  message("Processing: ", species_name)
  
  sub_df <- df %>%
    filter(str_trim(species) == str_trim(species_name)) %>%
    filter(!is.na(longitude), !is.na(latitude))
  
  if (nrow(sub_df) == 0) {
    message("No valid records for: ", species_name)
    return(NULL)
  }
  
  points_sf <- tryCatch({
    sf_obj <- st_as_sf(sub_df, coords = c("longitude", "latitude"), crs = 4326)
    sf_obj <- sf_obj[st_is_valid(sf_obj), ]
    if (nrow(sf_obj) == 0) stop("No valid geometries after conversion.")
    sf_obj
  }, error = function(e) {
    message("Failed to convert points to sf for ", species_name, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(points_sf) || nrow(points_sf) == 0) {
    message("Skipping ", species_name, " due to invalid geometry.")
    return(NULL)
  }
  
  soil_values <- map(1:nrow(points_sf), function(i) {
    get_soil_properties_for_point(points_sf[i, ])
  })
  
  soil_values_df <- bind_rows(soil_values)
  
  # Bind results and export
  species_soil <- bind_cols(sub_df, soil_values_df)
  species_file_name <- gsub(" ", "_", species_name)
  out_path <- file.path("processed_species", paste0(species_file_name, "_soil.csv"))
  write_csv(species_soil, out_path)
  
  message("Saved: ", species_file_name, "_soil.csv")
}


##### Run on subset of full dataset

# Clean species column and get unique names
gbif_10$species <- str_trim(gbif_10$species)
species_list <- unique(gbif_10$species)

# Loop over species
for (sp in species_list) {
  species_file_name <- gsub(" ", "_", sp)
  out_path <- file.path("processed_species", paste0(species_file_name, "_soil.csv"))
  
  if (!file.exists(out_path)) {
    process_species_soil(sp, gbif_10)
  } else {
    message("Already processed: ", sp)
  }
}





