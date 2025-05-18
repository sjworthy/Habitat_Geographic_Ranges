# Code to generate the biologically plausible null
# WITHIN species: Calculate geographic distance and microclimate distance
# Randomly sample 50 points WITHIN each species - preserves microclimate and geographic distance within species
# Combine samples across species
## can't combine distance matrices, but WITHIN species will be maintained after sampling since columns
## and rows of species will line up in the matrix, prior to distance calculations.
## Thus, will calculate distances WITHIN and among species to generate the null
# Run MRM - extract slopes, intercepts, R2 values
# Repeat 1000 times
# This null model does NOT break species-specific relationships between geographic and microclimate distance, 
# does NOT break biogeographic realms/distributions (i.e. southern distribution versus northern distribution), 
# accounts for dispersal limitation, but is constrained by species-specific adaptations to the environment.

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# read in complete dataset

all.data = read.csv("all.microclim.csv")
# Fraxinus albicans 34
# Quercus margaretta 93

# Initialize vector to store intercept and slope values
slopes <- numeric()
intercepts <- numeric()
R2 <- numeric()

# split all.data into species-specific dataframe

split.data = split(all.data, all.data$species)

# loop through each species dataframe, calculate microclimate and geographic distance
# randomly select 50 rows of each

# Prepare an empty list to geo_dist and microclim.dist object
geo.dist_values <- list()
microclim.dist_values <- list()

for (i in 1:1000) {
  for(j in length(split.data))
  sp = as.data.frame(split.data[[j]])
  
  # creating spatial matrix
  spat.dat = as.matrix(sp[,c("decimalLongitude", "decimalLatitude")])
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/1000 # convert to km
  
  # Sample 50 geo.dist rows/columns from each species
  geo.matrix = as.matrix(geo.dist.3)
  geo.rows <- sample(1:nrow(geo.matrix), size = 50, replace = TRUE)
  geo_samples <- geo.matrix[geo.rows, geo.rows] # subset the matrix
  
  # creating microclimate data
  microclim.dat = sp %>%
    dplyr::select(high_temp_C,low_temp_C,moisture_mm)
  
  # calculate gower distance for microclimate data
  microclim.dist = gowdis(as.matrix(microclim.dat))
  
  # Sample 50 microclim.dist rows/columns from each species
  microclim.matrix = as.matrix(microclim.dist)
  # sample the same rows for geographic distance and microclimate distance
  microclim_samples <- microclim.matrix[geo.rows, geo.rows] # subset the matrix

  # Combine the point data and the extracted values
  geo.dist_values[[length(geo.dist_values) + 1]] <- geo_samples
  microclim.dist_values[[length(microclim.dist_values) + 1]] <- microclim_samples
  
  for (i in seq_along(geo.dist_values)) {
    mat <- geo.dist_values[[i]]
    # Add a group prefix to row/col names to make them unique
    rownames(mat) <- colnames(mat) <- paste0("Group", i, "_", seq_len(nrow(mat)))
    geo.dist_values[[i]] <- mat
  }
  
  all_names <- unlist(lapply(geo.dist_values, rownames))
  
  for (mat in geo.dist_values) {
    idx <- rownames(mat)
    combined_matrix[idx, idx] <- mat
  }

  

  
  # Run the MRM model with the bootstrapped matrices
  model <- MRM(microclim.dist ~ geo.dist.3)
  
  # Extract the intercept and slope values
  intercept_value <- model$coef[1,1]
  slope_value <- model$coef[2,1]
  R2_value <- model$r.squared[1]
  
  intercepts[i] <- intercept_value
  slopes[i] <- slope_value
  R2[i] <- R2_value
  
  # Plot every 100th iteration
  if (i %% 100 == 0) {
    df_plot <- data.frame(
      GeoDist = as.vector(geo.dist.3),
      MicroclimDist = as.vector(microclim.dist)
    )
    
    plot_obj = ggplot(df_plot, aes(x = GeoDist, y = MicroclimDist)) +
      geom_point(shape = 21, fill = "grey", color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggtitle(paste("MRM Plot - Iteration", i)) +
      xlab("Geographic Distance (km)") +
      ylab("Microclimate Distance") +
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"))
    
    # Save the plot
    ggsave(filename = paste0("MRM_Global_Null_plot_", i, ".pdf"),
           plot = plot_obj,
           width = 5, height = 5)
  }
  
}

# combine intercept and slope vectors into a dataframe
null_output = as.data.frame(intercepts)
null_output$slopes = slopes
null_output$R2 = R2

write.csv(null_output, file = "global.null.results.csv")
