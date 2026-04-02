# Multiple Regression on Distance Matrices
# relationship between microclimate and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

#### Microclim MRM ####
# read in data
microclim.data = read.csv("All.Final.Data.csv", row.names = 1)
#microclim.data = read.csv("./Formatted.Data/All.Final.Data.csv", row.names = 1)


# split into species
species.list = split(microclim.data, microclim.data$species)

for(species.name in names(species.list)){
  
  species.data <- species.list[[species.name]]
  
  # randomly sample 48000 points (max vector length is 2.1, but have issues > 50,000)
  # this mean only Liquidambar styraciflua (n = 51647), Acer rubrum (n = 70995), 
  # and Quercus palustric (n = 90009) are randomly sampled
  
  if (nrow(species.data) > 48000) {
    set.seed(13)
    species.data <- species.data[sample(nrow(species.data), 48000), ]
  }
  
  # creating spatial matrix
  spat.dat = as.matrix(species.data[,c("decimalLongitude","decimalLatitude")])
  
  # creating microclimate data:
  # summer daytime maximum temperature (high_temp_C), 
  # winter nighttime minimum temperature (low_temp_C)
  # annual precipitation (moisture)
  
  microclimat.dat = species.data %>%
    dplyr::select(high_temp_C,low_temp_C,moisture_mm)
  
  # calculate gower distance for microclimate data
  microclim.dist = gowdis(microclimat.dat)
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/10000 # convert from m to hectometer, helps with small values of slope and intercept
  
  # Perform MRM
  MRM.microclim = MRM(microclim.dist ~ geo.dist.3)
  
  # Create a row with all the results for the current species
  clean_species_name = unique(species.data$species)
  
  # save the model
  saveRDS(MRM.microclim, file = paste0("./Microclim/MRM_microclim_model_", clean_species_name, ".RDS"))
  
  # Plot for each species
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    MicroClimDist = as.vector(microclim.dist)
  )
  
  # save the data for later plotting
  write.csv(df_plot, file = paste0("./Microclim/MRM_microclim_data_", clean_species_name, ".csv"))
  
  plot_obj = ggplot(df_plot, aes(x = GeoDist, y = MicroClimDist)) +
    geom_point(shape = 21, fill = "grey", color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggtitle(paste("MRM Microclim Plot", clean_species_name)) +
    xlab("Geographic Distance (hectometer)") +
    ylab("Microclimate Distance") +
    ylim(0,1)+
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"))
  
  # Save the plot
  ggsave(filename = paste0("./Microclim/MRM_Microclim_plot_", clean_species_name, ".png"),
         plot = plot_obj,
         width = 5, height = 5)
  
}

#### Getting results from models ####

microclim.results <- data.frame(
  species = character(),
  Intercept = numeric(),
  Intercept.p.value = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  R2.p.value = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  stringsAsFactors = FALSE
)

# List all .RDS files in your directory
setwd("./Results/Models/Microclim.models/")
rds_files <- list.files(path = ".", pattern = "\\.RDS$", full.names = TRUE)

for(i in rds_files) {
  
  MRM.microclim = readRDS(i)
  
  # Extract from model
  Intercept = MRM.microclim$coef[1,1]
  Intercept.p.value = MRM.microclim$coef[1,2]
  Slope = MRM.microclim$coef[2,1]
  Slope.p.value = MRM.microclim$coef[2,2]
  R2 = MRM.microclim$r.squared[1]
  R2.p.value = MRM.microclim$r.squared[2]
  F.value = MRM.microclim$F.test[[1]]
  F.test.p.value = MRM.microclim$F.test[[2]]
  
  # Create a row with all the results for the current species
  clean_species_name <- sub("^.*/MRM_microclim_model_(.*)\\.RDS$", "\\1", i)
  
  species.results <- data.frame(
    species = clean_species_name,
    Intercept = Intercept,
    Intercept.p.value = Intercept.p.value,
    Slope = Slope,
    Slope.p.value = Slope.p.value,
    R2 = R2,
    R2.p.value = R2.p.value,
    F.value = F.value,
    F.test.p.value = F.test.p.value,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  microclim.results <- rbind(microclim.results, species.results)
  
}

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
write.csv(microclim.results, file = "./Results/microclim.MRM.results.csv")

#### Parsing Results ####

microclim.results = read.csv("./Results/microclim.MRM.results.csv")

range(microclim.results$Intercept)
# 0.01594465 0.16815540
range(microclim.results$Slope)
# 0.0006953444 0.0101925208
range(microclim.results$R2)
# 0.0255713 0.8631044

intercept.plot = ggplot(microclim.results, aes(x = Intercept))+
  geom_density()+
  geom_vline(xintercept = 0.0601699, linewidth = 1.5) + # this is the global intercept
  theme_classic(base_size = 20) +
  ylab("Density")
intercept.plot

ggsave("./Plots/ESA.plots/ODS.intercept.png", height = 5, width = 5)

slope.plot = ggplot(microclim.results, aes(x = Slope))+
  geom_density()+
  geom_vline(xintercept = 0.0009340788, linewidth = 1.5) + # this is the global slope
  theme_classic(base_size = 20) +
  ylab("Density")

slope.plot

ggsave("./Plots/ESA.plots/ODS.slope.png", height = 5, width = 5)

ggplot(microclim.results, aes(x = R2))+
  geom_density()+
  geom_vline(xintercept = 0.3839845, size = 1.5) + # this is the global R2
  theme_classic(base_size = 20) +
  ylab("Density")+
  xlab("R-squared")

