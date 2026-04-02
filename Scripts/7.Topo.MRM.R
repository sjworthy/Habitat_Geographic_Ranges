# Multiple Regression on Distance Matrices
# relationship between topographic variables and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

#### Topographic Model ####

# read in topographic data
topo.data = read.csv("All.Final.Data.csv", row.names = 1)

# split into species
species.list = split(topo.data, topo.data$species)

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
  
  # creating topographic data
  topo.dat = species.data %>%
    dplyr::select(northness,eastness,mTPI,slope,elevation)
  
  # calculate gower distance for topographic data
  topo.dist = gowdis(topo.dat)

  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/10000 # convert from m to hectometer
  
  # Perform MRM
  MRM.topo = MRM(topo.dist ~ geo.dist.3)
  
  # Create a row with all the results for the current species
  clean_species_name = unique(species.data$species)
  
  # save the model
  saveRDS(MRM.topo, file = paste0("./Topo/MRM_topo_model_", clean_species_name, ".RDS"))
  
  # Plot for each species
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    TopoDist = as.vector(topo.dist)
  )
  
  # save the data for later plotting
  write.csv(df_plot, file = paste0("./Topo/MRM_Topo_data_", clean_species_name, ".csv"))
  
  plot_obj = ggplot(df_plot, aes(x = GeoDist, y = TopoDist)) +
    geom_point(shape = 21, fill = "grey", color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggtitle(paste("MRM Topo Plot", clean_species_name)) +
    xlab("Geographic Distance (hectometer)") +
    ylab("Topographic Distance") +
    ylim(0,1)+
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"))
  
  # Save the plot
  ggsave(filename = paste0("./Topo/MRM_Topo_plot_", clean_species_name, ".png"),
         plot = plot_obj,
         width = 5, height = 5)

}

#### Getting results from models ####

topo.results <- data.frame(
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
setwd("./Results/Models/Topo.models/")
rds_files <- list.files(path = ".", pattern = "\\.RDS$", full.names = TRUE)

for(i in rds_files) {
  
  MRM.topo = readRDS(i)
  
  # Extract from model
  Intercept = MRM.topo$coef[1,1]
  Intercept.p.value = MRM.topo$coef[1,2]
  Slope = MRM.topo$coef[2,1]
  Slope.p.value = MRM.topo$coef[2,2]
  R2 = MRM.topo$r.squared[1]
  R2.p.value = MRM.topo$r.squared[2]
  F.value = MRM.topo$F.test[[1]]
  F.test.p.value = MRM.topo$F.test[[2]]
  
  # Create a row with all the results for the current species
  clean_species_name <- sub("^.*/MRM_topo_model_(.*)\\.RDS$", "\\1", i)
  
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
  topo.results <- rbind(topo.results, species.results)
  
}

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
write.csv(topo.results, file = "./Results/topo.MRM.results.csv")

#### Parsing Results ####

topo.results = read.csv("./Results/topo.MRM.results.csv")

range(topo.results$Intercept)
# 0.008334895 0.202622090
range(topo.results$Slope)
# -0.0001720559  0.0033826190
range(topo.results$R2)
# 9.058681e-07 4.388925e-01

ggplot(topo.results, aes(x = Intercept))+
  geom_density()+
  geom_vline(xintercept = 0.07776802) + # global null value
  theme_classic()

ggplot(topo.results, aes(x = Slope))+
  geom_density()+
  geom_vline(xintercept = -6.118209e-05) + # global null value
  theme_classic()

ggplot(topo.results, aes(x = R2))+
  geom_density()+
  geom_vline(xintercept = 0.003036681) + # global null value
  theme_classic()

