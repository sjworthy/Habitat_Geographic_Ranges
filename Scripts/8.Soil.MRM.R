# Multiple Regression on Distance Matrices
# relationship between soil variables and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)
#library(StatMatch)

# read in soil data
# soil.data = read.csv("gbif.final.csv", row.names = 1)
soil.data = read.csv("./Formatted.Data/gbif.final.csv", row.names = 1)

# split into species
species.list = split(soil.data, soil.data$species)
species.list = species.list[-c(7,100)]

for(species.name in names(species.list)){
  
  species.data <- species.list[[species.name]]
  
  # creating spatial matrix
  spat.dat = as.matrix(species.data[,c("decimalLongitude","decimalLatitude")])
  
  # creating soil data
  soil.dat = species.data %>%
    dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)
  #soil.dat.2 = soil.dat
  #soil.dat.2$texture_d0_100 = as.factor(soil.dat.2$texture_d0_100)
  
  # calculate gower distance for scaled microclimate data
  soil.dist = gowdis(soil.dat)
  #soil.dist.2 = gower.dist(soil.dat.2)
  #soil.dist.3 = as.dist(soil.dist.2)
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  geo.dist.3 = geo.dist.2/10000 # convert from m to hectometer
  
  # Perform MRM
  MRM.soil = MRM(soil.dist ~ geo.dist.3)
  
  # Create a row with all the results for the current species
  clean_species_name = unique(species.data$species)

  # save the model
  saveRDS(MRM.soil, file = paste0("./Soil/MRM_soil_model_", clean_species_name, ".RDS"))
  
  # Plot for each species
  df_plot <- data.frame(
    GeoDist = as.vector(geo.dist.3),
    SoilDist = as.vector(soil.dist)
  )
  
    plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoilDist)) +
      geom_point(shape = 21, fill = "grey", color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggtitle(paste("MRM Soil Plot", clean_species_name)) +
      xlab("Geographic Distance (hectometer)") +
      ylab("Soil Distance") +
      xlim(0,1)+
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"))
    
    # Save the plot
    ggsave(filename = paste0("./Soil/MRM_Soil_plot_", clean_species_name, ".png"),
           plot = plot_obj,
           width = 5, height = 5)
  
    write.csv(df_plot, file = paste0("./Soil/MRM_Soil_data_", clean_species_name, ".csv"))
    
    }


#### Getting results from models ####

soil.results <- data.frame(
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
setwd("/Volumes/My Passport for Mac")
rds_files <- list.files(path = "./soil", pattern = "\\.RDS$", full.names = TRUE)

for(i in rds_files) {
  
  MRM.soil = readRDS(i)
  
# Extract from model
Intercept = MRM.soil$coef[1,1]
Intercept.p.value = MRM.soil$coef[1,2]
Slope = MRM.soil$coef[2,1]
Slope.p.value = MRM.soil$coef[2,2]
R2 = MRM.soil$r.squared[1]
R2.p.value = MRM.soil$r.squared[2]
F.value = MRM.soil$F.test[[1]]
F.test.p.value = MRM.soil$F.test[[2]]

# Create a row with all the results for the current species
clean_species_name <- sub("^.*/MRM_soil_model_(.*)\\.RDS$", "\\1", i)

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
soil.results <- rbind(soil.results, species.results)

}

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")
# write.csv(soil.results, file = "./Results/soil.MRM.results.csv")

#### Parsing Results ####

range(soil.results$Intercept)
# 0.03188286 0.29324038
range(soil.results$Slope)
# 0.0000116258 0.0008980342
range(soil.results$R2)
# 0.002657517 0.510592952

ggplot(soil.results, aes(x = Intercept))+
  geom_density()+
  theme_classic()

ggplot(soil.results, aes(x = Slope))+
  geom_density()+
  theme_classic()

ggplot(soil.results, aes(x = R2))+
  geom_density()+
  theme_classic()

#### PCA of slopes, intercepts, R2 from MRM soil models ####

soil.results = read.csv("./Results/soil.MRM.results.csv", row.names = 1)

cor.test(soil.results$Intercept, soil.results$Slope) # r = -0.42
cor.test(soil.results$Intercept, soil.results$R2) # r = -0.60
cor.test(soil.results$Slope, soil.results$R2) # r = 0.42

soil.results[110, c(2,4,6)] = c(0.23, 0.00004, 0.05)
soil.results[110,1] = "global"

# PCA of slope, intercept, and R2
pc = prcomp(soil.results[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc)
# PC1 66%, PC2 87%, PC3 100%
pc$rotation
# PC1: positively associated with Intercept and negatively associated with slope and R2
# PC2: negatively associated with Slope
# PC3: negatively associated with Intercept and R2


plot(pc$x[,1], soil.results$Intercept, type = "l")
plot(pc$x[,2], soil.results$Intercept, type = "l")
plot(pc$x[,3], soil.results$Intercept, type = "l")

plot(pc$x[,1], soil.results$Slope, type = "l")
plot(pc$x[,2], soil.results$Slope, type = "l")
plot(pc$x[,3], soil.results$Slope, type = "l")

plot(pc$x[,1], soil.results$R2, type = "l")
plot(pc$x[,2], soil.results$R2, type = "l")
plot(pc$x[,3], soil.results$R2, type = "l")

library(ggbiplot)

soil.pca.plot = ggbiplot(pc, varname.adjust = 1.1, varname.size = 6, varname.color = "blue")+
  geom_point()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.plot

ggsave("./Plots/soil.MRM.pca.png", width = 5, height = 5)

soil.results$species.2 = NA
soil.results[110,10] = "GLOBAL"

soil.pca.plot.labels = ggbiplot(pc, varname.adjust = 1.1, varname.size = 6, varname.color = "blue",
                                labels = soil.results$species.2, labels.size = 3)+
  geom_point(size = 0.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.plot.labels

ggsave("./Plots/soil.MRM.pca.labels.png", width = 5, height = 5)

## PC1 and PC3

soil.pca.2.3.plot = ggbiplot(pc,choices = c(2,3),varname.adjust = 1.1, varname.size = 6, varname.color = "blue")+
  geom_point()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.2.3.plot

#### PC of raw soil data ####

# PCA of slope, intercept, and R2
pc.soil = prcomp(soil.data[,c(4:9)], center = TRUE, scale = TRUE)
summary(pc.soil)
# PC1 43%, PC2 62%, PC3 80%, PC4 91%
pc.soil$rotation
# PC1: positively associated with clay, silt and negatively associated with sand
# PC2: negatively associated with ph and db
# PC3: positively associated with ec and ph

soil.data.pca = ggbiplot(pc.soil,choices = c(1,2),varname.adjust = 1.1, varname.size = 6, varname.color = "blue")+
  geom_point()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()
  #xlim(-4,4)
soil.data.pca

ggsave("./Plots/soil.data.pca.png", width = 5, height = 5)

test = lm(df_plot$SoilDist~df_plot$GeoDist)
test.2 = lm(combined_data$SoilDist ~ combined_data$GeoDist)


sub = soil.results %>%
  filter(Slope > 0.00004 | R2 > 0.05)

sub.2 = soil.results %>%
  filter(Slope > 0.00004 & R2 > 0.05)

par(mfrow = c(3,1))

hist(sub.2$Slope, nclass = 20)
hist(sub.2$Intercept, nclass = 20)
hist(sub.2$R2, nclass = 20)

