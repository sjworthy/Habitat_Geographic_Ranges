# Multiple Regression on Distance Matrices
# relationship between soil variables and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)
library(ggbiplot)
library(vegan)

#### Soil MRM ####
# read in soil data
# soil.data = read.csv("gbif.final.csv", row.names = 1)
soil.data = read.csv("./Formatted.Data/gbif.final.csv", row.names = 1)

# split into species
species.list = split(soil.data, soil.data$species)
species.list = species.list[c(7,8,48,100)]

for(species.name in names(species.list)){
  
  species.data <- species.list[[species.name]]
  
  # randomly sample 45000 points b/c that gives ~ 1 billion data points (max is 2.1, but have issues > 50,000)
  if (nrow(species.data) > 45000) {
    set.seed(13)
    species.data <- species.data[sample(nrow(species.data), 45000), ]
  }
  
  # creating spatial matrix
  spat.dat = as.matrix(species.data[,c("decimalLongitude","decimalLatitude")])
  
  # creating soil data
  soil.dat = species.data %>%
    dplyr::select(ph_d0_100,clay_d0_100,sand_d0_100,silt_d0_100,db_d0_100,ec_d0_100,texture_d0_100)

  # calculate gower distance for scaled microclimate data
  soil.dist = gowdis(soil.dat)

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
  
  write.csv(df_plot, file = paste0("./Soil/MRM_Soil_data_", clean_species_name, ".csv"))
  
    plot_obj = ggplot(df_plot, aes(x = GeoDist, y = SoilDist)) +
      geom_point(shape = 21, fill = "grey", color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      ggtitle(paste("MRM Soil Plot", clean_species_name)) +
      xlab("Geographic Distance (hectometer)") +
      ylab("Soil Distance") +
      ylim(0,1)+
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"))
    
    # Save the plot
    ggsave(filename = paste0("./Soil/MRM_Soil_plot_", clean_species_name, ".png"),
           plot = plot_obj,
           width = 5, height = 5)
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
# 0.000116258 0.008980342
range(soil.results$R2)
# 0.002657517 0.600395009

ggplot(soil.results, aes(x = Intercept))+
  geom_density()+
  geom_vline(xintercept = 0.2418721) +
  theme_classic()

ggplot(soil.results, aes(x = Slope))+
  geom_density()+
  geom_vline(xintercept = 0.0004634458) +
  theme_classic()

ggplot(soil.results, aes(x = R2))+
  geom_density()+
  geom_vline(xintercept = 0.05967644) +
  theme_classic()

#### PCA of slopes, intercepts, R2 from MRM soil models ####

soil.results = read.csv("./Results/soil.MRM.results.csv", row.names = 1)

cor.test(soil.results$Intercept, soil.results$Slope) # r = -0.44
cor.test(soil.results$Intercept, soil.results$R2) # r = -0.69
cor.test(soil.results$Slope, soil.results$R2) # r = 0.41

# PCA of slope, intercept, and R2
pc = prcomp(soil.results[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc)
# PC1 68%, PC2 90%, PC3 100%
pc$rotation
# PC1: positively associated with Intercept and negatively associated with slope and R2
# PC2: negatively associated with Slope
# PC3: positively associated with Intercept and R2

plot(pc$x[,1], soil.results$Intercept, type = "l")
plot(pc$x[,2], soil.results$Intercept, type = "l")
plot(pc$x[,3], soil.results$Intercept, type = "l")

plot(pc$x[,1], soil.results$Slope, type = "l")
plot(pc$x[,2], soil.results$Slope)
plot(pc$x[,3], soil.results$Slope, type = "l")

plot(pc$x[,1], soil.results$R2, type = "l")
plot(pc$x[,2], soil.results$R2, type = "l")
plot(pc$x[,3], soil.results$R2, type = "l")

soil.pca.plot = ggbiplot(pc, varname.adjust = 1.3, varname.size = 6, 
                         varname.color = "#D8511D", alpha = 0.3)+
  geom_point(alpha = 0.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.plot

# ggsave("./Plots/soil.MRM.pca.png", width = 5, height = 5)

soil.pca.plot.labels = ggbiplot(pc, varname.adjust = 1.1, varname.size = 6, varname.color = "blue",
                                labels = soil.results$species, labels.size = 3)+
  geom_point(size = 0.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.plot.labels

#ggsave("./Plots/soil.MRM.pca.labels.png", width = 5, height = 5)

# add the global values

soil.results[123,c(2,4,6)] = c(0.2418721,0.0004634458,0.05967644)
soil.results[123,1] = "Global"
soil.results$species.2 = NA
soil.results[123,10] = "Global"

# PCA of slope, intercept, and R2
pc = prcomp(soil.results[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc)
# PC1 68%, PC2 90%, PC3 100%
pc$rotation

soil.pca.plot.labels = ggbiplot(pc, varname.adjust = 1.1, varname.size = 6, varname.color = "blue",
                                labels = soil.results$species.2, labels.size = 3, 
                                alpha = 0.3)+
  geom_point(alpha = 0.3)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.plot.labels

#### PCA without Abies fraseri and F.albicans ####

soil.results.2 = soil.results %>%
  filter(!species %in% c("Abies fraseri"))

pc.sub = prcomp(soil.results.2[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc.sub)
# PC1 72%, PC2 90%, PC3 100%
pc.sub$rotation
# PC1: negatively associated with Intercept and positively associated with slope and R2
# PC2: negatively associated with Slope
# PC3: positively associated with Intercept and R2

soil.pca.sub.plot = ggbiplot(pc.sub, varname.adjust = 1.3, varname.size = 6, 
                         varname.color = "#D8511D", alpha = 0.3)+
  geom_point(alpha = 0.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.sub.plot

# removing Fraxinus albicans also

soil.results.3 = soil.results.2 %>%
  filter(!species %in% c("Fraxinus albicans"))

pc.sub = prcomp(soil.results.3[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc.sub)
# PC1 77%, PC2 90%, PC3 100%
pc.sub$rotation
# PC1: negatively associated with Intercept and positively associated with slope and R2
# PC2: positively associated with Slope
# PC3: positively associated with all three

soil.pca.sub.plot = ggbiplot(pc.sub, varname.adjust = 1.3, varname.size = 6, 
                             varname.color = "#D8511D", alpha = 0.3)+
  geom_point(alpha = 0.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlim(-4,4)
soil.pca.sub.plot

#### PC1 and PC3 ####

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

soil.results[110, c(2,4,6)] = c(0.23, 0.00004, 0.05)
soil.results[110,1] = "global"


#### Testing if soil values are significantly different from global null ####

# read in 999 nulls

nulls = read.csv("./Results/Soil.Global.Null/global.null.999.results.csv", row.names = 1)

intercept.null.means <- mean(nulls$intercepts)
intercept.nulls.sds <- sd(nulls$intercepts)
slope.null.means <- mean(nulls$slopes)
slope.nulls.sds <- sd(nulls$slopes)
R2.null.means <- mean(nulls$R2)
R2.nulls.sds <- sd(nulls$R2)

# read in the observed values

obs = read.csv("./Results/soil.MRM.results.csv", row.names = 1)

for(i in 1:nrow(obs)){
  obs.intercept = obs[i,2]
  obs.slope = obs[i,4]
  obs.R2 = obs[i,6]
  
  ses.intercet <- (obs.intercept - intercept.null.means) / intercept.nulls.sds
  obs[i,10] = ses.intercet
  ses.slope <- (obs.slope - slope.null.means) / slope.nulls.sds
  obs[i,11] = ses.slope
  ses.R2 <- (obs.R2 - R2.null.means) / R2.nulls.sds
  obs[i,12] = ses.R2
  
  rank.intercept = rank(c(obs.intercept,nulls$intercepts))[1]
  obs[i,13] = rank.intercept
  rank.slope = rank(c(obs.slope,nulls$slopes))[1]
  obs[i,14] = rank.slope
  rank.R2 = rank(c(obs.R2,nulls$R2))[1]
  obs[i,15] = rank.R2
  
  p.val.intercept = rank.intercept/1000
  obs[i,16] = p.val.intercept
  p.val.slope = rank.slope/1000
  obs[i,17] = p.val.slope
  p.val.R2 = rank.R2/1000
  obs[i,18] = p.val.R2
}

colnames(obs)[10:18] = c("SES.intercept","SES.slope","SES.R2",
                         "Rank.intercept","Rank.slope","Rank.R2",
                         "P.val.intercept","P.val.slope","P.val.R2")

#write.csv(obs, file = "./Results/soil.global.null.compare.results.csv")

#### Putting species into Categories: two tailed ####

obs = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope
shifting = obs %>%
  filter(P.val.intercept < 0.026 & P.val.slope > 0.974)
# Strong Shifter, low intercept, high slope, high R2
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.026 & P.val.slope > 0.974 & P.val.R2 > 0.974)
# Specialist, low intercept, low slope
specialist = obs %>%
  filter(P.val.intercept < 0.026 & P.val.slope < 0.026)
# Strong Specialist, low intercept, low slope, high R2
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.026 & P.val.slope < 0.026 & P.val.R2 > 0.974)
# Overdisperser, high intercept, low slope
overdisperser = obs %>%
  filter(P.val.intercept > 0.974 & P.val.slope < 0.026)
# Strong Overdisperser, high intercept, low slope, low R2
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.974 & P.val.slope < 0.026 & P.val.R2 < 0.026)
# Overdispersed shifters, high intercept, high slope
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.974 & P.val.slope > 0.974)
# Strong Overdispersed shifters, high intercept, high slope, high R2
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.974 & P.val.slope > 0.974 & P.val.R2 > 0.974)

# get remaining species
# should remove 77 species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 45 species left

# Uncategorized: either slope or intercept is significant, but not both
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.026 | P.val.intercept > 0.974 |
           P.val.slope < 0.026 | P.val.slope > 0.974)

# True generalists: non-significant, intercept, slope, R2
true.generalists = remain.sp %>%
  filter(P.val.intercept >= 0.026 & P.val.intercept <= 0.974 &
           P.val.slope >= 0.026 & P.val.slope <= 0.974 &
           P.val.R2 >= 0.026 & P.val.R2 <= 0.974)
# R2.generalist: non-significant intercpet, slope but significant R2
R2.generalists = remain.sp %>%
  filter((P.val.intercept >= 0.026 & P.val.intercept <= 0.974) &
           (P.val.slope >= 0.026 & P.val.slope <= 0.974) &
           (P.val.R2 < 0.026 | P.val.R2 > 0.974))

# Add categories to soil dataframe

obs$Category = dplyr::case_when(
  obs$species %in% specialist$species ~ "specialists",
  obs$species %in% shifting$species ~ "shifting",
  obs$species %in% overdisperser$species ~ "overdisperser",
  obs$species %in% overdisperse.shifter$species ~ "overdisper.shifter",
  obs$species %in% c(true.generalists$species,R2.generalists$species,No.Cat$species) ~ "generalists",
  TRUE ~ NA_character_
)

write.csv(obs, "./Results/soil.global.null.compare.results.csv")


#### Putting species into Categories: one tailed: USED THIS ####

obs = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95 & P.val.R2 > 0.95)
# Specialist, low intercept, low slope
specialist = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05)
# Strong Specialist, low intercept, low slope, high R2
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05 & P.val.R2 > 0.95)
# Overdisperser, high intercept, low slope
overdisperser = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05)
# Strong Overdisperser, high intercept, low slope, low R2
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 < 0.05)
# Overdispersed shifters, high intercept, high slope
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 36 species left

# Uncategorized: either slope or intercept is significant, but not both
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.05 | P.val.intercept > 0.95 |
           P.val.slope < 0.05 | P.val.slope > 0.95)

# True generalists: non-significant, intercept, slope, R2
true.generalists = remain.sp %>%
  filter(P.val.intercept >= 0.05 & P.val.intercept <= 0.95 &
           P.val.slope >= 0.05 & P.val.slope <= 0.95 &
           P.val.R2 >= 0.05 & P.val.R2 <= 0.95)
# R2.generalist: non-significant intercpet, slope but significant R2
R2.generalists = remain.sp %>%
  filter((P.val.intercept >= 0.05 & P.val.intercept <= 0.95) &
           (P.val.slope >= 0.05 & P.val.slope <= 0.95) &
           (P.val.R2 < 0.05 | P.val.R2 > 0.95))

# Add categories to soil dataframe
obs$Category = dplyr::case_when(
  obs$species %in% specialist$species ~ "specialists",
  obs$species %in% shifting$species ~ "shifting",
  obs$species %in% overdisperser$species ~ "overdisperser",
  obs$species %in% overdisperse.shifter$species ~ "overdisper.shifter",
  obs$species %in% c(true.generalists$species,R2.generalists$species,No.Cat$species) ~ "generalists",
  TRUE ~ NA_character_
)

# Add significance for strength to dataframe
obs$significant = dplyr::case_when(
  obs$species %in% c(strong.shifting$species,strong.overdisperser$species,
                     strong.specialist$species, strong.overdisperse.shifter$species) ~ "significant",
  TRUE ~ "non-significant"
)


write.csv(obs, "./Results/soil.global.null.compare.results.csv")


#### Putting species into Categories: one tailed R2 only ####

obs = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope
shifting = obs %>%
  filter(P.val.intercept < 0.026 & P.val.slope > 0.974)
# Strong Shifter, low intercept, high slope, high R2
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.026 & P.val.slope > 0.974 & P.val.R2 > 0.95)
# Specialist, low intercept, low slope
specialist = obs %>%
  filter(P.val.intercept < 0.026 & P.val.slope < 0.026)
# Strong Specialist, low intercept, low slope, high R2
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.026 & P.val.slope < 0.026 & P.val.R2 > 0.95)
# Overdisperser, high intercept, low slope
overdisperser = obs %>%
  filter(P.val.intercept > 0.974 & P.val.slope < 0.026)
# Strong Overdisperser, high intercept, low slope, low R2
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.974 & P.val.slope < 0.026 & P.val.R2 < 0.05)
# Overdispersed shifters, high intercept, high slope
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.974 & P.val.slope > 0.974)
# Strong Overdispersed shifters, high intercept, high slope, high R2
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.974 & P.val.slope > 0.974 & P.val.R2 > 0.95)

# get remaining species
# should remove 77 species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 45 species left

# Uncategorized: either slope or intercept is significant, but not both
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.026 | P.val.intercept > 0.974 |
           P.val.slope < 0.026 | P.val.slope > 0.974)

# True generalists: non-significant, intercept, slope, R2
true.generalists = remain.sp %>%
  filter(P.val.intercept >= 0.026 & P.val.intercept <= 0.974 &
           P.val.slope >= 0.026 & P.val.slope <= 0.974 &
           P.val.R2 >= 0.026 & P.val.R2 <= 0.974)
# R2.generalist: non-significant intercpet, slope but significant R2
R2.generalists = remain.sp %>%
  filter((P.val.intercept >= 0.026 & P.val.intercept <= 0.974) &
           (P.val.slope >= 0.026 & P.val.slope <= 0.974) &
           (P.val.R2 < 0.026 | P.val.R2 > 0.974))

### NMDS of slopes, intercepts, R2 from MRM soil models for ESA ####

# read in soil data with categories
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

soil[123,c(2,4,6)] = c(0.2418721,0.0004634458,0.05967644)
soil[123,1] = "Null"
soil[123,19] = "Null"
soil[123,20] = "significant"

# data for NMDS
soil.2 = soil[,c(2,4,6)]

soil.nmds = metaMDS(soil.2, distance = "bray")
soil.nmds
# stress = 0.05272845, this is good

# plotting
nmds_scores <- as.data.frame(scores(soil.nmds)$sites)
nmds_scores$Category <- soil$Category
nmds_scores$shape = soil$significant


# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
soil.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Category, fill = Category, shape = shape),
             size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("significant" = 16, "non-significant" = 1),
                     name = "Strong",
                     labels = c("significant" = "Significant", "non-significant" = "Non-significant")) +
  scale_fill_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Category",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Category",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  theme_classic() +
  stat_ellipse(aes(color = Category), linewidth = 1) 
  #theme(legend.position = "none")

soil.nmds

ggsave("./Plots/ESA.plots/soil.MNDS.ellipses.legend.png", width = 5, height = 5)

### PCA for ESA ####

soil.results = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# add the global values

soil.results[123,c(2,4,6)] = c(0.2418721,0.0004634458,0.05967644)
soil.results[123,1] = "Global"
soil.results[123,19] = "Global"


# adding which points were significantly different from the null
soil.results$sig.diff.null = ifelse(soil.results$Category == "", "not.sig", "sig")

# Create a color column: orange for Global, black for others
#soil.results$point.color <- ifelse(soil.results$species == "Lower", "#D8511D", "black")

# Assign shape: open (shape 1) for not.sig, closed (shape 16) for sig
soil.results$point.shape <- ifelse(soil.results$sig.diff.null == "not.sig", 1, 16)

soil.results$GlobalFlag = ifelse(soil.results$species == "Global", "Global", "Other")

soil.results$point.color <- dplyr::case_when(
  soil.results$Category == "Lower" ~ "blue",
  soil.results$Category == "Shifting" ~ "green",
  soil.results$Category == "Generalist" ~ "black",
  TRUE ~ "grey"  # fallback if needed
)

# PCA of slope, intercept, and R2
pc = prcomp(soil.results[,c(2,4,6)], center = TRUE, scale = TRUE)
summary(pc)
# PC1 68%, PC2 90%, PC3 100%
pc$rotation

# Base plot with ggbiplot, no points (alpha = 0 hides them)
soil.pca.plot <- ggbiplot(pc, varname.adjust = 1.3, varname.size = 6, 
                          varname.color = "#D8511D", alpha = 0) +
  geom_vline(xintercept = 0.6802697) +
  geom_hline(yintercept = 0.03648963) +
  theme_classic() +
  xlim(-4, 4)

# Extract ggbiplot-built data
built_plot <- ggplot_build(ggbiplot(pc, varname.adjust = 1.3, varname.size = 6, varname.color = "#D8511D"))

# Get PC1 and PC2 coordinates
point_data <- built_plot$data[[1]]
point_data$species <- soil.results$species
point_data$point.color <- soil.results$point.color

# Add GlobalFlag and sig.diff.null back using row order
point_data$GlobalFlag <- soil.results$species == "Global"
point_data$sig.diff.null <- ifelse(soil.results$Category == "", "not.sig", "sig")

# Re-plot manually on top of ggbiplot
soil.pca.plot.2 <- soil.pca.plot +
  geom_point(data = point_data,
             aes(x = x, y = y, 
                 color = factor(point.color),
                 shape = factor(sig.diff.null),
                 alpha = factor(GlobalFlag)),
             size = 3, stroke = 1) +
  scale_shape_manual(values = c("sig" = 16, "not.sig" = 1),
                     name = "Different from Null",
                     labels = c("sig" = "Yes", "not.sig" = "No")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  guides(alpha = "none")+
  theme(legend.position = "right") +
  xlab("Standardized PC1 (67.9%)") +
  ylab("Standardized PC2 (21.7%)")

soil.pca.plot.2

# ggsave("./Plots/ESA.plots/soil.MRM.pca.png", width = 6, height = 6)


# Define x range for the plot
x_vals <- seq(0, 200, length.out = 50)

# Expand the data for plotting
plot_data <- soil %>%
  mutate(id = row_number()) %>%
  tidyr::crossing(x = x_vals) %>%
  mutate(y = Intercept + 
           Slope * x)

plot_data_2 = plot_data %>%
  filter(!species %in% "Abies fraseri")

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, group = id, color = Category)) +
  geom_line() +
  theme_minimal()+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#7C873E",
    "overdisper.shifter" = "#6E687E",
    "Global" = "black"))
  #scale_linetype_manual(values = c("significant" = "solid", "non-significant" = "dashed"))

### Macro versus Microclimate figure ESA ####

library(topoclimate.pred)
setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges")

# using Niobrara as an example
elev <- raster("./Raw.Data/USGS_13_n43w101_20220726.tif")
names(elev) <- "elevation"

# Crop to focus on Niobrara plot  
ext <- extent(-100.03, -100.02, 42.78, 42.79)
ext <- extent(-100.03, -100.01, 42.77, 42.79)
et <- crop(elev, ext)

# Terrian map of cropped area
hillshade.crop <- hillShade(terrain(et, "slope"), terrain(et, "aspect"))
plot(hillshade.crop, col = colorRampPalette(c("black", "white"))(200), legend = F)

# Get microclimate for Niobrara. 
# Set include_inputs = TRUE to get macroclimate variables to output. 
# Took 10 minutes to run. Use functions.R that adapts bioclimate function.

setwd("/Users/samanthaworthy/Documents/GitHub/Habitat_Geographic_Ranges/Scripts/topo_data")
terr <- terrain(et, c("slope", "aspect"))
projection(terr) = projection(et)

ne <- northeast(terr)
wind <- windex(ne, terr)
tpi <- mtpi(et)
macro <- macroclimate(et)
d <- stack(setNames(tpi, "tpi"), ne, setNames(wind, "wind"), macro) %>%
  rasterToPoints() %>% as.data.frame() %>% as_tibble() %>% mutate(id = 1:nrow(.))
md <- readRDS("model_metadata.rds")[2,]
deltas <- get_deltas(md, d, "samples_full.csv")
topo <- microclimate(md, d, deltas, macro)
topo.2 <- stack(topo, ne, wind, tpi, terr,
              setNames(macro, paste0("macro_", names(macro))))


plot(topo.2$moisture, col = viridis::viridis_pal()(25))[1]

x_breaks <- pretty(range(topo.data$x), n = 3)

# Add -100.0200 to breaks
x_breaks <- c(-100.0275,-100.0250,-100.0225,-100.0200)


topo.data = as.data.frame(rasterToPoints(hillshade.crop))
topo = ggplot(topo.data, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradient(low = "black", high = "white")+
  scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
topo

ggsave("./Plots/ESA.plots/Niobrara.topo.png", width = 7, height = 5)

data <- as.data.frame(rasterToPoints(topo.2))

precip.micro = ggplot(data, aes(x = x, y = y, fill = moisture)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Annual Precipitation (mm)", direction = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude")
  #theme(legend.position = "none")
precip.micro

ggsave("./Plots/ESA.plots/Niobrara.micro.precip.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.micro.precip.legend.png", width = 7, height = 5)

precip.macro = ggplot(data, aes(x = x, y = y, fill = macro_bio12)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Annual Precipitation (mm)", direction = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
precip.macro

ggsave("./Plots/ESA.plots/Niobrara.macro.precip.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.macro.precip.legend.png", width = 7, height = 5)

high.temp.micro = ggplot(data, aes(x = x, y = y, fill = high_temp)) +
  geom_raster() +
  scale_fill_viridis_c(name = "High Temperature (°C)", option = "C") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
high.temp.micro

ggsave("./Plots/ESA.plots/Niobrara.micro.high.temp.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.micro.high.temp.legend.png", width = 7, height = 5)

high.temp.macro = ggplot(data, aes(x = x, y = y, fill = macro_bio5)) +
  geom_raster() +
  scale_fill_viridis_c(name = "High Temperature (°C)", option = "C") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") 
  #theme(legend.position = "none")
high.temp.macro

ggsave("./Plots/ESA.plots/Niobrara.macro.high.temp.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.macro.high.temp.legend.png", width = 7, height = 5)

cold.temp.micro = ggplot(data, aes(x = x, y = y, fill = low_temp)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Low Temperature (°C)", option = "C") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
cold.temp.micro

ggsave("./Plots/ESA.plots/Niobrara.micro.low.temp.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.micro.low.temp.legend.png", width = 7, height = 5)

cold.temp.macro = ggplot(data, aes(x = x, y = y, fill = macro_bio6)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Low Temperature (°C)", option = "C") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude")
  #theme(legend.position = "none")
cold.temp.macro

ggsave("./Plots/ESA.plots/Niobrara.macro.cold.temp.png", width = 7, height = 5)
ggsave("./Plots/ESA.plots/Niobrara.macro.cold.temp.legend.png", width = 7, height = 5)

#### Populus deltoides map ####

sub = all.data %>%
  filter(decimalLongitude < -100 & decimalLatitude > 42)

all.data = read.csv("./Formatted.Data/gbif.final.all.csv", row.names = 1)
sp = all.data %>%
  filter(species %in% ("Populus deltoides"))
sp = sp %>%
  filter(!is.na(macro_bio5_max_temp_warm_month_C))
sp = sp %>%
  filter(macro_bio5_max_temp_warm_month_C < 37)

points_sf <- sf::st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

library(maps)
maps::map(database = "state")
library(tigris)
us_states <- states(cb = TRUE)
continental_states <- us_states %>%
  filter(!NAME %in% (c("Alaska","American Samoa","Guam","Commonwealth of the Northern Mariana Islands","Hawaii","United States Virgin Islands",
                       "Puerto Rico")))
states.map = continental_states %>%
  st_as_sf %>%
  st_transform(st_crs(points_sf))

# so both plots have the same scale
temp_range <- range(c(points_sf$moisture_mm, points_sf$macro_bio12_total_annual_precip_mm), na.rm = TRUE)

dis.plot = ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = points_sf, aes(color = macro_bio5_max_temp_warm_month_C), size = 2) +
  scale_color_viridis(name = "High Temperature (°C)", option = "C") +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()
dis.plot

ggsave("./Plots/ESA.plots/P.deltoides.hightemp.png", width = 5, height = 3)

p2 = ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = points_sf, aes(color = macro_bio12_total_annual_precip_mm), size = 2) +
  scale_color_viridis(name = "Macro", option = "C", limits = temp_range) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic() +
  ggtitle("Temperature at 10m Scale")

library(cowplot)
plot_grid(p1,p2)
  
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

ggplot(sp, aes(high_temp_C)) +
  geom_density()

# Asimina triloba - moisture

sp = all.data %>%
  filter(species %in% ("Ulmus thomasii"))

points_sf <- sf::st_as_sf(sp, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

range(sp$high_temp_C, na.rm = TRUE)
range(sp$macro_bio5_max_temp_warm_month_C, na.rm = TRUE)
range(sp$low_temp_C, na.rm = TRUE)
range(sp$macro_bi06_min_temp_cold_month_C, na.rm = TRUE)
range(sp$moisture_mm, na.rm = TRUE)
range(sp$macro_bio12_total_annual_precip_mm, na.rm = TRUE)

# soil plot

# Define bounding box
lon_min <- -100.5
lon_max <- -100
lat_min <-  42.8
lat_max <-  43

# Subset soil.data
soil.sub <- soil.data[
  soil.data$decimalLongitude >= lon_min &
    soil.data$decimalLongitude <= lon_max &
    soil.data$decimalLatitude  >= lat_min &
    soil.data$decimalLatitude  <= lat_max, ]

soil.plot = ggplot(soil.data, aes(x = decimalLongitude, y = decimalLatitude, fill = ph_d0_100)) +
  geom_tile() +
  scale_fill_gradient(low = "black", high = "white")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")
soil.plot

#### Category plots ####

x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.25      
slope <- 0.0002       
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df.generalist <- data.frame(x = x_vals, y = y_vals)


# Shifter

set.seed(123)

# Simulate realistic data
x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.05      # Low intercept
slope <- 0.0009        # High slope (relative to small y-scale)
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df <- data.frame(x = x_vals, y = y_vals)

# Plot
shifter.plot = ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "#5495CF", size = 2) +
  #geom_smooth(data = df.generalist, aes(x = x, y = y), method = "lm", 
              #se = FALSE, color = "darkgray", size = 2) + 
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 20) +
  labs(
       x = "Geographic Distance (km)",
       y = "Soil Dissimilarity")
  #theme(axis.title.y = element_text(color = "#B4674E"))
shifter.plot

ggsave("./Plots/ESA.plots/soil.shifter.png", width = 5, height = 5)

# specialists

set.seed(123)

# Simulate realistic data
x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.05      # Still low
slope <- 0.00001        # Lower slope than before
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df <- data.frame(x = x_vals, y = y_vals)

# Plot
specialist.plot <- ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "#DB4743", size = 2) +
  #geom_smooth(data = df.generalist, aes(x = x, y = y), method = "lm", 
              #se = FALSE, color = "darkgray", size = 2) + 
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 20) +
  labs(
       x = "Geographic Distance (km)",
       y = "Microclimate Dissimilarity")
  #theme(axis.title.y = element_text(color = "#2E8289"))

specialist.plot

ggsave("./Plots/ESA.plots/specialist.png", width = 5, height = 5)


# overdisperser

set.seed(123)

# Simulate realistic data
x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.75      # High
slope <- 0.00001        # Lower slope than before
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df <- data.frame(x = x_vals, y = y_vals)

# Plot
overdisperser.plot <- ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "#548F01", size = 2) +
  #geom_smooth(data = df.generalist, aes(x = x, y = y), method = "lm", 
              #se = FALSE, color = "darkgray", size = 2) + 
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 15) +
  labs(
       x = "Geographic Distance (km)",
       y = "Microclimate Dissimilarity")

overdisperser.plot

ggsave("./Plots/ESA.plots/overdisperser.png", width = 5, height = 5)


# overdispersed shifter

set.seed(123)

# Simulate realistic data
x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.75      # High
slope <- 0.0003        # high slope
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df <- data.frame(x = x_vals, y = y_vals)

# Plot
overdisperse.shifter.plot <- ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "#B46DB3", size = 2) +
  #geom_smooth(data = df.generalist, aes(x = x, y = y), method = "lm", 
              #se = FALSE, color = "darkgray", size = 2) + 
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 20) +
  labs( 
       x = "Geographic Distance (km)",
       y = "Microclimate Dissimilarity")

overdisperse.shifter.plot

ggsave("./Plots/ESA.plots/overdisperse.shifter.png", width = 5, height = 5)

# Generalist

set.seed(123)

# Simulate realistic data
x_vals <- seq(0, 1000, length.out = 100)
intercept <- 0.27      
slope <- 0.00020       
noise <- rnorm(100, mean = 0, sd = 0.02)

y_vals <- intercept + slope * x_vals + noise

df <- data.frame(x = x_vals, y = y_vals)

# Plot
generalist.plot <- ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "#F5AF4D", size = 55) +
  geom_smooth(data = df.generalist, aes(x = x, y = y), method = "lm", 
              se = FALSE, color = "darkgray", size = 2) + 
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 15) +
  labs(
       x = "Geographic Distance (km)",
       y = "Microclimate Dissimilarity")

generalist.plot

ggsave("./Plots/ESA.plots/generalist.null.png", width = 5, height = 5)

# all patterns on one plot

set.seed(123)

# Simulate realistic data
gen_x_vals <- seq(0, 1000, length.out = 100)
gen_intercept <- 0.27      
gen_slope <- 0.00024       
gen_noise <- rnorm(100, mean = 0, sd = 0.02)
gen_y_vals <- gen_intercept + gen_slope * gen_x_vals + gen_noise
gen_df <- data.frame(x = gen_x_vals, y = gen_y_vals)

# Simulate realistic data
ods_x_vals <- seq(0, 1000, length.out = 100)
ods_intercept <- 0.75      # High
ods_slope <- 0.0003        # high slope
ods_noise <- rnorm(100, mean = 0, sd = 0.02)
ods_y_vals <- ods_intercept + ods_slope * ods_x_vals + ods_noise
ods_df <- data.frame(x = ods_x_vals, y = ods_y_vals)

od_x_vals <- seq(0, 1000, length.out = 100)
od_intercept <- 0.75      # High
od_slope <- 0.00001        # Lower slope than before
od_noise <- rnorm(100, mean = 0, sd = 0.02)
od_y_vals <- od_intercept + od_slope * od_x_vals + od_noise
od_df <- data.frame(x = od_x_vals, y = od_y_vals)

spec_x_vals <- seq(0, 1000, length.out = 100)
spec_intercept <- 0.05      # Still low
spec_slope <- 0.00001        # Lower slope than before
spec_noise <- rnorm(100, mean = 0, sd = 0.02)
spec_y_vals <- spec_intercept + spec_slope * spec_x_vals + spec_noise
spec_df <- data.frame(x = spec_x_vals, y = spec_y_vals)

shift_x_vals <- seq(0, 1000, length.out = 100)
shift_intercept <- 0.05      # Low intercept
shift_slope <- 0.0009        # High slope (relative to small y-scale)
shift_noise <- rnorm(100, mean = 0, sd = 0.02)
shift_y_vals <- shift_intercept + shift_slope * shift_x_vals + shift_noise
shift_df <- data.frame(x = shift_x_vals, y = shift_y_vals)

# Plot
all.plot <- ggplot() +
  geom_smooth(data = gen_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#F5AF4D", size = 2) +
  #geom_smooth(data = ods_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#B46DB3", size = 2) +
  geom_smooth(data = od_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#548F01", size = 2)+ 
  geom_smooth(data = spec_df, aes(x = x, y = y),method = "lm", se = FALSE, color = "#DB4743", size = 2) +
  geom_smooth(data = shift_df, aes(x = x, y = y),method = "lm", se = FALSE, color = "#5495CF", size = 2) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 15) +
  labs(
       x = "Geographic Distance (km)",
       y = "Microclimate Dissimilarity")

all.plot

ggsave("./Plots/ESA.plots/all.patterns.png", width = 5, height = 5)


