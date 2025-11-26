# Multiple Regression on Distance Matrices
# relationship between topographic variables and geographic distance
# code written to run on UNL HCC SWAN

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

#### Topographic Model ####

# read in topographic data
topo.data = read.csv("gbif.final.all.csv", row.names = 1)

# split into species
species.list = split(topo.data, topo.data$species)
species.list = species.list[c(1,2,3,20,23,25,34,36,38,40,50,52,54,57,59,66,
                              69,70,75,78,88,91,93,94,104,105,106,
                              109,117)]

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
  topo.dat = species.data %>%
    dplyr::select(northness,eastness,mTPI,slope,elevation)
  
  # calculate gower distance for scaled microclimate data
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
setwd("/Volumes/My Passport for Mac")
rds_files <- list.files(path = "./topo", pattern = "\\.RDS$", full.names = TRUE)

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

range(topo.results$Intercept)
# 0.008332238 0.202622090
range(topo.results$Slope)
# -0.0001677807  0.0033826190
range(topo.results$R2)
# 9.058681e-07 4.378792e-01

ggplot(topo.results, aes(x = Intercept))+
  geom_density()+
  geom_vline(xintercept = 0.09715561) +
  theme_classic()

ggplot(topo.results, aes(x = Slope))+
  geom_density()+
  geom_vline(xintercept = -8.065893e-05) +
  theme_classic()

ggplot(topo.results, aes(x = R2))+
  geom_density()+
  geom_vline(xintercept = 0.003881603) +
  theme_classic()

#### Testing if soil values are significantly different from global null ####

# read in 999 nulls

nulls = read.csv("./Results/Topo.Global.Null/global.topo.null.999.results.csv", row.names = 1)

intercept.null.means <- mean(nulls$intercepts)
intercept.nulls.sds <- sd(nulls$intercepts)
slope.null.means <- mean(nulls$slopes)
slope.nulls.sds <- sd(nulls$slopes)
R2.null.means <- mean(nulls$R2)
R2.nulls.sds <- sd(nulls$R2)

# read in the observed values

obs = read.csv("./Results/topo.MRM.results.csv", row.names = 1)

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

write.csv(obs, file = "./Results/topo.global.null.compare.results.csv")

#### Putting species into Categories: one tailed: USED THIS ####

obs = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

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
# 28 species left

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


write.csv(obs, "./Results/topo.global.null.compare.results.csv")

### NMDS of slopes, intercepts, R2 from MRM soil models for ESA ####

# read in soil data with categories
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

topo[123,c(2,4,6)] = c(0.09715561,-8.065893e-05,0.003881603)
topo[123,1] = "Null"
topo[123,19] = "Null"
topo[123,20] = "significant"

# data for NMDS
topo.2 = topo[,c(2,4,6)]

# some values are negative so transform to all positive
# Shift all values to be positive by adding the absolute minimum
topo.trans <- topo.2 + abs(min(topo.2))

topo.nmds = metaMDS(topo.trans, distance = "bray")
topo.nmds
# stress = 0.05708112, this is good

# plotting
nmds_scores <- as.data.frame(scores(topo.nmds)$sites)
nmds_scores$Category <- topo$Category
nmds_scores$shape = topo$significant

# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
topo.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
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

topo.nmds

ggsave("./Plots/ESA.plots/topo.MNDS.ellipses.legend.png", width = 5, height = 5)





