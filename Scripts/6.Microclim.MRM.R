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

# split into species
species.list = split(microclim.data, microclim.data$species)
# Round 1: complete
species.list = species.list[c(1,2,3,20,23,25,34,36,38,40,50,52,54,57,59,66,
                              69,70,75,78,88,91,93,94,104,105,106,
                              109,117)]
# Completed in Round 1:
# A.fraseri,A.barbatum,A.leucoderme,C.aquatica,
# C.illinoinensis,C.texana,F.albicans,F.caroliniana,F.profunda,G.aquatica,
# M.pomifera,M.fraseri,M.macrophylla,N.aquatica,N.ogeche,P.echinata,P.palustris,
# P.pungens,P.aquatica,P.heterophylla,Q.incana,Q.lyrata,Q.margaretta,
# Q.marilandica,Q.sinuata,Q.stellata,Q.texana,R.pseudoacacia,T.caroliniana

# Completed in Round 2: 
# A.glabra,B.nigra,C.alba,Q.imbricaria,

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
setwd("/Volumes/My Passport for Mac")
rds_files <- list.files(path = "./microclim", pattern = "\\.RDS$", full.names = TRUE)

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
# 0.01791455 0.16815540
range(microclim.results$Slope)
# 0.0006949658 0.0101925208
range(microclim.results$R2)
# 0.0255713 0.8708292

intercept.plot = ggplot(microclim.results, aes(x = Intercept))+
  geom_density()+
  geom_vline(xintercept = 0.07125867, size = 1.5) +
  theme_classic(base_size = 20) +
  ylab("Density")
intercept.plot

ggsave("./Plots/ESA.plots/ODS.intercept.png", height = 5, width = 5)

slope.plot = ggplot(microclim.results, aes(x = Slope))+
  geom_density()+
  geom_vline(xintercept = 0.001053699, size = 1.5) +
  theme_classic(base_size = 20) +
  ylab("Density")

slope.plot

ggsave("./Plots/ESA.plots/ODS.slope.png", height = 5, width = 5)

ggplot(microclim.results, aes(x = R2))+
  geom_density()+
  geom_vline(xintercept = 0.3762577, size = 1.5) +
  theme_classic(base_size = 20) +
  ylab("Density")+
  xlab("R-squared")

#### Testing if soil values are significantly different from global null ####

# read in 999 nulls

nulls = read.csv("./Results/Microclim.Global.Null/global.microclim.null.999.results.csv", row.names = 1)

intercept.null.means <- mean(nulls$intercepts)
intercept.nulls.sds <- sd(nulls$intercepts)
slope.null.means <- mean(nulls$slopes)
slope.nulls.sds <- sd(nulls$slopes)
R2.null.means <- mean(nulls$R2)
R2.nulls.sds <- sd(nulls$R2)

# read in the observed values

obs = read.csv("./Results/microclim.MRM.results.csv", row.names = 1)

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

write.csv(obs, file = "./Results/microclim.global.null.compare.results.csv")

#### Putting species into Categories: one tailed: USED THIS ####

obs = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

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
# 35 species left

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


write.csv(obs, "./Results/microclim.global.null.compare.results.csv")

### NMDS of slopes, intercepts, R2 from MRM soil models for ESA ####

# read in soil data with categories
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

microclim[123,c(2,4,6)] = c(0.07125867,0.001053699,0.3762577)
microclim[123,1] = "Null"
microclim[123,19] = "Null"
microclim[123,20] = "significant"

# data for NMDS
microclim.2 = microclim[,c(2,4,6)]

microclim.nmds = metaMDS(microclim.2, distance = "bray")
microclim.nmds
# stress = 0.03192277, this is good

# plotting
nmds_scores <- as.data.frame(scores(microclim.nmds)$sites)
nmds_scores$Category <- microclim$Category
nmds_scores$shape = microclim$significant

# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
microclim.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Category, fill = Category, shape = shape),
             size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("significant" = 16, "non-significant" = 1),
                     name = "Strength",
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

microclim.nmds

ggsave("./Plots/ESA.plots/microclim.MNDS.ellipses.legend.png", width = 5, height = 5)




