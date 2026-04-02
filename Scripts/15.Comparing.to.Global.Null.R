# Script to test if microclim, topo, or soil values of intercept, slope, and R2
# for each species significantly different from global null estimates of the values.

library(tidyverse)

#### Microclim #####

# read in 999 nulls

nulls = read.csv("./Results/global.microclim.null.999.results.csv", row.names = 1)

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

# write.csv(obs, file = "./Results/microclim.global.null.compare.results.csv")

mean(obs$Intercept) # 0.08361417
mean(obs$Slope) # 0.001727047
mean(obs$R2) # 0.3914622

mean(nulls$intercepts) # 0.0601699
mean(nulls$slopes) #0.0009340788
mean(nulls$R2) # 0.3839845

#### Topography ####
# read in 999 nulls

nulls = read.csv("./Results/global.topo.null.999.results.csv", row.names = 1)

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

# write.csv(obs, file = "./Results/topo.global.null.compare.results.csv")

mean(obs$Intercept)
mean(obs$Slope)
mean(obs$R2)

mean(nulls$intercepts) # 0.07776802
mean(nulls$slopes) # -6.118209e-05
mean(nulls$R2) # 0.003036681

#### Soil ####
# read in 999 nulls

nulls = read.csv("./Results/global.soil.null.999.results.csv", row.names = 1)

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

# write.csv(obs, file = "./Results/soil.global.null.compare.results.csv")

mean(obs$Intercept) # 0.1948436
mean(obs$Slope) # 0.0009679467
mean(obs$R2) # 0.101624

mean(nulls$intercepts) # 0.2306483
mean(nulls$slopes) # 0.0004231967
mean(nulls$R2) # 0.05483943

#### Putting species into Microclimate Categories: one tailed ####

obs = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope, 19 sp.
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2, 19 sp.
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# Specialist, low intercept, low slope, 2 sp.
specialist = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05)
# Strong Specialist, low intercept, low slope, high R2, 2 sp.
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdisperser, high intercept, low slope, 15 sp.
overdisperser = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05)
# Strong Overdisperser, high intercept, low slope, high R2, 1 sp.
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdispersed shifters, high intercept, high slope, 68 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 23 sp.
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 18 species left

# Uncategorized: either slope or intercept is significant, but not both
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.05 | P.val.intercept > 0.95 |
           P.val.slope < 0.05 | P.val.slope > 0.95)
# 18 species

# True generalists: non-significant, intercept, slope, R2, 0 sp. 
true.generalists = remain.sp %>%
  filter(P.val.intercept >= 0.05 & P.val.intercept <= 0.95 &
           P.val.slope >= 0.05 & P.val.slope <= 0.95 &
           P.val.R2 >= 0.05 & P.val.R2 <= 0.95)

# R2.generalist: non-significant intercept, slope but significant R2, 0 sp. 
R2.generalists = remain.sp %>%
  filter((P.val.intercept >= 0.05 & P.val.intercept <= 0.95) &
           (P.val.slope >= 0.05 & P.val.slope <= 0.95) &
           (P.val.R2 < 0.05 | P.val.R2 > 0.95))

# Add categories to microclim dataframe
obs$Category = dplyr::case_when(
  obs$species %in% specialist$species ~ "specialists",
  obs$species %in% shifting$species ~ "shifting",
  obs$species %in% overdisperser$species ~ "overdisperser",
  obs$species %in% overdisperse.shifter$species ~ "overdisper.shifter",
  obs$species %in% c(true.generalists$species,R2.generalists$species,No.Cat$species) ~ "generalists",
  TRUE ~ NA_character_)

# Add significance for strength to dataframe
obs$significant = dplyr::case_when(
  obs$species %in% c(strong.shifting$species,strong.overdisperser$species,
                     strong.specialist$species, strong.overdisperse.shifter$species) ~ "significant",
  TRUE ~ "non-significant")

#write.csv(obs, "./Results/microclim.global.null.compare.results.csv")

#### Putting species into Topography Categories: one tailed ####

obs = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope, 60 sp.
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2, 49 sp. 
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# Specialist, low intercept, low slope, 0 sp. 
specialist = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05)
# Strong Specialist, low intercept, low slope, high R2, 0 sp.
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdisperser, high intercept, low slope, 3 sp.
overdisperser = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05)
# Strong Overdisperser, high intercept, low slope, high R2, 1 sp.
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdispersed shifters, high intercept, high slope, 35 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 21 sp.
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 24 species left

# Uncategorized: either slope or intercept is significant, but not both, 21 sp.
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.05 | P.val.intercept > 0.95 |
           P.val.slope < 0.05 | P.val.slope > 0.95)

# True generalists: non-significant, intercept, slope, R2, 3 sp.
true.generalists = remain.sp %>%
  filter(P.val.intercept >= 0.05 & P.val.intercept <= 0.95 &
           P.val.slope >= 0.05 & P.val.slope <= 0.95 &
           P.val.R2 >= 0.05 & P.val.R2 <= 0.95)

# R2.generalist: non-significant intercpet, slope but significant R2
R2.generalists = remain.sp %>%
  filter((P.val.intercept >= 0.05 & P.val.intercept <= 0.95) &
           (P.val.slope >= 0.05 & P.val.slope <= 0.95) &
           (P.val.R2 < 0.05 | P.val.R2 > 0.95))

# Add categories to topography dataframe
obs$Category = dplyr::case_when(
  obs$species %in% specialist$species ~ "specialists",
  obs$species %in% shifting$species ~ "shifting",
  obs$species %in% overdisperser$species ~ "overdisperser",
  obs$species %in% overdisperse.shifter$species ~ "overdisper.shifter",
  obs$species %in% c(true.generalists$species,R2.generalists$species,No.Cat$species) ~ "generalists",
  TRUE ~ NA_character_)

# Add significance for strength to dataframe
obs$significant = dplyr::case_when(
  obs$species %in% c(strong.shifting$species,strong.overdisperser$species,
                     strong.specialist$species, strong.overdisperse.shifter$species,
                     true.generalists$species) ~ "significant",
  TRUE ~ "non-significant")

#write.csv(obs, "./Results/topo.global.null.compare.results.csv")

#### Putting species into Soil Categories: one tailed ####

obs = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# Shifter, low intercept, high slope, 66 sp. 
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2, 47 sp.
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# Specialist, low intercept, low slope, 12 sp. 
specialist = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05)
# Strong Specialist, low intercept, low slope, high R2, 0 sp.
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdisperser, high intercept, low slope, 12 sp.
overdisperser = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05)
# Strong Overdisperser, high intercept, low slope, high R2, 0 sp.
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdispersed shifters, high intercept, high slope, 10 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 7 sp.
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 22 species left

# Uncategorized: either slope or intercept is significant, but not both, 20 sp.
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.05 | P.val.intercept > 0.95 |
           P.val.slope < 0.05 | P.val.slope > 0.95)

# True generalists: non-significant, intercept, slope, R2, 1 sp.
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
  TRUE ~ NA_character_)

# Add significance for strength to dataframe
obs$significant = dplyr::case_when(
  obs$species %in% c(strong.shifting$species,strong.overdisperser$species,
                     strong.specialist$species, strong.overdisperse.shifter$species,
                     true.generalists$species) ~ "significant",
  TRUE ~ "non-significant")

# write.csv(obs, "./Results/soil.global.null.compare.results.csv")


