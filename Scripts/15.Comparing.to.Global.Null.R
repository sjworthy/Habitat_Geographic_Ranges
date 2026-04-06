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

obs = read.csv("./Results/microclim.MRM.results.csv", row.names = 1) %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

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

mean(obs$Intercept) # 0.08386907
mean(obs$Slope) # 0.001750449
mean(obs$R2) # 0.3921366

mean(nulls$intercepts) # 0.06071191
mean(nulls$slopes) # 0.0009420475
mean(nulls$R2) # 0.3821377

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

obs = read.csv("./Results/topo.MRM.results.csv", row.names = 1) %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

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

mean(obs$Intercept) # 0.07699623
mean(obs$Slope) # 0.0002556053
mean(obs$R2) # 0.03882695

mean(nulls$intercepts) # 0.07825661
mean(nulls$slopes) # -6.267296e-05
mean(nulls$R2) # 0.003126151

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

obs = read.csv("./Results/soil.MRM.results.csv", row.names = 1) %>% 
  filter(!species %in% c("Ailanthus altissima","Paulownia tomentosa","Triadica sebifera"))

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

mean(obs$Intercept) # 0.1941249
mean(obs$Slope) # 0.0009784194
mean(obs$R2) # 0.1022176

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

# Specialist, low intercept, low slope, 1 sp.
specialist = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05)
# Strong Specialist, low intercept, low slope, high R2, 2 sp.
strong.specialist = specialist %>%
  filter(P.val.intercept < 0.05 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdisperser, high intercept, low slope, 14 sp.
overdisperser = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05)
# Strong Overdisperser, high intercept, low slope, high R2, 0 sp.
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdispersed shifters, high intercept, high slope, 68 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 24 sp.
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 17 species left

# Uncategorized: either slope or intercept is significant, but not both
No.Cat = remain.sp %>%
  filter(P.val.intercept < 0.05 | P.val.intercept > 0.95 |
           P.val.slope < 0.05 | P.val.slope > 0.95)
# 17 species

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

# Shifter, low intercept, high slope, 58 sp.
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2, 47 sp. 
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
# Strong Overdisperser, high intercept, low slope, high R2, 0 sp.
strong.overdisperser = overdisperser %>%
  filter(P.val.intercept > 0.95 & P.val.slope < 0.05 & P.val.R2 > 0.95)

# Overdispersed shifters, high intercept, high slope, 34 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 20 sp.
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

# Shifter, low intercept, high slope, 64 sp. 
shifting = obs %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95)
# Strong Shifter, low intercept, high slope, high R2, 46 sp.
strong.shifting = shifting %>%
  filter(P.val.intercept < 0.05 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# Specialist, low intercept, low slope, 13 sp. 
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

# Overdispersed shifters, high intercept, high slope, 9 sp.
overdisperse.shifter = obs %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95)
# Strong Overdispersed shifters, high intercept, high slope, high R2, 6 sp.
strong.overdisperse.shifter = overdisperse.shifter %>%
  filter(P.val.intercept > 0.95 & P.val.slope > 0.95 & P.val.R2 > 0.95)

# get remaining species
remain.sp = obs %>%
  filter(!species %in% c(shifting$species, specialist$species, overdisperser$species, 
                         overdisperse.shifter$species))
# 21 species left

# Uncategorized: either slope or intercept is significant, but not both, 19 sp.
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

#write.csv(obs, "./Results/soil.global.null.compare.results.csv")


