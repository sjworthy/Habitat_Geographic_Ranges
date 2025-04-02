# Multiple Regression on Distance Matrices
# Relationship between topography and macroclimate

# https://github.com/csdambros/BioGeoAmazonia

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

species.list = list.files("./Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  macro.corr = numeric(),
  geo.temp.corr = numeric(),
  geo.ppt.corr = numeric(),
  Intercept = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  Intercept.temp = numeric(),
  Slope.temp = numeric(),
  Slope.p.value.temp = numeric(),
  R2.temp = numeric(),
  F.value.temp = numeric(),
  F.test.p.value.temp = numeric(),
  Intercept.ppt = numeric(),
  Slope.ppt = numeric(),
  Slope.p.value.ppt = numeric(),
  R2.ppt = numeric(),
  F.value.ppt = numeric(),
  F.test.p.value.ppt = numeric(),
  Intercept.geo = numeric(),
  Slope.geo = numeric(),
  Slope.p.value.geo = numeric(),
  R2.geo = numeric(),
  F.value.geo = numeric(),
  F.test.p.value.geo = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # creating spatial matrix
  macroclim.dat = species %>%
    select(macro_bio1_mean_annual_temp_C,macro_bio12_total_annual_precip_mm,) %>%
    scale()
  
  macro.temp = macroclim.dat[,1]
  macro.ppt = macroclim.dat[,2]
  
  macro.corr = cor(macro.temp,macro.ppt)
  
  # creating microclim data, scaled
  topo.dat = species %>%
    select(slope,aspect,mTPI) %>%
    scale()
  
  # calculate gower distance for scaled microclimate data
  topo.dist = gowdis(as.matrix(topo.dat))
  
  # calculate gower distance for macroclimate data
  macroclim.dist = gowdis(as.matrix(macroclim.dat))
  macro.temp.dist = gowdis(as.matrix(macro.temp))
  macro.ppt.dist = gowdis(as.matrix(macro.ppt))
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  
  geo.temp.corr = cor(geo.dist.2,macro.temp.dist)
  geo.ppt.corr = cor(geo.dist.2,macro.ppt.dist)
  
  # Perform MRM
  MRM.mod = MRM(topo.dist ~ macroclim.dist)
  
  # Extract from model
  Intercept = MRM.mod$coef[1,1]
  Slope = MRM.mod$coef[2,1]
  Slope.p.value = MRM.mod$coef[2,2]
  R2 = MRM.mod$r.squared[1]
  F.value = MRM.mod$F.test[1]
  F.test.p.value = MRM.mod$F.test[2]
  
  # Perform MRM with temp alone
  MRM.mod.temp = MRM(topo.dist ~ macro.temp.dist)
  
  # Extract from model
  Intercept.temp = MRM.mod.temp$coef[1,1]
  Slope.temp = MRM.mod.temp$coef[2,1]
  Slope.p.value.temp = MRM.mod.temp$coef[2,2]
  R2.temp = MRM.mod.temp$r.squared[1]
  F.value.temp = MRM.mod.temp$F.test[1]
  F.test.p.value.temp = MRM.mod.temp$F.test[2]
  
  # Perform MRM with ppt alone
  MRM.mod.ppt = MRM(topo.dist ~ macro.ppt.dist)
  
  # Extract from model
  Intercept.ppt = MRM.mod.ppt$coef[1,1]
  Slope.ppt = MRM.mod.ppt$coef[2,1]
  Slope.p.value.ppt = MRM.mod.ppt$coef[2,2]
  R2.ppt = MRM.mod.ppt$r.squared[1]
  F.value.ppt = MRM.mod.ppt$F.test[1]
  F.test.p.value.ppt = MRM.mod.ppt$F.test[2]
  
  # Perform MRM with geo dist
  MRM.mod.geo = MRM(topo.dist ~ geo.dist.2)
  
  # Extract from model
  Intercept.geo = MRM.mod.geo$coef[1,1]
  Slope.geo = MRM.mod.geo$coef[2,1]
  Slope.p.value.geo = MRM.mod.geo$coef[2,2]
  R2.geo = MRM.mod.geo$r.squared[1]
  F.value.geo = MRM.mod.geo$F.test[1]
  F.test.p.value.geo = MRM.mod.geo$F.test[2]
  
  # Create a row with all the results for the current species
  species_name <- gsub("^clim_data_|\\.csv$", "", basename(species_file)) %>%
    gsub(" ", "_", .)  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    n = number.individuals,
    macro.corr = macro.corr,
    geo.temp.corr = geo.temp.corr,
    geo.ppt.corr = geo.ppt.corr,
    Intercept = Intercept,
    Slope = Slope,
    Slope.p.value = Slope.p.value,
    R2 = R2,
    F.value = F.value,
    F.test.p.value = F.test.p.value,
    Intercept.temp = Intercept.temp,
    Slope.temp = Slope.temp,
    Slope.p.value.temp = Slope.p.value.temp,
    R2.temp = R2.temp,
    F.value.temp = F.value.temp,
    F.test.p.value.temp = F.test.p.value.temp,
    Intercept.ppt = Intercept.ppt,
    Slope.ppt = Slope.ppt,
    Slope.p.value.ppt = Slope.p.value.ppt,
    R2.ppt = R2.ppt,
    F.value.ppt = F.value.ppt,
    F.test.p.value.ppt = F.test.p.value.ppt,
    Intercept.geo = Intercept.geo,
    Slope.geo = Slope.geo,
    Slope.p.value.geo = Slope.p.value.geo,
    R2.geo = R2.geo,
    F.value.geo = F.value.geo,
    F.test.p.value.geo = F.test.p.value.geo,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}

# write.csv(results, file = "./Results/Q3.MRM.results.csv")

#### Summarizing the Results ####

dat = read.csv("./Results/Q3.MRM.results.csv")

##### considering MAT and MAP as macroclimate variables #####
range(dat$R2)
# 4.693801e-06 5.014648e-02
range(dat$Slope)
# -0.1722752  0.1657351

range(dat$Intercept)
# 0.1225016 0.4593709

ggplot(dat, aes(Slope))+
  geom_density()

# how to plot from Chatgpt
dat.2 = dat %>%
  select(Intercept,Slope) %>%
  mutate(Line.id = 1:107)

# Set up x-axis values for plotting
x_vals <- seq(0, 0.5, length.out = 100)

# Create a function to get the y-values for each line
get_line_y <- function(Intercept, Slope, x_vals) {
  Intercept + Slope * x_vals
}

# Generate the data for all lines
lines_df <- do.call(rbind, lapply(1:107, function(i) {
  data.frame(
    Line.id = i,
    x = x_vals,
    y = get_line_y(dat.2$Intercept[i], dat.2$Slope[i], x_vals)
  )
}))

# Plot all the lines using ggplot2
ggplot(lines_df, aes(x = x, y = y, group = Line.id, color = factor(Line.id))) +
  geom_line() +
  labs(x = "Macroclimate Distance", y = "Topographic Distance") +
  theme_classic() +
  theme(legend.position = "none")  # Optional: hides the legend if you don't need it

range(dat$Slope.p.value)
# 0.001 0.978

dat %>%
  filter(Slope.p.value > 0.05) %>%
  tally()
# 33 species non-significant p > 0.05

# number of species with positive slope
dat %>%
  filter(Slope > 0) %>%
  tally()
# 88 species 

# number of species with positive and significant slope
dat %>%
  filter(Slope > 0) %>%
  filter(Slope.p.value < 0.05) %>%
  tally()
# 67 species 


##### considering MAT ONLY as macroclimate variables #####
range(dat$R2.temp)
# 1.911568e-06 4.012910e-02
range(dat$Slope.temp)
# -0.08082809  0.09245877

range(dat$Intercept.temp)
# 0.1250809 0.4999203

ggplot(dat, aes(Slope.temp))+
  geom_density()

# how to plot from Chatgpt
dat.2 = dat %>%
  select(Intercept.temp,Slope.temp) %>%
  mutate(Line.id = 1:107)

# Set up x-axis values for plotting
x_vals <- seq(0, 0.5, length.out = 100)

# Create a function to get the y-values for each line
get_line_y <- function(Intercept, Slope, x_vals) {
  Intercept + Slope * x_vals
}

# Generate the data for all lines
lines_df <- do.call(rbind, lapply(1:107, function(i) {
  data.frame(
    Line.id = i,
    x = x_vals,
    y = get_line_y(dat.2$Intercept.temp[i], dat.2$Slope.temp[i], x_vals)
  )
}))

# Plot all the lines using ggplot2
ggplot(lines_df, aes(x = x, y = y, group = Line.id, color = factor(Line.id))) +
  geom_line() +
  labs(x = "Temperature Distance", y = "Topographic Distance") +
  theme_classic() +
  theme(legend.position = "none")  # Optional: hides the legend if you don't need it

range(dat$Slope.p.value.temp)
# 0.001 0.955

dat %>%
  filter(Slope.p.value.temp > 0.05) %>%
  tally()
# 36 species non-significant p > 0.05

##### considering MAP ONLY as macroclimate variables #####
range(dat$R2.ppt)
# 5.176774e-08 6.843002e-02
range(dat$Slope.ppt)
# -0.1024382  0.2221649

range(dat$Intercept.ppt)
# 0.1267673 0.4468580

ggplot(dat, aes(Slope.ppt))+
  geom_density()

# how to plot from Chatgpt
dat.2 = dat %>%
  select(Intercept.ppt,Slope.ppt) %>%
  mutate(Line.id = 1:107)

# Set up x-axis values for plotting
x_vals <- seq(0, 0.5, length.out = 100)

# Create a function to get the y-values for each line
get_line_y <- function(Intercept, Slope, x_vals) {
  Intercept + Slope * x_vals
}

# Generate the data for all lines
lines_df <- do.call(rbind, lapply(1:107, function(i) {
  data.frame(
    Line.id = i,
    x = x_vals,
    y = get_line_y(dat.2$Intercept.ppt[i], dat.2$Slope.ppt[i], x_vals)
  )
}))

# Plot all the lines using ggplot2
ggplot(lines_df, aes(x = x, y = y, group = Line.id, color = factor(Line.id))) +
  geom_line() +
  labs(x = "Precipitation Distance", y = "Topographic Distance") +
  theme_classic() +
  theme(legend.position = "none")  # Optional: hides the legend if you don't need it

range(dat$Slope.p.value.ppt)
# 0.001 0.973

dat %>%
  filter(Slope.p.value.ppt > 0.05) %>%
  tally()
# 38 species non-significant p > 0.05


##### considering GEO ONLY as macroclimate variables #####
range(dat$R2.geo)
# 1.947437e-08 1.979446e-01
range(dat$Slope.geo)
# -5.699940e-07  2.564288e-07

range(dat$Intercept.geo)
# 0.1231646 0.6654824

ggplot(dat, aes(Slope.geo))+
  geom_density()

# how to plot from Chatgpt
dat.2 = dat %>%
  select(Intercept.geo,Slope.geo) %>%
  mutate(Line.id = 1:107)

# Set up x-axis values for plotting
x_vals <- seq(0, 0.7, length.out = 100)

# Create a function to get the y-values for each line
get_line_y <- function(Intercept, Slope, x_vals) {
  Intercept + Slope * x_vals
}

# Generate the data for all lines
lines_df <- do.call(rbind, lapply(1:107, function(i) {
  data.frame(
    Line.id = i,
    x = x_vals,
    y = get_line_y(dat.2$Intercept.geo[i], dat.2$Slope.geo[i], x_vals)
  )
}))

# Plot all the lines using ggplot2
ggplot(lines_df, aes(x = x, y = y, group = Line.id, color = factor(Line.id))) +
  geom_line() +
  labs(x = "Geographic Distance", y = "Topographic Distance") +
  theme_classic() +
  theme(legend.position = "none")  # Optional: hides the legend if you don't need it

range(dat$Slope.p.value.geo)
# 0.001 0.990

dat %>%
  filter(Slope.p.value.geo > 0.05) %>%
  tally()
# 37 species non-significant p > 0.05


