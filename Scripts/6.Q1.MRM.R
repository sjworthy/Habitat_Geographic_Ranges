# Multiple Regression on Distance Matrices
# relationship between microclimate and geogrpahic distance

# https://github.com/csdambros/BioGeoAmazonia

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

species.list = list.files("./Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  Intercept = numeric(),
  Slope = numeric(),
  Slope.p.value = numeric(),
  R2 = numeric(),
  F.value = numeric(),
  F.test.p.value = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # creating spatial matrix
  spat.dat = as.matrix(species[,c(2,3)])
  
  # creating microclim data, scaled
  microclim.dat = species %>%
    dplyr::select(high_temp_C,low_temp_C,moisture_mm) %>%
    scale()
  
  # calculate gower distance for scaled microclimate data
  microclim.dist = gowdis(as.matrix(microclim.dat))
  
  # calculate Haversine distance for spatial data
  geo.dist = distm(spat.dat, fun = distHaversine)
  geo.dist.2 = as.dist(geo.dist) # convert to dist object
  
  # Perform MRM
  MRM.mod = MRM(microclim.dist ~ geo.dist.2)
  
  # Extract from model
  Intercept = MRM.mod$coef[1,1]
  Slope = MRM.mod$coef[2,1]
  Slope.p.value = MRM.mod$coef[2,2]
  R2 = MRM.mod$r.squared[1]
  F.value = MRM.mod$F.test[1]
  F.test.p.value = MRM.mod$F.test[2]
  
  # Create a row with all the results for the current species
  species_name <- gsub("^clim_data_|\\.csv$", "", basename(species_file)) %>%
    gsub(" ", "_", .)  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    n = number.individuals,
    Intercept = Intercept,
    Slope = Slope,
    Slope.p.value = Slope.p.value,
    R2 = R2,
    F.value = F.value,
    F.test.p.value = F.test.p.value,
    stringsAsFactors = FALSE)
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}

# write.csv(results, file = "./Results/Q1.MRM.results.csv")

#### Summarizing the Results ####

dat = read.csv("./Results/Q1.MRM.results.csv")

range(dat$R2)
# 0.00000239 0.79067390
range(dat$Slope)
# 4.48e-10 9.52e-07

ggplot(dat, aes(Slope))+
  geom_density()

# how to plot from Chatgpt
dat.2 = dat %>%
  dplyr::select(Intercept,Slope) %>%
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
  labs(x = "Geographic Distance", y = "Microclimate Distance") +
  theme_classic() +
  theme(legend.position = "none")  # Optional: hides the legend if you don't need it


range(dat$Intercept)
# 0.02150594 0.31811891

range(dat$Slope.p.value)
# 0.001 0.953
# 4 species non-significant p > 0.05

#### Looking at Generalists versus Specialists ####

dat = read.csv("./Results/Q1.MRM.results.csv")

range(dat$R2)
# 0.00000239 0.79067390
range(dat$Slope)
# 4.48e-10 9.52e-07
range(dat$Intercept)
# 0.02150594 0.31811891

quantile(dat$R2)
#0%        25%        50%        75%       100% 
#0.00000239 0.10298696 0.24938713 0.41764822 0.79067390 
quantile(dat$Slope)
#0%       25%       50%       75%      100% 
#4.480e-10 8.925e-08 1.340e-07 2.290e-07 9.520e-07 
quantile(dat$Intercept)
#         0%        25%        50%        75%       100% 
#0.02150594 0.08141044 0.10206210 0.13081068 0.31811891 

# Classify true generalists as < 25% R2 and > 75% Intercept
true.gen = dat %>%
  filter(R2 < 0.10298696) %>%
  filter(Intercept > 0.13081068)
# 15 species

# Classify true specialists as R2 > 75% and Intercept < 25%
true.special = dat %>%
  filter(R2 > 0.41764822) %>%
  filter(Intercept < 0.08141044)
# 19 species

# Classify combo Intercept > 75 % and R2 > 75%
combo = dat %>%
  filter(R2 > 0.41764822) %>%
  filter(Intercept > 0.13081068)
# 0 species

# distribution of R2 values and intercept values

ggplot(dat, aes(y = Intercept, x = R2))+
  geom_point()
# strong negative relationship, R2 increases as intercept gets lower

# PCA of slope and R2

pc = prcomp(dat[,c(5,7)], scale = TRUE, center = TRUE)
summary(pc)
pc$rotation

library(ggbiplot)
ggbiplot(pc, labels = dat$species, varname.adjust = 1)+
  theme_classic()

# adding intercept
pc.all = prcomp(dat[,c(4,5,7)], scale = TRUE, center = TRUE)
summary(pc.all)
pc.all$rotation

ggbiplot(pc.all, labels = dat$species, varname.adjust = 1)+
  theme_classic()

### Heat Map ####

all.dat = read.csv("./Results/climate.data.output.1.csv")

pc.all = prcomp(all.dat[,c()])


##### Scaling


