# RDA models to address question 1

# microclimate response variables: high temp, low temp, moisture
# microclimate predictors: lat and long

# split into species dataframes
# scale and center responses variables since temp and moisture on very differenct scales

library(vegan)
library(tidyverse)
library(ggpubr)

#### RDA ####

# significance testing using anova.cca, automatically does 999 permutations.
# F statistic corresponds to an overall test of sign. of an RDA by comparing the computed model to
# a null model. This test is based on the null hypothesis that the strength of the linear relationship
# calculated by the R2 is not larger than the value that would be obtained for unrelated Y and X matrices
# of the same size. 

# constrained proportion: variance of Y explained by X
# unconstrained proportion: unexplained variance in Y

species.list = list.files("./Results/Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  full_r2 = numeric(),
  full_constrained_var = numeric(),
  full_unconstrained_var = numeric(),
  geo_cor = numeric(),
  lat_r2 = numeric(),
  lat_constrained_var = numeric(),
  lat_unconstrained_var = numeric(),
  long_r2 = numeric(),
  long_constrained_var = numeric(),
  long_unconstrained_var = numeric(),
  full_anova_p_value = numeric(),
  full_anova_lat_p_value = numeric(),
  full_anova_long_p_value = numeric(),
  lat_anova_p_value = numeric(),
  long_anova_p_value = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # create response matrix, scaling variables b/c different units
  response = species %>%
    select(high_temp_C,low_temp_C,moisture_mm) %>%
    scale()  
  
  # Full rda with both latitude and longitude, get adjusted R2, constrained and unconstrained variance
  full.rda = rda(response ~ decimalLatitude + decimalLongitude, species)
  full_r2 <- RsquareAdj(full.rda)$adj.r.squared
  full.rda.summary = summary(full.rda)
  full.constrained.var = full.rda.summary$constr.chi/full.rda.summary$tot.chi
  full.unconstrained.var = full.rda.summary$unconst.chi/full.rda.summary$tot.chi
  
  ad = adonis2(response ~ decimalLatitude + decimalLongitude, data = species)
  
  # significance of rda
  full.anova = anova.cca(full.rda)
  full_anova_p_value <- full.anova$`Pr(>F)`[1]
  full.anova.terms = anova.cca(full.rda, by = "terms")
  full_anova_lat_p_value = full.anova.terms$`Pr(>F)`[1]
  full_anova_long_p_value = full.anova.terms$`Pr(>F)`[2]
  
  # get correlation between lat and long for the species
  geo.cor = cor(species$decimalLatitude, species$decimalLongitude)
  
  # rda with latitude, get adjusted R2, constrained and unconstrained variance
  lat.rda = rda(response ~ decimalLatitude, species)
  lat_r2 <- RsquareAdj(lat.rda)$adj.r.squared
  lat.rda.summary = summary(lat.rda)
  lat.constrained.var = lat.rda.summary$constr.chi/lat.rda.summary$tot.chi
  lat.unconstrained.var = lat.rda.summary$unconst.chi/lat.rda.summary$tot.chi
  
  # significance of rda
  lat.anova = anova.cca(lat.rda)
  lat_anova_p_value <- lat.anova$`Pr(>F)`[1]
  
  # rda with longitude, get adjusted R2, constrained and unconstrained variance
  long.rda = rda(response ~ decimalLongitude, species)
  long_r2 <- RsquareAdj(long.rda)$adj.r.squared
  long.rda.summary = summary(long.rda)
  long.constrained.var = long.rda.summary$constr.chi/long.rda.summary$tot.chi
  long.unconstrained.var = long.rda.summary$unconst.chi/long.rda.summary$tot.chi
  
  # significance of rda
  long.anova = anova.cca(long.rda)
  long_anova_p_value <- long.anova$`Pr(>F)`[1]
  
  # Create a row with all the results for the current species
  species_name <- gsub("^clim_data_|\\.csv$", "", basename(species_file)) %>%
    gsub(" ", "_", .)  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    n = number.individuals,
    full_r2 = full_r2,
    full_constrained_var = full.constrained.var,
    full_unconstrained_var = full.unconstrained.var,
    geo_cor = geo.cor,
    lat_r2 = lat_r2,
    lat_constrained_var = lat.constrained.var,
    lat_unconstrained_var = lat.unconstrained.var,
    long_r2 = long_r2,
    long_constrained_var = long.constrained.var,
    long_unconstrained_var = long.unconstrained.var,
    full_anova_p_value = full_anova_p_value,
    full_anova_lat_p_value = full_anova_lat_p_value,
    full_anova_long_p_value = full_anova_long_p_value,
    lat_anova_p_value = lat_anova_p_value,
    long_anova_p_value = long_anova_p_value,
    stringsAsFactors = FALSE
  )
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}

write.csv(results, file = "./Results/Q1.RDA.results.csv")
# remove species that didn't have at least 100 individuals

results.2 = results %>%
  filter(n > 99)

# number of rows where geo.cor > 0.5
results.2 %>%
  filter(geo_cor > 0.5 | geo_cor < -0.5) %>%
  nrow()
# 40 out of 125

range(results.2$full_r2)
range(results.2$full_constrained_var)
quantile(results.2$full_constrained_var)
#0%       25%       50%       75%      100% 
#0.1373631 0.4387693 0.5442747 0.6039492 0.7930756 

# number of species with less than half variance explained
results.2 %>%
  filter(full_constrained_var < 0.50) %>%
  nrow()
# 55 out of 125 (44%)

# number of species with less than 30% variance explained
results.2 %>%
  filter(full_constrained_var < 0.30) %>%
  nrow()
# 13 out of 125 (10%)

range(results.2$lat_constrained_var)
range(results.2$long_constrained_var)

# number of species where latitude explains > 50% of variation
results.2 %>%
  filter(lat_constrained_var > 0.50) %>%
  nrow()
# 14 out of 125 (11%)

# number of species where longitude explains > 50% of variation
results.2 %>%
  filter(long_constrained_var > 0.50) %>%
  nrow()
# 2 out of 125 (1.6%)

# subset to models with lower geo.cor
results.3 = results.2 %>%
  filter(geo_cor < 0.5 & geo_cor > -0.5)
# 85 species

# number of species with less than half variance explained
results.3 %>%
  filter(full_constrained_var < 0.50) %>%
  nrow()
# 35 out of 85 (41%)

# number of species with less than 30% variance explained
results.3 %>%
  filter(full_constrained_var < 0.30) %>%
  nrow()
# 4 out of 85 (5%)

# number of species where latitude explains > 50% of variation
results.3 %>%
  filter(lat_constrained_var > 0.50) %>%
  nrow()
# 5 out of 85 (6%)

# number of species where longitude explains > 50% of variation
results.3 %>%
  filter(long_constrained_var > 0.50) %>%
  nrow()
# 0 out of 85 (0%)

#### PCA of microclimate data ####

species.list = list.files("./Results/Species.Climate/", pattern = "*.csv", full.names = TRUE)

species.1 = read.csv(species.list[[1]])

pca.sp1 = prcomp(species.1[,c(4:6)], scale  = TRUE, center = TRUE)
summary(pca.sp1)

biplot(pca.sp1)

species.1$micro.clim.PC1 = pca.sp1$x[,1]

species.1.mod = lm(micro.clim.PC1 ~ decimalLatitude + decimalLongitude, data = species.1)
summary(species.1.mod)

cor.test(species.1$decimalLatitude, species.1$decimalLongitude)

species.1.mod.lat = lm(micro.clim.PC1 ~ decimalLatitude, data = species.1)
species.1.mod.long = lm(micro.clim.PC1 ~ decimalLongitude, data = species.1)
summary(species.1.mod.lat)
summary(species.1.mod.long)

ggplot(species.1, aes(x = decimalLatitude, y = micro.clim.PC1)) +
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 2, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 1.5, aes(label = ..rr.label..), size = 3)+
  theme_classic()


ggplot(species.1, aes(x = decimalLongitude, y = micro.clim.PC1)) +
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_regline_equation(label.y = 2, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 1.5, aes(label = ..rr.label..), size = 3)+
  theme_classic()

#### PCA loop ####

species.list = list.files("./Results/Species.Climate/", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
  n = numeric(),
  pc1.var = numeric(),
  pc2.var = numeric(),
  pc1.high.temp.loading = numeric(),
  pc1.low.temp.loading = numeric(),
  pc1.moisture.loading = numeric(),
  pc2.high.temp.loading = numeric(),
  pc2.low.temp.loading = numeric(),
  pc2.moisture.loading = numeric(),
  pc1.full.mod.intercept = numeric(),
  pc1.full.mod.lat.slope = numeric(),
  pc1.full.mod.long.slope = numeric(),
  pc1.full.mod.lat.p.value = numeric(),
  pc1.full.mod.long.p.value = numeric(),
  pc1.full.mod.adj.R2 = numeric(),
  pc2.full.mod.intercept = numeric(),
  pc2.full.mod.lat.slope = numeric(),
  pc2.full.mod.long.slope = numeric(),
  pc2.full.mod.lat.p.value = numeric(),
  pc2.full.mod.long.p.value = numeric(),
  pc2.full.mod.adj.R2 = numeric(),
  pc1.lat.mod.intercept = numeric(),
  pc1.lat.mod.lat.slope = numeric(),
  pc1.lat.mod.lat.p.value = numeric(),
  pc1.lat.mod.adj.R2 = numeric(),
  pc2.lat.mod.intercept = numeric(),
  pc2.lat.mod.lat.slope = numeric(),
  pc2.lat.mod.lat.p.value = numeric(),
  pc2.lat.mod.adj.R2 = numeric(),
  pc1.long.mod.intercept = numeric(),
  pc1.long.mod.long.slope = numeric(),
  pc1.long.mod.long.p.value = numeric(),
  pc1.long.mod.adj.R2 = numeric(),
  pc2.long.mod.intercept = numeric(),
  pc2.long.mod.long.slope = numeric(),
  pc2.long.mod.long.p.value = numeric(),
  pc2.long.mod.adj.R2 = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  number.individuals = nrow(species)
  
  # performing pca on microclimate data
  microclim.pc = prcomp(species[,c(4:6)], scale  = TRUE, center = TRUE)
  pc.summary = summary(microclim.pc)
  pc1.var = pc.summary$importance[3,1]
  pc2.var = pc.summary$importance[3,2]
  pc1.high.temp.loading = microclim.pc$rotation[1,1]
  pc1.low.temp.loading = microclim.pc$rotation[2,1]
  pc1.moisture.loading = microclim.pc$rotation[3,1]
  pc2.high.temp.loading = microclim.pc$rotation[1,2]
  pc2.low.temp.loading = microclim.pc$rotation[2,2]
  pc2.moisture.loading = microclim.pc$rotation[3,2]
  
  # add pc scores to species dataframe
  species$micro.clim.PC1 = microclim.pc$x[,1]
  species$micro.clim.PC2 = microclim.pc$x[,2]
  
  # full model PC1
  pc1.full.mod = lm(micro.clim.PC1 ~ decimalLatitude + decimalLongitude, data = species)
  pc1.full.mod.sum = summary(pc1.full.mod)
  pc1.full.mod.intercept = pc1.full.mod.sum$coefficients[1,1]
  pc1.full.mod.lat.slope = pc1.full.mod.sum$coefficients[2,1]
  pc1.full.mod.long.slope = pc1.full.mod.sum$coefficients[3,1]
  pc1.full.mod.lat.p.value = pc1.full.mod.sum$coefficients[2,4]
  pc1.full.mod.long.p.value = pc1.full.mod.sum$coefficients[3,4]
  pc1.full.mod.adj.R2 = pc1.full.mod.sum$adj.r.squared
  
  # full model PC2
  pc2.full.mod = lm(micro.clim.PC2 ~ decimalLatitude + decimalLongitude, data = species)
  pc2.full.mod.sum = summary(pc2.full.mod)
  pc2.full.mod.intercept = pc2.full.mod.sum$coefficients[1,1]
  pc2.full.mod.lat.slope = pc2.full.mod.sum$coefficients[2,1]
  pc2.full.mod.long.slope = pc2.full.mod.sum$coefficients[3,1]
  pc2.full.mod.lat.p.value = pc2.full.mod.sum$coefficients[2,4]
  pc2.full.mod.long.p.value = pc2.full.mod.sum$coefficients[3,4]
  pc2.full.mod.adj.R2 = pc2.full.mod.sum$adj.r.squared
  
  # lat models
  pc1.lat.mod = lm(micro.clim.PC1 ~ decimalLatitude, data = species)
  pc1.lat.mod.sum = summary(pc1.lat.mod)
  pc1.lat.mod.intercept = pc1.lat.mod.sum$coefficients[1,1]
  pc1.lat.mod.lat.slope = pc1.lat.mod.sum$coefficients[2,1]
  pc1.lat.mod.lat.p.value = pc1.lat.mod.sum$coefficients[2,4]
  pc1.lat.mod.adj.R2 = pc1.lat.mod.sum$adj.r.squared
  
  # lat model PC2
  pc2.lat.mod = lm(micro.clim.PC2 ~ decimalLatitude, data = species)
  pc2.lat.mod.sum = summary(pc2.lat.mod)
  pc2.lat.mod.intercept = pc2.lat.mod.sum$coefficients[1,1]
  pc2.lat.mod.lat.slope = pc2.lat.mod.sum$coefficients[2,1]
  pc2.lat.mod.lat.p.value = pc2.lat.mod.sum$coefficients[2,4]
  pc2.lat.mod.adj.R2 = pc2.lat.mod.sum$adj.r.squared
  
  # long models
  pc1.long.mod = lm(micro.clim.PC1 ~ decimalLongitude, data = species)
  pc1.long.mod.sum = summary(pc1.long.mod)
  pc1.long.mod.intercept = pc1.long.mod.sum$coefficients[1,1]
  pc1.long.mod.long.slope = pc1.long.mod.sum$coefficients[2,1]
  pc1.long.mod.long.p.value = pc1.long.mod.sum$coefficients[2,4]
  pc1.long.mod.adj.R2 = pc1.long.mod.sum$adj.r.squared
  
  # long model PC2
  pc2.long.mod = lm(micro.clim.PC2 ~ decimalLongitude, data = species)
  pc2.long.mod.sum = summary(pc2.long.mod)
  pc2.long.mod.intercept = pc2.long.mod.sum$coefficients[1,1]
  pc2.long.mod.long.slope = pc2.long.mod.sum$coefficients[2,1]
  pc2.long.mod.long.p.value = pc2.long.mod.sum$coefficients[2,4]
  pc2.long.mod.adj.R2 = pc2.long.mod.sum$adj.r.squared
  
  # Create a row with all the results for the current species
  species_name <- gsub("^clim_data_|\\.csv$", "", basename(species_file)) %>%
    gsub(" ", "_", .)  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    n = number.individuals,
    pc1.var = pc1.var,
    pc2.var = pc2.var,
    pc1.high.temp.loading = pc1.high.temp.loading,
    pc1.low.temp.loading = pc1.low.temp.loading,
    pc1.moisture.loading = pc1.moisture.loading,
    pc2.high.temp.loading = pc2.high.temp.loading,
    pc2.low.temp.loading = pc2.low.temp.loading,
    pc2.moisture.loading = pc2.moisture.loading,
    pc1.full.mod.intercept = pc1.full.mod.intercept,
    pc1.full.mod.lat.slope = pc1.full.mod.lat.slope,
    pc1.full.mod.long.slope = pc1.full.mod.long.slope,
    pc1.full.mod.lat.p.value = pc1.full.mod.lat.p.value,
    pc1.full.mod.long.p.value = pc1.full.mod.long.p.value,
    pc1.full.mod.adj.R2 = pc1.full.mod.adj.R2,
    pc2.full.mod.intercept = pc2.full.mod.intercept,
    pc2.full.mod.lat.slope = pc2.full.mod.lat.slope,
    pc2.full.mod.long.slope = pc2.full.mod.long.slope,
    pc2.full.mod.lat.p.value = pc2.full.mod.lat.p.value,
    pc2.full.mod.long.p.value = pc2.full.mod.long.p.value,
    pc2.full.mod.adj.R2 = pc2.full.mod.adj.R2,
    pc1.lat.mod.intercept = pc1.lat.mod.intercept,
    pc1.lat.mod.lat.slope = pc1.lat.mod.lat.slope,
    pc1.lat.mod.lat.p.value = pc1.lat.mod.lat.p.value,
    pc1.lat.mod.adj.R2 = pc1.lat.mod.adj.R2,
    pc2.lat.mod.intercept = pc2.lat.mod.intercept,
    pc2.lat.mod.lat.slope = pc2.lat.mod.lat.slope,
    pc2.lat.mod.lat.p.value = pc2.lat.mod.lat.p.value,
    pc2.lat.mod.adj.R2 = pc2.lat.mod.adj.R2,
    pc1.long.mod.intercept = pc1.long.mod.intercept,
    pc1.long.mod.long.slope = pc1.long.mod.long.slope,
    pc1.long.mod.long.p.value = pc1.long.mod.long.p.value,
    pc1.long.mod.adj.R2 = pc1.long.mod.adj.R2,
    pc2.long.mod.intercept = pc2.long.mod.intercept,
    pc2.long.mod.long.slope = pc2.long.mod.long.slope,
    pc2.long.mod.long.p.value = pc2.long.mod.long.p.value,
    pc2.long.mod.adj.R2 = pc2.long.mod.adj.R2,
    stringsAsFactors = FALSE
  )
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}

write.csv(results, file = "./Results/Q1.PCA.results.csv")

results.2 = results %>%
  filter(n > 99)

# Range of R2 for full, lat, and long models

range(results.2$pc1.full.mod.adj.R2)
# 0.008429213 0.932228418
range(results.2$pc2.full.mod.adj.R2)
# -0.004588984  0.948623628
range(results.2$pc1.lat.mod.adj.R2)
# -0.002049352  0.927826058
range(results.2$pc1.long.mod.adj.R2)
# -0.00718477  0.74997710

# just focus on PC1 for now
# number of species with greater than 50% variance explained
results.2 %>%
  filter(pc1.full.mod.adj.R2 > 0.50) %>%
  nrow()
# 76 out of 124 (61%)

# number of species with less than 30% variance explained
results.2 %>%
  filter(pc1.full.mod.adj.R2 < 0.30) %>%
  nrow()
# 31 out of 124 (25%)

results.2 %>%
  filter(pc1.lat.mod.adj.R2 > 0.50) %>%
  nrow()
# 58 out of 124 (47%)

results.2 %>%
  filter(pc1.long.mod.adj.R2 > 0.50) %>%
  nrow()
# 20 out of 124 (16%)


# number of species with significant lat in full model with pc1
results.2 %>%
  filter(pc1.full.mod.lat.p.value < 0.05) %>%
  nrow()
# 121 out of 124 (98%) - even though very little of the variation is explained

# number of species with significant long in full model with pc1
results.2 %>%
  filter(pc1.full.mod.long.p.value < 0.05) %>%
  nrow()
# 114 out of 124 (92%)

# number of species with significant lat in lat model with pc1
results.2 %>%
  filter(pc1.lat.mod.lat.p.value < 0.05) %>%
  nrow()
# 119 out of 124 (96%) - even though very little of the variation is explained

# number of species with significant long in long model with pc1
results.2 %>%
  filter(pc1.long.mod.long.p.value < 0.05) %>%
  nrow()
# 112 out of 124 (90%)

# range of intercept and slope values of full model
range(results.2$pc1.full.mod.intercept)
# -258.4939  115.3228
range(results.2$pc1.full.mod.lat.slope)
#-1.416994  2.787474
range(results.2$pc1.full.mod.long.slope)
#-1.9336017  0.8641082

ggplot(results.2, aes(pc1.full.mod.intercept))+
  geom_density()

ggplot(results.2, aes(pc1.full.mod.lat.slope))+
  geom_density()
ggplot(results.2, aes(pc1.full.mod.long.slope))+
  geom_density()
  
ggplot(results.2, aes(pc1.lat.mod.lat.slope))+
  geom_density()
ggplot(results.2, aes(pc1.long.mod.long.slope))+
  geom_density()

#### Multivariate Multiple Regression ####
# https://library.virginia.edu/data/articles/getting-started-with-multivariate-multiple-regression
# better than multiple single regression b/c it modifies hypothesis tests for regression parameters 
# and our confidence intervals for predictions.

species.list = list.files("./Results/Species.Climate/", pattern = "*.csv", full.names = TRUE)

species.1 = read.csv(species.list[[1]])
species.1.cols = species.1[,c(2:6)]
species.1.v2 = species.1.cols %>%
  mutate(scale.high.temp = scale(high_temp_C),
         scale.low.temp = scale(low_temp_C),
         scale.moisture = scale(moisture_mm))

mod.1 = lm(cbind(scale.high.temp,scale.low.temp,scale.moisture) ~ decimalLongitude + decimalLatitude, 
           data = species.1.v2)
summary(mod.1)

Y = as.matrix(species.1.v2[,c("scale.high.temp","scale.low.temp","scale.moisture")])
mod.2 = lm(Y ~ decimalLongitude + decimalLatitude, 
           data = species.1.v2)
summary(mod.2)

plot(resid(mod.1)[,1]~fitted(mod.1)[,1])
lines(lowess(resid(mod.1)[,1]~fitted(mod.1)[,1]), col="red")
plot(resid(mod.1)[,2]~fitted(mod.1)[,2])
lines(lowess(resid(mod.1)[,2]~fitted(mod.1)[,2]), col="red")
plot(resid(mod.1)[,3]~fitted(mod.1)[,3])
lines(lowess(resid(mod.1)[,3]~fitted(mod.1)[,3]), col="red")

rst = rstandard(mod.1)

qqnorm(rst[,1])
qqline(rst[,1])
qqnorm(rst[,2])
qqline(rst[,2])
qqnorm(rst[,3])
qqline(rst[,3])

plot(sqrt(abs(rst[,1]))~fitted(mod.1)[,1])
lines(lowess(sqrt(abs(rst[,1]))~fitted(mod.1)[,1]), col="red")
plot(sqrt(abs(rst[,2]))~fitted(mod.1)[,2])
lines(lowess(sqrt(abs(rst[,2]))~fitted(mod.1)[,2]), col="red")
plot(sqrt(abs(rst[,3]))~fitted(mod.1)[,3])
lines(lowess(sqrt(abs(rst[,3]))~fitted(mod.1)[,3]), col="red")
             
library(car)
Anova(mod.1)
confint(mod.1)
manova(mod.1)

anova(mod.1, test = "Wilks")
Anova(mod.1, test = "Wilks")

adonis(full.rda)
             
