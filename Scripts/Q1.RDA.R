# RDA models to address question 1

# microclimate response variables: high temp, low temp, moisture
# microclimate predictors: lat and long

# split into species dataframes
# scale and center responses variables since temp and moisture on very differenct scales

library(vegan)
library(tidyverse)

#### RDA ####

species.list = list.files(".", pattern = "*.csv", full.names = TRUE)

results <- data.frame(
  species = character(),
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
  lat_anova_p_value = numeric(),
  long_anova_p_value = numeric(),
  stringsAsFactors = FALSE
)

for(species_file in species.list){
  
  # Read in the species data
  species = read.csv(species_file)
  
  # create response matrix, scaling variables b/c different units
  response = species %>%
    select(high_temp_C,low_temp_C,moisture_mm) %>%
    scale()  
  
  # Full rda with both latitude and longitude, get adjusted R2, constrained and unconstrained variance
  full.rda = rda(response ~ decimalLatitude + decimalLongitude, species)
  full_r2 <- RsquareAdj(full.rda)$adj.r.squared
  full.rda.summary = summary(full.rda)
  full.constrained.var = full.rda.summary$cont$importance[2, ]
  full.unconstrained.var = full.rda.summary$cont$importance[3, ]
  
  # significance of rda
  full.anova = anova.cca(full.rda)
  full_anova_p_value <- full_anova$`Pr(>F)`[1]
###########  full.anova.terms = anova.cca(full.rda, by = "terms")
  
  # get correlation between lat and long for the species
  geo.cor = cor(species$decimalLatitude, species$decimalLongitude)
  
  # rda with latitude, get adjusted R2, constrained and unconstrained variance
  lat.rda = rda(response ~ decimalLatitude, species)
  lat_r2 <- RsquareAdj(lat.rda)$adj.r.squared
  lat.rda.summary = summary(lat.rda)
  lat.constrained.var = lat.rda.summary$cont$importance[2, ]
  lat.unconstrained.var = lat.rda.summary$cont$importance[3, ]
  
  # significance of rda
  lat.anova = anova.cca(lat.rda)
  lat_anova_p_value <- lat_anova$`Pr(>F)`[1]
  
  # rda with longitude, get adjusted R2, constrained and unconstrained variance
  long.rda = rda(response ~ decimalLongitude, species)
  long_r2 <- RsquareAdj(long.rda)$adj.r.squared
  long.rda.summary = summary(long.rda)
  long.constrained.var = long.rda.summary$cont$importance[2, ]
  long.unconstrained.var = long.rda.summary$cont$importance[3, ]
  
  # significance of rda
  long.anova = anova.cca(long.rda)
  long_anova_p_value <- long_anova$`Pr(>F)`[1]
  
  # Create a row with all the results for the current species
  species_name <- gsub(".csv", "", basename(species_file))  # Extract species name from file name
  
  species_results <- data.frame(
    species = species_name,
    full_r2 = full_r2,
    full_constrained_var = full_constrained_var,
    full_unconstrained_var = full_unconstrained_var,
    geo_cor = geo_cor,
    lat_r2 = lat_r2,
    lat_constrained_var = lat_constrained_var,
    lat_unconstrained_var = lat_unconstrained_var,
    long_r2 = long_r2,
    long_constrained_var = long_constrained_var,
    long_unconstrained_var = long_unconstrained_var,
    full_anova_p_value = full_anova_p_value,
    lat_anova_p_value = lat_anova_p_value,
    long_anova_p_value = long_anova_p_value,
    stringsAsFactors = FALSE
  )
  
  # Add the results for this species to the results data frame
  results <- rbind(results, species_results)
}





A.fraseri.rda = rda(response.2 ~ decimalLatitude + decimalLongitude, A.fraseri)
full_r2 <- RsquareAdj(A.fraseri.rda)$adj.r.squared # 0.13
summary(A.fraseri.rda)
# constrained proportion: variance of Y explained by X
# unconstrained proportion: unexplained variance in Y

vif.cca(A.fraseri.rda) # need to be below 3

# significance testing using anova.cca, automatically does 999 permutations.
# F statistic correspons to an overall test of sign. of an RDA by comparing the computed model to
# a null model. This test is based on the null hypothesis that the strength of the linear relationship
# calculated by the R2 is not larger than the value that would be obtained for unrelated Y and X matrices
# of the same size. 
anova.cca(A.fraseri.rda)
anova.cca(A.fraseri.rda, by = "terms")

ordiplot(A.fraseri.rda, scaling = 1, type = "text")
ordiplot(A.fraseri.rda, scaling = 2, type = "text")

A.fraseri.rda.lat = rda(response.2 ~ decimalLatitude, A.fraseri)
RsquareAdj(A.fraseri.rda.lat)$adj.r.squared
A.fraseri.rda.long = rda(response.2 ~ decimalLongitude, A.fraseri)
RsquareAdj(A.fraseri.rda.long)$adj.r.squared
anova.cca(A.fraseri.rda.long)

test = as.data.frame(unique(clim.data$species))

