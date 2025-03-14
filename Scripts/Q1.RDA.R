# RDA models to address question 1

# microclimate response variables: high temp, low temp, moisture
# microclimate predictors: lat and long

# split into species dataframes
# scale and center responses variables since temp and moisture on very differenct scales

library(vegan)
library(tidyverse)

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



