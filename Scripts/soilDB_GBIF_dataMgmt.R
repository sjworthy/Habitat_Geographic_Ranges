
### Read in packages
library(soilDB)
library(dplyr)
library(sf)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

### Load Lat/Lon CSV 
gbif = read.csv("good.gbif.data.csv", row.names = 1)

# Set your folder path
folder_path <- "/Users/bhm11001/Library/CloudStorage/OneDrive-MichiganStateUniversity/Post-Doc/Microclimate Habitat Project/processed_species"

# List all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read all CSVs into a list
csv_list <- lapply(csv_files, read.csv)

# Name the files
names(csv_list) <- tools::file_path_sans_ext(basename(csv_files))

# Combine csvs into a dataframe
allSpp_soil <- do.call(rbind, csv_list)

# Remove row names to avoid long names like 'Species_soil.1'
rownames(allSpp_soil) <- NULL

# Just pull out columns associated with the 0-100 cm depth
allSpp_soil_final <- allSpp_soil[,c(1:3,25:31)]

# Write csv
write.csv(allSpp_soil_final, "gbif.data.soils.0.100cm.csv", row.names = FALSE)

# Pull out columns with soil variables
cols <- names(allSpp_soil)
soil_cols <- cols[match("ph_d0_20", cols):match("texture_d0_100", cols)]

# Define the soil column range
# soil_cols <- names(allSpp_soil)[match("ph_d0_20", names(allSpp_soil)):match("texture_d0_100", names(allSpp_soil))]

# Calculate percent of rows with all soil NAs per species
na_summary_by_species <- allSpp_soil %>%
  group_by(species) %>%
  summarise(
    n_obs = n(),
    n_all_soilNA = sum(if_all(all_of(soil_cols), is.na)),
    percent_all_soilNA = 100 * n_all_soilNA / n_obs
  ) %>%
  arrange(desc(percent_all_soilNA))

#View(na_summary_by_species)
write.csv(na_summary_by_species, "sp_soilDat_summary.csv")

# All soil columns are NA
n_all_na <- allSpp_soil %>%
  filter(if_all(all_of(soil_cols), is.na)) %>%
  nrow()

n_all_na # 63781 out of 1336865 (~5%)

# 2. Any soil column is NA
some_na_rows <- allSpp_soil %>%
  filter(if_any(all_of(soil_cols), is.na) & !if_all(all_of(soil_cols), is.na))
n_some_na <- nrow(some_na_rows)
n_some_na # 19647 out of 1336865 (1.5%)

# First, filter the rows where some but not all soil values are NA
some_na_only <- allSpp_soil %>%
  filter(if_any(all_of(soil_cols), is.na) & !if_all(all_of(soil_cols), is.na))
# Now summarize NA counts for only those rows
na_summary_some_only <- some_na_only %>%
  summarise(across(all_of(soil_cols), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "num_NA") %>%
  filter(num_NA > 0) %>%
  arrange(desc(num_NA))
print(na_summary_some_only, n = 28)

## EC - lots of NAs and missing data at all depths
## pH - least NAs at 0-100 cm
## Texture, clay, sand, silt - 0-100 cm
# Bulk density (db) - 0-50 or 0-100 cm, no difference

# Now summarize NA counts for only those rows
na_summary <- allSpp_soil %>%
  summarise(across(all_of(soil_cols), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "num_NA") %>%
  filter(num_NA > 0) %>%
  arrange(desc(num_NA))
print(na_summary, n = 28)
