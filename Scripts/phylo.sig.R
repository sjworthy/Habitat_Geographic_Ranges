# Testing for phylogenetic signal
# https://www.sciencedirect.com/science/article/pii/S2468265922001329
# https://github.com/jinyizju/U.PhyloMaker?tab=readme-ov-file
# https://github.com/megatrees

library(devtools)
devtools::install_github("jinyizju/U.PhyloMaker")
library(U.PhyloMaker)
library(tidyverse)
devtools::install_github("ecoinfor/U.Taxonstand")
library(U.Taxonstand)

# species list including species, genus, family, species.relative, and genus.relative
# first column is mandatory and others are optional
sp.list <- read.csv("./Raw.Data/species.csv")

# Verify taxonomy
taxon.list = sp.list %>%
  select(Species,Scientific.Name)

taxon.list$Name <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", taxon.list$Scientific.Name)
taxon.list$Author <- sub("^[A-Za-z]+ [A-Za-z]+ (.*)", "\\1", taxon.list$Scientific.Name)

taxon.list.2 = taxon.list %>%
  select(Name) %>%
  mutate(Rank = 2)

# https://github.com/nameMatch/Database/tree/main/Plants_WFO
load("Plants_WFO.rdata") # World Flora Online

res <- nameMatch(spList=taxon.list.2, spSource=database, author = TRUE, max.distance= 1)

# Change Quercus margaretta to Quercus margarettae
# Change Quercus sinuata var. sinuata to Quercus sinuata




sp.list.2 = sp.list %>%
  filter(US.extent != "BAD") %>%
  select("Species","Genus")

colnames(sp.list.2) = c("species","genus")

# megatree for plants
megatree <- read.tree("https://raw.githubusercontent.com/megatrees/plant_20221117/refs/heads/main/plant_megatree.tre")

# genus family relationship file
gen.list <- read.csv("./Formatted.Data/plant_genus_list.csv")

# generate a phylogeny for the sample species list
result <- phylo.maker(sp.list.2, megatree, gen.list, nodes.type = 1, scenario = 3)

