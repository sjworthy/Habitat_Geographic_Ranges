# Testing for phylogenetic signal
# https://www.sciencedirect.com/science/article/pii/S2468265922001329
# https://github.com/jinyizju/U.PhyloMaker?tab=readme-ov-file
# https://github.com/megatrees

library(devtools)
#devtools::install_github("jinyizju/U.PhyloMaker")
library(U.PhyloMaker)
library(tidyverse)
#devtools::install_github("ecoinfor/U.Taxonstand")
library(U.Taxonstand)
library(phytools)
library(ape)
library(caper)
library(phangorn)
library(picante)

#### Building the Phylogeny ####
# species list including species, genus, family, species.relative, and genus.relative
# first column is mandatory and others are optional
sp.list <- read.csv("./Raw.Data/species.csv")

# Verify taxonomy
taxon.list = sp.list %>%
  filter(US.extent != "BAD") %>%
  dplyr::select(Species) %>%
  mutate(Rank = 2)

# make changes to taxonomy based on known 
# Fraxinus texensis = Fraxinus albicans
# Quercus prinus = Quercus michauxii so prinus is removed
# Taxodium ascendens = Taxodium distichum, ascendens is removed

taxon.list[119,1] = "Fraxinus albicans"
taxon.list.2 = taxon.list %>%
  filter(!Species %in% c("Quercus prinus","Taxodium ascendens","Aesculus flava"))
# remove A. flava b/c of low sample size
# 122 species
colnames(taxon.list.2)[1] = "Name"

# https://github.com/nameMatch/Database/tree/main/Plants_WFO
load("./Raw.Data/Plants_WFO.rdata") # World Flora Online

res <- nameMatch(spList=taxon.list.2, spSource=database, author = TRUE, max.distance= 1)

# Change Quercus margaretta to Quercus margarettae
# Change Quercus sinuata var. sinuata to Quercus sinuata

res.2 = res %>%
  select(Submitted_Name,Accepted_SPNAME)

# corrected species list
sp.list.2 = sp.list %>%
  filter(US.extent != "BAD") %>%
  dplyr::select("Species","Genus")

sp.list.2 = sp.list.2 %>%
  filter(!Species %in% c("Quercus prinus","Taxodium ascendens","Aesculus flava"))

sp.list.2[122,1] = "Quercus margarettae"
sp.list.2[100,1] = "Quercus sinuata"
sp.list.2[119,1] = "Fraxinus albicans"
# 122 species

colnames(sp.list.2) = c("species","genus")

# megatree for plants
megatree <- read.tree("https://raw.githubusercontent.com/megatrees/plant_20221117/refs/heads/main/plant_megatree.tre")

# genus family relationship file
gen.list <- read.csv("./Formatted.Data/plant_genus_list.csv")

# generate a phylogeny for the sample species list
result <- phylo.maker(sp.list.2, megatree, gen.list, nodes.type = 1, scenario = 3)

#write.tree(result$phylo, file = "./Results/phylo.tre")
#write.csv(result$sp.list, file = "./Results/phylo_splist.csv")

plot.phylo(result$phylo, cex = 0.5)

#### Testing for Phylogenetic Signal ####

# read in the data
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# add new species column with _ between genus and species to match phylogeny
microclim = microclim %>%
  mutate(species.2 = str_replace(species, " ", "_"))
topo = topo %>%
  mutate(species.2 = str_replace(species, " ", "_"))
soil = soil %>%
  mutate(species.2 = str_replace(species, " ", "_"))

# Change Quercus margaretta to Quercus margarettae
microclim <- microclim %>%
  mutate(species.2 = recode(species.2, "Quercus_margaretta" = "Quercus_margarettae"))
topo <- topo %>%
  mutate(species.2 = recode(species.2, "Quercus_margaretta" = "Quercus_margarettae"))
soil <- soil %>%
  mutate(species.2 = recode(species.2, "Quercus_margaretta" = "Quercus_margarettae"))

# read in the phylo
phylo = read.tree("./Results/phylo.tre")

# phylogenetic signal as category
# determined by quantifying the parsimony Sankoff score calculated from the distribution of trait
# categories on the phylogeny using the phangorn package. The significance of the score determined
# by randomly shuffling the species on the tips of the phylo 999 times to generate a null
# distribution that is compared to the observed parsimony score to calculate a P. value.
# P < 0.05 indicative of closely related species having similar traits.

# microclim
cat.microclim = as.data.frame(microclim$Category)
rownames(cat.microclim) = microclim$species.2
colnames(cat.microclim) = "microclim.category"

microclim.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  microclim.obs = t(data.frame(cat.microclim))
  microclim.obs2 = phyDat(t(microclim.obs), type = "USER", levels = attributes(factor(microclim.obs))$levels)
  microclim.null[i] = parsimony(rand.tree, microclim.obs2, method = "sankoff")
}

microclim.obs.real = parsimony(phylo, microclim.obs2, method = "sankoff")
# 48
microclim.p.value = (rank(c(microclim.obs.real,microclim.null))[1])/1000
# 0.0645

# topo
cat.topo = as.data.frame(topo$Category)
rownames(cat.topo) = topo$species.2
colnames(cat.topo) = "topo.category"

topo.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  topo.obs = t(data.frame(cat.topo))
  topo.obs2 = phyDat(t(topo.obs), type = "USER", levels = attributes(factor(topo.obs))$levels)
  topo.null[i] = parsimony(rand.tree, topo.obs2, method = "sankoff")
}

topo.obs.real = parsimony(phylo, topo.obs2, method = "sankoff")
# 51
topo.p.value = (rank(c(topo.obs.real,topo.null))[1])/1000
# 0.1015

# soil
cat.soil = as.data.frame(soil$Category)
rownames(cat.soil) = soil$species.2
colnames(cat.soil) = "soil.category"

soil.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  soil.obs = t(data.frame(cat.soil))
  soil.obs2 = phyDat(t(soil.obs), type = "USER", levels = attributes(factor(soil.obs))$levels)
  soil.null[i] = parsimony(rand.tree, soil.obs2, method = "sankoff")
}

soil.obs.real = parsimony(phylo, soil.obs2, method = "sankoff")
# 53
soil.p.value = (rank(c(soil.obs.real,soil.null))[1])/1000
# 0.4635

#### Testing for Phylogenetic Signal combo patterns ####
# read in the data
dat = read.csv("./Results/occupancy.patterns.csv")

# add new species column with _ between genus and species to match phylogeny
dat = dat %>%
  mutate(species.2 = str_replace(X, " ", "_"))

# Change Quercus margaretta to Quercus margarettae
dat <- dat %>%
  mutate(species.2 = recode(species.2, "Quercus_margaretta" = "Quercus_margarettae"))

# read in the phylo
phylo = read.tree("./Results/phylo.tre")

# phylogenetic signal as category
# determined by quantifying the parsimony Sankoff score calculated from the distribution of trait
# categories on the phylogeny using the phangorn package. The significance of the score determined
# by randomly shuffling the species on the tips of the phylo 999 times to generate a null
# distribution that is compared to the observed parsimony score to calculate a P. value.
# P < 0.05 indicative of closely related species having similar traits.

# make data frames for the groups
combo.letter = as.data.frame(dat$combo.letter)
rownames(combo.letter) = dat$species.2
colnames(combo.letter) = "combo.letter"

clim.topo.letter = as.data.frame(dat$clim.topo.letter)
rownames(clim.topo.letter) = dat$species.2
colnames(clim.topo.letter) = "clim.topo.letter"

clim.soil.letter = as.data.frame(dat$micro.soil.letter)
rownames(clim.soil.letter) = dat$species.2
colnames(clim.soil.letter) = "clim.soil.letter"

topo.soil.letter = as.data.frame(dat$topo.soil.letter)
rownames(topo.soil.letter) = dat$species.2
colnames(topo.soil.letter) = "topo.soil.letter"

combo.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  combo.obs = t(data.frame(combo.letter))
  combo.obs2 = phyDat(t(combo.obs), type = "USER", levels = attributes(factor(combo.obs))$levels)
  combo.null[i] = parsimony(rand.tree, combo.obs2, method = "sankoff")
}

combo.obs.real = parsimony(phylo, combo.obs2, method = "sankoff")
# 95
combo.p.value = (rank(c(combo.obs.real,combo.null))[1])/1000
# 0.1765

clim.topo.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  clim.topo.obs = t(data.frame(clim.topo.letter))
  clim.topo.obs2 = phyDat(t(clim.topo.obs), type = "USER", levels = attributes(factor(clim.topo.obs))$levels)
  clim.topo.null[i] = parsimony(rand.tree, clim.topo.obs2, method = "sankoff")
}

clim.topo.obs.real = parsimony(phylo, clim.topo.obs2, method = "sankoff")
# 81
clim.topo.p.value = (rank(c(clim.topo.obs.real,clim.topo.null))[1])/1000
# 0.1625

clim.soil.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  clim.soil.obs = t(data.frame(clim.soil.letter))
  clim.soil.obs2 = phyDat(t(clim.soil.obs), type = "USER", levels = attributes(factor(clim.soil.obs))$levels)
  clim.soil.null[i] = parsimony(rand.tree, clim.soil.obs2, method = "sankoff")
}

clim.soil.obs.real = parsimony(phylo, clim.soil.obs2, method = "sankoff")
# 77
clim.soil.p.value = (rank(c(clim.soil.obs.real,clim.soil.null))[1])/1000
# 0.189

topo.soil.null = c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(phylo)
  topo.soil.obs = t(data.frame(topo.soil.letter))
  topo.soil.obs2 = phyDat(t(topo.soil.obs), type = "USER", levels = attributes(factor(topo.soil.obs))$levels)
  topo.soil.null[i] = parsimony(rand.tree, topo.soil.obs2, method = "sankoff")
}

topo.soil.obs.real = parsimony(phylo, topo.soil.obs2, method = "sankoff")
# 78
topo.soil.p.value = (rank(c(topo.soil.obs.real,topo.soil.null))[1])/1000
# 0.122
