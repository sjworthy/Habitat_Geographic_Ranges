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


### Building the Phylogeny ####
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
# 66
microclim.p.value = (rank(c(microclim.obs.real,microclim.null))[1])/1000
# 0.442

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
# 39
topo.p.value = (rank(c(topo.obs.real,topo.null))[1])/1000
# 0.6715

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
# 56
soil.p.value = (rank(c(soil.obs.real,soil.null))[1])/1000
# 0.605


### Plotting Categories on the phylogeny ####

library(ggtree)

phylo <- ggtree::read.tree("./Results/phylo.tre")
phylo$tip.label = gsub("_", " ", phylo$tip.label)

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo) +
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree()
tree

# combine the categories from all three variables
all.cats = as.data.frame(microclim$Category)
colnames(all.cats)[1] = "microclim.cat"
all.cats$topo.cat = topo$Category
all.cats$soil.cat = soil$Category
all.cats$species = microclim$species

all.cats <- all.cats %>%
  mutate(species = recode(species, "Quercus margaretta" = "Quercus margarettae"))

row.names(all.cats) = all.cats$species
all.cats = all.cats[,c(1:3)]

cat.phylo = gheatmap(tree, all.cats, offset=35, width=0.3, 
         colnames=FALSE, legend_title="Category") +
  scale_fill_manual(breaks=c("generalists", "overdisper.shifter", "overdisperser",
                             "shifting", "specialists"), 
                    values=c("#F5AF4D","#B46DB3","#548F01","#5495CF","#DB4743"),
                    name="Pattern",
                    labels = c(
                      "shifting" = "Shifter",
                      "specialists" = "Specialist",
                      "generalists" = "Generalist",
                      "overdisperser" = "Overdisperser",
                      "overdisper.shifter" = "Overdispersed Shifter"))
cat.phylo

ggsave("./Plots/ESA.plots/phylo.png", height = 10, width = 12)

### Heatmap for ESA ####

# read in the data
all.cats = read.csv("./Formatted.Data/cat.summary.csv")

# Add row ID to preserve order (useful since we’re dropping species names)
all.cats$RowID <- factor(seq_len(nrow(all.cats)))  # as factor to control plotting order

# Pivot longer for plotting
cats.long <- all.cats %>%
  pivot_longer(cols = c(Microclimate, Topography, Soil),
               names_to = "CategoryType",
               values_to = "Category") %>%
  mutate(Species = NA)  # placeholder to match value.column structure

# Create the fourth "Species" column as fake heatmap column
value.column <- data.frame(
  RowID = all.cats$RowID,
  CategoryType = "Species",  # this becomes a 4th column
  Category = NA,
  Species = all.cats$Species
)

# Combine data
cats.long <- bind_rows(cats.long, value.column)

# order of columns for plotting
cats.long$CategoryType <- factor(cats.long$CategoryType, levels = c("Microclimate", "Topography", 
                                                                    "Soil","Species"))

# Plot
heat.plot =ggplot() +
  # Only plot tiles for the heatmap columns
  geom_tile(
    data = subset(cats.long, CategoryType != "Species"),
    aes(x = CategoryType, y = RowID, fill = Category),
    color = "white"
  ) +
  # Add text for Species values in the last column
  geom_text(
    data = subset(cats.long, CategoryType == "Species"),
    aes(x = CategoryType, y = RowID, label = Species),
    size = 7
  ) +
  scale_fill_manual(
    breaks = c("generalists", "overdisper.shifter", "overdisperser", "shifting", "specialists"),
    values = c("#F5AF4D", "#B46DB3", "#548F01", "#5495CF", "#DB4743"),
    name = "Pattern",
    na.value = "transparent",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter"
    )) + 
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

heat.plot

ggsave("./Plots/ESA.plots/heatmap.png", height = 10, width = 12)
