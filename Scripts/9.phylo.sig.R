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


# species list including species, genus, family, species.relative, and genus.relative
# first column is mandatory and others are optional
sp.list <- read.csv("./Raw.Data/species.csv")

# Verify taxonomy
taxon.list = sp.list %>%
  filter(US.extent != "BAD") %>%
  select(Species) %>%
  mutate(Rank = 2)

# make changes to taxonomy based on known 
# Fraxinus texensis = Fraxinus albicans
# Quercus prinus = Quercus michauxii so prinus is removed
# Taxodium ascendens = Taxodium distichum, ascendens is removed

taxon.list[40,1] = "Fraxinus albicans"
taxon.list.2 = taxon.list[-c(102,114),]
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
  select("Species","Genus")

sp.list.2[93,1] = "Quercus margarettae"
sp.list.2[105,1] = "Quercus sinuata"
sp.list.2[40,1] = "Fraxinus albicans"
sp.list.3 = sp.list.2[-c(102,114),]
# 122 species

colnames(sp.list.3) = c("species","genus")

# megatree for plants
megatree <- read.tree("https://raw.githubusercontent.com/megatrees/plant_20221117/refs/heads/main/plant_megatree.tre")

# genus family relationship file
gen.list <- read.csv("./Formatted.Data/plant_genus_list.csv")

# generate a phylogeny for the sample species list
result <- phylo.maker(sp.list.3, megatree, gen.list, nodes.type = 1, scenario = 3)

write.tree(result$phylo, file = "./Results/phylo.tre")
write.csv(result$sp.list, file = "./Results/phylo_splist.csv")

plot.phylo(result$phylo, cex = 0.5)

#### Testing for Phylogenetic Signal ####

microclim.dat = read.csv("./Results/Q1.MRM.results.csv")
topo.dat = read.csv("./Results/Q3.MRM.results.csv")

# Change Quercus margaretta to Quercus margarettae

microclim.dat[79,2] = "Quercus_margarettae"
topo.dat[79,2] = "Quercus_margarettae"

# remove the BAD range species

microclim.dat.2 = microclim.dat %>%
  filter(!species %in% c("Populus_deltoides","Thuja_occidentalis"))

topo.dat.2 = topo.dat %>%
  filter(!species %in% c("Populus_deltoides","Thuja_occidentalis"))

phylo = read.tree("./Results/phylo.tre")

# drop dips for species we don't have data for yet

pruned.tree = drop.tip(phylo, phylo$tip.label[-match(microclim.dat.2$species, phylo$tip.label)])

slope.names = setNames(microclim.dat.2$Slope,microclim.dat.2$species)

microclim_trait=phylosig(pruned.tree, slope.names, method = "K", test = TRUE, nsim = 10000)
microclim_trait
# 0.0325744, p = 0.0789

# test if significant K (0.0325744 is significantly different from 1)
nullX<-fastBM(pruned.tree, nsim=10000) # simulate 10000 datasets
# for each, carry out a test for phylogenetic signal
# and accumulate these into a vector using sapply
nullK<-apply(nullX, 2, phylosig, tree=pruned.tree)

# if K_trait less than 1
# calculate p-values by counting the proportion of times our simulated values of K were greater than our observed values.
pval.less = mean(nullK<=microclim_trait$K) 
pval.less # significantly less than 1

plotTree.barplot(pruned.tree, slope.names,args.plotTree=list(fsize=0.5),
                 args.barplot=list(xlab="microclim slope"))

topo.slope.names = setNames(topo.dat.2$Slope,topo.dat.2$species)

topo_trait=phylosig(pruned.tree, topo.slope.names, method = "K", test = TRUE, nsim = 10000)
topo_trait
# 0.00382547, p = 0.7594

plotTree.barplot(pruned.tree, topo.slope.names,args.plotTree=list(fsize=0.5),
                 args.barplot=list(xlab="topo slope"))

### Testing for phylogenetic signal in range size ####

range.size = read.csv("./Results/Range.Size.csv")

range.size.2 = na.omit(range.size)
range.size.2[102,1] = "Quercus_margarettae"

range.size.2$Species <- gsub(" ", "_", range.size.2$Species)

microclim.dat = read.csv("./Results/Q1.MRM.results.csv")

# Change Quercus margaretta to Quercus margarettae
microclim.dat[79,2] = "Quercus_margarettae"

# remove the BAD range species

microclim.dat.2 = microclim.dat %>%
  filter(!species %in% c("Populus_deltoides","Thuja_occidentalis"))

phylo = read.tree("./Results/phylo.tre")

# drop dips for species we don't have data for yet

pruned.tree = drop.tip(phylo, phylo$tip.label[-match(microclim.dat.2$species, phylo$tip.label)])

# prune range data by tip labels

range.size.3 <- range.size.2[range.size.2$Species%in% pruned.tree$tip.label, ]

cor.test(range.size.3$latRange,range.size.3$greatCircDist)
cor.test(range.size.3$meanPairDist,range.size.3$greatCircDist)
cor.test(range.size.3$minSpanTree,range.size.3$greatCircDist)
# All really correlated

GRD.range.names = setNames(range.size.3$greatCircDist,range.size.3$Species)

range_trait=phylosig(pruned.tree, GRD.range.names, method = "K", test = TRUE, nsim = 10000)
range_trait
# 0.00485223, p = 0.668

#### relationship between range size and slopes ####

pruned.tree.2 = drop.tip(pruned.tree, pruned.tree$tip.label[-match(range.size.3$Species, pruned.tree$tip.label)])

microclim.dat.3 = microclim.dat.2[microclim.dat.2$species%in% range.size.3$Species, ]
colnames(microclim.dat.3)[2] = "sp"
colnames(range.size.3)[1] = "sp"

# merge data

all.dat = merge(microclim.dat.3,range.size.3, by = "sp")

# getting error so removing node labels
pruned.tree.2$node.label<-NULL

comp.data<-comparative.data(pruned.tree.2, all.dat, names.col="sp", vcv.dim=2, warn.dropped=TRUE, vcv = TRUE)
micro.range.pgls = pgls(Slope~greatCircDist, data=comp.data)
summary(micro.range.pgls) # significant, negative p = 4.155e-08

ggplot(all.dat, aes(x = greatCircDist, y = Slope)) +
  geom_point()
