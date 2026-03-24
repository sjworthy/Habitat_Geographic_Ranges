# code to evaluate relationships among occupancy patterns and range size

library(tidyverse)
library(rstatix)
library(vegan)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggpubr)
library(ape)
library(phytools)

# ISSUES: EOO and CentoridX = equal variance, CentroidY unequal variance
# shifters and specialists not normally distributed (can't use Welch's test)

#### Combining all the ranges ####
# combining all the range size estimates together and running analyses of range size

# EOO normal range size, Extent of Occurrence
# area contained within the shortest continuous imaginary boundary which can be drawn to 
# emcompass all the sites of present occurrence of a taxon
# square kilometers
EOO.norm = read.csv("./Formatted.Data/EOO.norm.csv", row.names = 1)

# EOO alpha range size, Extent of Occurrence with alpha hull
# Uses an alpha hull instead of a convex hull
EOO.alpha = read.csv("./Formatted.Data/EOO.alpha.csv", row.names = 1)

# all other range size metrics
range.size = read.csv("./Results/Range.Size.csv")
# NA for largest number of occurrence species for last two range metrics - too large to compute

# merge all together
merged_df <- EOO.norm %>%
  full_join(EOO.alpha, by = "tax") %>%
  full_join(range.size, by = c("tax" = "Species"))

colnames(merged_df)[1] = "Species"
colnames(merged_df)[2] = "EOO.norm"
colnames(merged_df)[3] = "EOO.norm.issues"
colnames(merged_df)[4] = "EOO.alpha"
colnames(merged_df)[5] = "EOO.alpha.issue"

#write.csv(merged_df, file = "./Formatted.Data/all.range.estimates.csv")

#### Visualizing correlation among range estimates ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)


cor(ranges[, c("EOO.norm", "EOO.alpha", "centroidX", "centroidY", "latRange",
               "greatCircDist")], use = "pairwise.complete.obs", method = "pearson")

# EOO.norm and EOO alpha essentially the same and highly correlated with latRange and greatCircDist
# centroidX and centroidY not highly correlated with anything

# test EOO.alpha, centroidX, and centroidY since least correlated

### Box plots of range and microclimate patterns ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data

microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df

all.dat = inner_join(microclim, ranges, by = c("species" = "Species"))

# Plotting
ggplot(all.dat, aes(x = Category, y = EOO.alpha)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "EOO.alpha")
ggplot(all.dat, aes(x = Category, y = centroidX)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidX")
ggplot(all.dat, aes(x = Category, y = centroidY)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidY")

### Box plots of range and topography patterns ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data

topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df

all.dat = inner_join(topo, ranges, by = c("species" = "Species"))

# Plotting
ggplot(all.dat, aes(x = Category, y = EOO.alpha)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "EOO.alpha")
ggplot(all.dat, aes(x = Category, y = centroidX)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidX")
ggplot(all.dat, aes(x = Category, y = centroidY)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidY")

### Box plots of range and soil patterns ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data

soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df

all.dat = inner_join(soil, ranges, by = c("species" = "Species"))

# Plotting
ggplot(all.dat, aes(x = Category, y = EOO.alpha)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "EOO.alpha")
ggplot(all.dat, aes(x = Category, y = centroidX)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidX")
ggplot(all.dat, aes(x = Category, y = centroidY)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Category", y = "centroidY")

#### Kruskal-Wallis for group differences microclim ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in microclim data
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(microclim, ranges, by = c("species" = "Species"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, n = 122, df = 4, p = 0.0000127
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# generalists significantly larger overdispersed shifter, p adj = 0.000908
# overdispersed shifter significantly smaller overdispersed, p adj = 0.0425
# overdispersed shifter significantly smaller shifting, p adj = 0.0193
# overdispersed shifter significantly smaller specialists, only for p value = 0.00841, p adj = 0.0589

EOO.alpha.posthoc.2 <- EOO.alpha.posthoc %>% add_xy_position(x = "Category")

# EOO plot
EOO.plot = ggplot(all.dat, aes(x = Category, y = EOO.alpha, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3")) +
  theme_classic(base_size = 15) +
  labs(y = expression("Extent of Occurrence (km"^2*")")) +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalist", 
                            "overdisper.shifter"= "Overdispersed\nShifter",
                            "overdisperser" = "Overdisperser", "shifting" = "Shifter",
                            "specialists" = "Specialist"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)+
  ggtitle("Microclimate")
EOO.plot

ggsave("./Plots/microclim.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # not significant, p = 0.899
centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # significant, p = 0.395

#### Krustal-Wallis for group differences topography ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in microclim data
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(topo, ranges, by = c("species" = "Species"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, df = 3, n = 122, p = 0.00496
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# overdisper.shifter lower than shifter, p adj = 0.0111

EOO.alpha.posthoc.2 <- EOO.alpha.posthoc %>% add_xy_position(x = "Category")

# EOO plot
EOO.plot = ggplot(all.dat, aes(x = Category, y = EOO.alpha, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3")) +
  theme_classic(base_size = 15) +
  labs(y = expression("Extent of Occurrence (km"^2*")")) +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalist", 
                            "overdisper.shifter"= "Overdispersed\nShifter",
                            "overdisperser" = "Overdisperser", "shifting" = "Shifter",
                            "specialists" = "Specialist"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)+
  ggtitle("Topography")
EOO.plot

ggsave("./Plots/topo.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # significant, p = 0.848
centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # not significant, p = 0.21

#### Krustal-Wallis for group differences soil ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(soil, ranges, by = c("species" = "Species"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, p = 0.0196
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# none significantly different after p-value adjustment
# generalists higher than overdispersed shifter p = 0.00633
# generalists higher than shifters p = 0.0407
# overdispersed shifter less than overdisperser p = 0.0277
# overdispersed shifter less than specialists p = 0.0107

EOO.alpha.posthoc.2 <- EOO.alpha.posthoc %>% add_xy_position(x = "Category")

# EOO plot
EOO.plot = ggplot(all.dat, aes(x = Category, y = EOO.alpha, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3")) +
  theme_classic(base_size = 15) +
  labs(y = expression("Extent of Occurrence (km"^2*")")) +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalist", 
                            "overdisper.shifter"= "Overdispersed\nShifter",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifter",
                            "specialists" = "Specialist"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)+
  ggtitle("Soil")
EOO.plot

ggsave("./Plots/soil.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # significant, p = 1.93e-10
centroidX.posthoc = dunn_test(centroidX ~ Category, data = all.dat)
centroidX.posthoc

centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # not significant, p = 0.0117
centroidY.posthoc = dunn_test(centroidY ~ Category, data = all.dat)
centroidY.posthoc

#### Testing for phylogenetic signal in range size ####

# ranges
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in the phylo
phylo = read.tree("./Results/phylo.tre")

ranges[93,1] = "Quercus margarettae"
ranges = ranges %>%
  mutate(species.2 = str_replace(Species, " ", "_"))

EOO.range = setNames(ranges$EOO.norm, ranges$species.2)

set.seed(13) # set the seed so K is always the same

K_EOO = phylosig(phylo, EOO.range, method = "K", test = TRUE, nsim = 10000)
K_EOO # K = 0.00410193, p = 0.7507


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




