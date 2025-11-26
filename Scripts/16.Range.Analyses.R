library(tidyverse)
library(rstatix)
library(vegan)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggpubr)

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

write.csv(merged_df, file = "./Formatted.Data/all.range.estimates.csv")

#### Visualizing correlation among range estimates ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)


cor(ranges[, c("EOO.norm", "EOO.alpha", "centroidX", "centroidY", "latRange",
               "greatCircDist")], use = "pairwise.complete.obs", method = "pearson")

# EOO.norm and EOO alpha essentially the same and highly correlated with latRange and greatCircDist
# centroidX and centroidY not highly correlated with anything

# test EOO.alpha, centroidX, and centroidY since least correlated

### box plots of range and categories ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data

soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# for soil, just two categories

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

#### Krustal-Wallis for group differences soil ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(soil, ranges, by = c("species" = "Species"))

# get rid of overdisperser b/c only 2

all.dat = all.dat %>%
  filter(!Category %in% c("overdisperser"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, p = 0.000714
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# overdisper.shift sig. diff from specialists
# shifters sig. diff from specialists

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
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdispersed\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)
EOO.plot

ggsave("./Plots/ESA.plots/soil.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # significant, p = 0.0000949
centroidX.posthoc = dunn_test(centroidX ~ Category, data = all.dat)
centroidX.posthoc
# generalists sig. diff from shifting
# overdisper.shift diff from shifting
# overdisper.shift sig. diff from specialists

centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # not significant, p = 0.0974

#### Krustal-Wallis for group differences microclim ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in microclim data
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(microclim, ranges, by = c("species" = "Species"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, p = 0.00000000914
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# overdisperse shifters different from all other categories

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
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdispersed\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)
EOO.plot

ggsave("./Plots/ESA.plots/microclim.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # not significant, p = 0.268

centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # significant, p = 0.00752
centroidY.posthoc = dunn_test(centroidY ~ Category, data = all.dat)
centroidY.posthoc
# overdisper shifters sig. diff. from generalists

#### Krustal-Wallis for group differences topography ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in microclim data
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(topo, ranges, by = c("species" = "Species"))

# Remove overdisperser since only 1 species

all.dat = all.dat %>%
  filter(!Category %in% c("overdisperser"))

all.dat$Category = as.factor(all.dat$Category)

# using Krustal-Wallis 
EOO.alpha.test = kruskal_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.test # significant, p = 0.000013
EOO.alpha.posthoc = dunn_test(EOO.alpha ~ Category, data = all.dat)
EOO.alpha.posthoc
# overdisper.shift sig. diff from generalists and shifters

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
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdispersed\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none") +
  stat_pvalue_manual(EOO.alpha.posthoc.2, hide.ns = TRUE)
EOO.plot

ggsave("./Plots/ESA.plots/topo.EOO.png", width = 7, height = 5)

centroidX.test = kruskal_test(centroidX ~ Category, data = all.dat)
centroidX.test # significant, p = 0.0497
centroidX.posthoc = dunn_test(centroidX ~ Category, data = all.dat)
centroidX.posthoc
# none are different so p ~ 0.05

centroidY.test = kruskal_test(centroidY ~ Category, data = all.dat)
centroidY.test # not significant, p = 0.475

#### PERMANOVA for group differences ####
# can't use bray b/c of negative values

all.dat$Category = as.factor(all.dat$Category)

range.vars <- cbind(all.dat$EOO.alpha,all.dat$centroidX,all.dat$centroidY)
adonis_result <- adonis2(range.vars ~ Category, data = all.dat, perm = 999, method = "euclidean")
adonis_result # significant p = 0.005


pairwise.adonis(range.vars, all.dat$Category, sim.method = "euclidian",
                sim.function = "vegdist", p.adjust.m = "holm")

# shifters and specialists significantly differ

#### PERMANOVA without overdispersers ####

range.vars <- cbind(all.dat.2$EOO.alpha,all.dat.2$centroidX,all.dat.2$centroidY)
adonis_result <- adonis2(range.vars ~ Category, data = all.dat.2, perm = 999, method = "euclidean")
adonis_result # significant p = 0.001

pairwise.adonis(range.vars, all.dat.2$Category, sim.method = "euclidian",
                sim.function = "vegdist", p.adjust.m = "holm")

# shifters and specialists significantly differ
# overdisper shifter and specialists significanlly differ

### Welch's ANOVA ####
# handles unequal variances better, but need normal distributions

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(soil, ranges, by = c("species" = "Species"))

all.dat.2 = all.dat %>%
  filter(!Category %in% "overdisperser")

# test for normality
all.dat.2 %>%
  group_by(Category) %>%
  summarise(EOO_p_value = shapiro.test(EOO.alpha)$p.value)
# shifters and specialists not normal
all.dat.2 %>%
  group_by(Category) %>%
  summarise(centroidX_p_value = shapiro.test(centroidX)$p.value)
# shifters and specialists not normal
all.dat.2 %>%
  group_by(Category) %>%
  summarise(centroidY_p_value = shapiro.test(centroidY)$p.value)
# shifters not normal

# Test for equal variances
bartlett.test(EOO.alpha ~ Category, data = all.dat.2) # NS
bartlett.test(centroidX ~ Category, data = all.dat.2) # NS
bartlett.test(centroidY ~ Category, data = all.dat.2) # significant

oneway.test(EOO.alpha ~ Category, data = all.dat.2, var.equal = FALSE)
# Perform the Games-Howell post hoc test
games_howell_test(data = all.dat.2, formula = EOO.alpha ~ Category)
# overdisper.shifter dif. specialists
# shifters dif. specialists
oneway.test(centroidX ~ Category, data = all.dat.2, var.equal = FALSE)
games_howell_test(data = all.dat.2, formula = centroidX ~ Category)
# generalists dif. shifting
# overdisp.shifter dif. shifting
# overdisp.shifter dif. specialists
oneway.test(centroidY ~ Category, data = all.dat.2, var.equal = FALSE)
games_howell_test(data = all.dat.2, formula = centroidY ~ Category)
# generalists dif. overdisper.shifting
# overdisp.shifter dif. shifting
# overdisp.shifter dif. specialists

#### Box plots for ESA ####

# read in range data
ranges = read.csv("./Formatted.Data/all.range.estimates.csv", row.names = 1)

# read in soil data
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

# merge these together to get ranges and quadrants into the same df
all.dat = inner_join(soil, ranges, by = c("species" = "Species"))

# EOO plot
EOO.plot = ggplot(all.dat, aes(x = Category, y = EOO.alpha, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#7C873E",
    "overdisper.shifter" = "#6E687E")) +
  theme_classic(base_size = 15) +
  labs(y = expression("Extent of Occurrence (km"^2*")")) +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdispersed\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none")
EOO.plot

ggsave("./Plots/ESA.plots/soil.EOO.png", width = 7, height = 5)

# CentroidX plot
centroidX.plot = ggplot(all.dat, aes(x = Category, y = centroidX, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#7C873E",
    "overdisper.shifter" = "#6E687E")) +
  theme_classic(base_size = 15) +
  labs(y = "Centroid of Longitudinal Range (°)") +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdisperse\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none")
centroidX.plot

ggsave("./Plots/ESA.plots/soil.centroidX.png", width = 7, height = 5)

# CentroidY plot
centroidY.plot = ggplot(all.dat, aes(x = Category, y = centroidY, color = Category)) +
  geom_boxplot() +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#7C873E",
    "overdisper.shifter" = "#6E687E")) +
  theme_classic(base_size = 15) +
  labs(y = "Centroid of Latitudinal Range (°)") +
  scale_x_discrete(name ="Category", 
                   labels=c("generalists" = "Generalists", 
                            "overdisper.shifter"= "Overdisperse\nShifters",
                            "overdisperser" = "Overdispersers", "shifting" = "Shifters",
                            "specialists" = "Specialists"))+
  theme(legend.position="none")
centroidY.plot

ggsave("./Plots/ESA.plots/soil.centroidY.png", width = 7, height = 5)

"#EAAE37"

#### Explanatory figure for ESA ####

# Create dummy data for line and confidence envelopes
x_vals <- seq(0, 1000, length.out = 100)

plot_data <- data.frame(
  geographic_distance = x_vals,
  fit = 0.1 + 0.0008 * x_vals,                       # Example line equation
  lwr_wide = 0.1 + 0.0008 * x_vals - 0.10,          # Wide envelope lower bound
  upr_wide = 0.1 + 0.0008 * x_vals + 0.10,          # Wide envelope upper bound
  lwr_tight = 0.1 + 0.0008 * x_vals - 0.05,         # Tight envelope lower bound
  upr_tight = 0.1 + 0.0008 * x_vals + 0.05          # Tight envelope upper bound
)

# Plot
ggplot() +
  geom_ribbon(data = plot_data, aes(x = geographic_distance, ymin = lwr_wide, ymax = upr_wide),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(data = plot_data, aes(x = geographic_distance, ymin = lwr_tight, ymax = upr_tight),
              fill = "gray90") +
  geom_line(data = plot_data, aes(x = geographic_distance, y = fit),
            color = "#5495CF", size = 2) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic(base_size = 15) +
  labs(x = "Geographic Distance (hm)", y = "Microclimate Distance", title = "SHIFTER")



#### Testing for phylogenetic signal in range size ####

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




