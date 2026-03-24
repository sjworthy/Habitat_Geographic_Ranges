# Code to generate figures 

library(vegan)
library(ggtree)
library(tidyverse)

### NMDS of slopes, intercepts, R2 from MRM microclim models ####

# read in soil data with categories
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)

microclim[123,c(2,4,6)] = c(0.08361417,0.001727047,0.3914622)
microclim[123,1] = "Null"
microclim[123,19] = "Null"
microclim[123,20] = "significant"

# data for NMDS
microclim.2 = microclim[,c(2,4,6)]

microclim.nmds = metaMDS(microclim.2, distance = "bray")
microclim.nmds
# stress = 0.03166897, this is good

# plotting
nmds_scores <- as.data.frame(scores(microclim.nmds)$sites)
nmds_scores$Category <- microclim$Category
nmds_scores$shape = microclim$significant

# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
microclim.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Category, fill = Category, shape = shape),
             size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("significant" = 16, "non-significant" = 1),
                     name = "Strength",
                     labels = c("significant" = "Significant", "non-significant" = "Non-significant")) +
  scale_fill_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  theme_classic() +
  stat_ellipse(aes(color = Category), linewidth = 1)+
  ggtitle("Microclimate")+
  theme(legend.position = "none")

microclim.nmds

#ggsave("./Plots/microclim.MNDS.ellipses.legend.png", width = 5, height = 5)
#ggsave("./Plots/microclim.MNDS.ellipses.png", width = 5, height = 5)

### NMDS of slopes, intercepts, R2 from MRM topography models ####

# read in soil data with categories
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)

topo[123,c(2,4,6)] = c(0.07670251,0.0002537613,0.03900371)
topo[123,1] = "Null"
topo[123,19] = "Null"
topo[123,20] = "significant"

# data for NMDS
topo.2 = topo[,c(2,4,6)]

# some values are negative so transform to all positive
# Shift all values to be positive by adding the absolute minimum
topo.trans <- topo.2 + abs(min(topo.2))

topo.nmds = metaMDS(topo.trans, distance = "bray")
topo.nmds
# stress = 0.05739008, this is good

# plotting
nmds_scores <- as.data.frame(scores(topo.nmds)$sites)
nmds_scores$Category <- topo$Category
nmds_scores$shape = topo$significant

# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
topo.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Category, fill = Category, shape = shape),
             size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("significant" = 16, "non-significant" = 1),
                     name = "Strong",
                     labels = c("significant" = "Significant", "non-significant" = "Non-significant")) +
  scale_fill_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  theme_classic() +
  stat_ellipse(aes(color = Category), linewidth = 1)+
  ggtitle("Topography")+
  theme(legend.position = "none")

topo.nmds

ggsave("./Plots/topo.MNDS.ellipses.png", width = 5, height = 5)

### NMDS of slopes, intercepts, R2 from MRM soil models ####

# read in soil data with categories
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

soil[123,c(2,4,6)] = c(0.1948436,0.0009679467,0.101624)
soil[123,1] = "Null"
soil[123,19] = "Null"
soil[123,20] = "significant"

# data for NMDS
soil.2 = soil[,c(2,4,6)]

soil.nmds = metaMDS(soil.2, distance = "bray")
soil.nmds
# stress = 0.05266094, this is good

# plotting
nmds_scores <- as.data.frame(scores(soil.nmds)$sites)
nmds_scores$Category <- soil$Category
nmds_scores$shape = soil$significant

# colors
# shifter = "#5495CF",
# specialist = "#DB4743",
# generalist = "#F5AF4D",
# overdisperser = "#548F01",
# overdisperser shifter ="#B46DB3"
# Null = "black"

# Plot with shape mapped to combined variable  
soil.nmds = ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Category, fill = Category, shape = shape),
             size = 3, stroke = 1.2) +
  scale_shape_manual(values = c("significant" = 16, "non-significant" = 1),
                     name = "Strong",
                     labels = c("significant" = "Significant", "non-significant" = "Non-significant")) +
  scale_fill_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  scale_color_manual(values = c(
    "shifting" = "#5495CF",
    "specialists" = "#DB4743",
    "generalists" = "#F5AF4D",
    "overdisperser" = "#548F01",
    "overdisper.shifter" = "#B46DB3",
    "Null" = "black"),
    name = "Pattern",
    labels = c(
      "shifting" = "Shifter",
      "specialists" = "Specialist",
      "generalists" = "Generalist",
      "overdisperser" = "Overdisperser",
      "overdisper.shifter" = "Overdispersed Shifter",
      "Null" = "Null"
    )) +
  theme_classic() +
  stat_ellipse(aes(color = Category), linewidth = 1)+
  ggtitle("Soil")+
  theme(legend.position = "none")

soil.nmds

# ggsave("./Plots/soil.MNDS.ellipses.png", width = 5, height = 5)

### Plotting Categories on the phylogeny ####

phylo <- ggtree::read.tree("./Results/phylo.tre")
phylo$tip.label = gsub("_", " ", phylo$tip.label)

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo) +
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree()
tree

# read in the category data
microclim = read.csv("./Results/microclim.global.null.compare.results.csv", row.names = 1)
topo = read.csv("./Results/topo.global.null.compare.results.csv", row.names = 1)
soil = read.csv("./Results/soil.global.null.compare.results.csv", row.names = 1)

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
colnames(all.cats) = c("Microclimate","Topography","Soil")

write.csv(all.cats, file = "./Results/occupancy.patterns.csv")

cat.phylo = gheatmap(tree, all.cats, offset=22, width=0.4, 
                     colnames=TRUE, colnames_position = "top",
                     colnames_offset_y = 0.5, legend_title="Category") +
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
  #theme(legend.position = "none")
cat.phylo

#ggsave("./Plots/phylo.png", height = 12, width = 12)
#ggsave("./Plots/phylo.legend.png", height = 12, width = 12)

#### Explanatory Figure ####
# all patterns on one plot

set.seed(123)

# Simulate realistic data
gen_x_vals <- seq(0, 1000, length.out = 100)
gen_intercept <- 0.27      
gen_slope <- 0.00024       
gen_noise <- rnorm(100, mean = 0, sd = 0.02)
gen_y_vals <- gen_intercept + gen_slope * gen_x_vals + gen_noise
gen_df <- data.frame(x = gen_x_vals, y = gen_y_vals)

# Simulate realistic data
ods_x_vals <- seq(0, 1000, length.out = 100)
ods_intercept <- 0.75      # High
ods_slope <- 0.0003        # high slope
ods_noise <- rnorm(100, mean = 0, sd = 0.02)
ods_y_vals <- ods_intercept + ods_slope * ods_x_vals + ods_noise
ods_df <- data.frame(x = ods_x_vals, y = ods_y_vals)

od_x_vals <- seq(0, 1000, length.out = 100)
od_intercept <- 0.75      # High
od_slope <- 0.00001        # Lower slope than before
od_noise <- rnorm(100, mean = 0, sd = 0.02)
od_y_vals <- od_intercept + od_slope * od_x_vals + od_noise
od_df <- data.frame(x = od_x_vals, y = od_y_vals)

spec_x_vals <- seq(0, 1000, length.out = 100)
spec_intercept <- 0.05      # Still low
spec_slope <- 0.00001        # Lower slope than before
spec_noise <- rnorm(100, mean = 0, sd = 0.02)
spec_y_vals <- spec_intercept + spec_slope * spec_x_vals + spec_noise
spec_df <- data.frame(x = spec_x_vals, y = spec_y_vals)

shift_x_vals <- seq(0, 1000, length.out = 100)
shift_intercept <- 0.05      # Low intercept
shift_slope <- 0.0009        # High slope (relative to small y-scale)
shift_noise <- rnorm(100, mean = 0, sd = 0.02)
shift_y_vals <- shift_intercept + shift_slope * shift_x_vals + shift_noise
shift_df <- data.frame(x = shift_x_vals, y = shift_y_vals)

# Plot
all.plot <- ggplot() +
  geom_smooth(data = gen_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#F5AF4D", size = 2) +
  geom_smooth(data = ods_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#B46DB3", size = 2) +
  geom_smooth(data = od_df, aes(x = x, y = y), method = "lm", se = FALSE, color = "#548F01", size = 2)+ 
  geom_smooth(data = spec_df, aes(x = x, y = y),method = "lm", se = FALSE, color = "#DB4743", size = 2) +
  geom_smooth(data = shift_df, aes(x = x, y = y),method = "lm", se = FALSE, color = "#5495CF", size = 2) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1)) +
  theme_classic(base_size = 22) +
  labs(
    x = "Geographic Distance (km)",
    y = "Envirionmental Dissimilarity")

all.plot

ggsave("./Plots/all.patterns.png", width = 5, height = 5)

#### Old ESA figure ####

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








