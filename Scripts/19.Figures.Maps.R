# Code to generate map figures

library(maps)
library(tigris)
options(tigris_use_cache = TRUE)
library(sf)
library(BIEN)
library(viridis)
library(cowplot)
library(elevatr)
library(terra)
library(patchwork)

#### Macro minus Micro ####

# read in data
data = read.csv("./Formatted.Data/All.Final.Data.csv", row.names = 1)

# convert to sf object
points_sf = sf::st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) # WG 84

# calculate the difference between macro and micro
points_sf$ppt.diff = points_sf$macro_bio12_total_annual_precip_mm - points_sf$moisture_mm
points_sf$high.temp.diff = points_sf$macro_bio5_max_temp_warm_month_C - points_sf$high_temp_C
points_sf$low.temp.diff = points_sf$macro_bi06_min_temp_cold_month_C - points_sf$low_temp_C

# Get ranges/limits for plotting
ppt_lim_max = round(max(points_sf$ppt.diff, na.rm = TRUE))
ppt_lim_min = round(min(points_sf$ppt.diff, na.rm = TRUE))
high_temp_lim_max = round(max(points_sf$high.temp.diff, na.rm = TRUE))
high_temp_lim_min = round(min(points_sf$high.temp.diff, na.rm = TRUE))
low_temp_lim_max = round(max(points_sf$low.temp.diff, na.rm = TRUE))
low_temp_lim_min = round(min(points_sf$low.temp.diff, na.rm = TRUE))

# get base map of US states
maps::map(database = "state")
us_states = states(cb = TRUE)
continental_states = us_states %>%
  filter(!NAME %in% (c("Alaska","American Samoa","Guam","Commonwealth of the Northern Mariana Islands","Hawaii","United States Virgin Islands",
                       "Puerto Rico")))
states.map = continental_states %>%
  st_as_sf %>%
  st_transform(st_crs(points_sf))

# Make bounding box to keep everything east of -113.5
bbox = st_bbox(c(
  xmin = -113.5,
  xmax = -66,
  ymin = 24,
  ymax = 50), crs = st_crs(states.map))

# Crop the states map based on bounding box
states.map.east = st_crop(states.map, bbox)

# Get elevation raster for US map
dem = get_elev_raster(locations = states.map, 
                      z = 5, # resolution (higher = more detail, slower)
                      clip = "locations")

dem = rast(dem)  # convert to terra raster
dem = project(dem, crs(states.map)) # match projections
dem_crop = crop(dem, states.map.east) # crop to east states map
contours = as.contour(dem_crop, levels = seq(0, 3000, by = 200)) # create contour lines
contours_sf = st_as_sf(contours) # convert to sf for plotting

# get species list
species_list = unique(points_sf$species)

# open PDF to save all plots
pdf("./Plots/Species_Plots.pdf", width = 10, height = 12)

for(sp in species_list) {
  
  # subset species
  sp_data <- points_sf[points_sf$species == sp, ]
  
  ppt.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = sp_data, aes(color = ppt.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "Precipitation Difference (mm)",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 6),
        legend.spacing.y = unit(0.1, "cm"))+
  ggtitle(paste("Precipitation Difference -", sp))

high.temp.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = sp_data, aes(color = high.temp.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "Temperature Difference (°C)",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 6),
        legend.spacing.y = unit(0.1, "cm"))+
  ggtitle(paste("High Temperature Difference -", sp))

low.temp.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = sp_data, aes(color = low.temp.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "Temperature Difference  (°C)",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 6),
        legend.spacing.y = unit(0.1, "cm"))+
  ggtitle(paste("Low Temperature Difference -", sp))

# combine plots (one page per species)
combined_plot = ppt.diff.plot / high.temp.diff.plot / low.temp.diff.plot

print(combined_plot)
}

# close PDF
dev.off()

##### Map of Distances ####

cercis.data = read.csv("./Mean_Dist/Mean.Dist.Cercis canadensis.csv")
cercis_sf <- sf::st_as_sf(cercis.data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) # WG 84

micro.clim.dist.plot = ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = cercis_sf, aes(color = microclim.dist.mean), size = 2) +
  scale_color_viridis(name = "Mean Microclimate Distance", option = "C", limits = c(0,1)) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Mean Microclimate Distance")
micro.clim.dist.plot

ggsave("./Plots/Cercis.microclim.dist.limits.png", width = 8, height = 6)
ggsave("./Plots/Cercis.microclim.dist.png", width = 8, height = 6)

topo.dist.plot = ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = cercis_sf, aes(color = topo.dist.mean), size = 2) +
  scale_color_viridis(name = "Mean Topography Distance", option = "C", limits = c(0,1)) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Mean Topography Distance")
topo.dist.plot

ggsave("./Plots/Cercis.topo.dist.limits.png", width = 8, height = 6)
ggsave("./Plots/Cercis.topo.dist.png", width = 8, height = 6)

soil.dist.plot = ggplot()+
  geom_sf(data = states.map, fill = "white")+
  geom_sf(data = cercis_sf, aes(color = soil.dist.mean), size = 2) +
  scale_color_viridis(name = "Mean Soil Distance", option = "C") +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Mean Soil Distance")
soil.dist.plot

ggsave("./Plots/Cercis.soil.dist.limits.png", width = 8, height = 6)
ggsave("./Plots/Cercis.soil.dist.png", width = 8, height = 6)







range(points_sf$ppt.diff) # -335.6558  244.8182
range(points_sf$high.temp.diff) # -1.410587  1.444130
range(points_sf$low.temp.diff) # -1.739132  2.456522


ppt.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = points_sf, aes(color = ppt.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "Precipitation Difference",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Difference between Macro- and Micro- Precipitation")
ppt.diff.plot

high.temp.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = points_sf, aes(color = high.temp.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "High Temperature Difference",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Difference between Macro- and Micro- High Temperature")
high.temp.diff.plot

low.temp.diff.plot = ggplot()+
  geom_sf(data = states.map.east, fill = "white")+
  geom_sf(data = points_sf, aes(color = low.temp.diff), size = 0.5) +
  geom_sf(data = contours_sf,color = "grey40", size = 0.2, alpha = 0.6)+
  scale_color_gradient2(
    name = "Low Temperature Difference",
    low = "#3B6FB6",     # negative values
    mid = "#F0F0F0",   # zero
    high = "#C23B33",     # positive values
    midpoint = 0)+
  coord_sf(xlim = c(-115, -66), ylim = c(24, 50), expand = FALSE) +
  theme_classic()+
  ggtitle("Difference between Macro- and Micro- Low Temperature")
low.temp.diff.plot


