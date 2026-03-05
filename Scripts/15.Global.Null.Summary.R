# code to process and summarize global nulls.

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

mean(null_output$intercepts)
# 0.07125867
mean(null_output$slopes)
# 0.001053699
mean(null_output$R2)
# 0.3762577

ggplot(null_output, aes(intercepts))+
  geom_density()
ggplot(null_output, aes(slopes))+
  geom_density()
ggplot(null_output, aes(R2))+
  geom_density()

#### Geographic Distance Bins ####

range(df_plot$GeoDist)
# 0 to 3538.93

# Create quantile bins (5 bins) and convert to factor
df_plot$Geographic.Distance.Quantile <- cut(
  df_plot$GeoDist,
  breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0  524.6038  908.4138  1273.9590  1646.7984 3538.9299

Global.Null.Quant.Plot = ggplot(df_plot, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–524",
                              "2" = "525–908",
                              "3" = "909–1273",
                              "4" = "1274–1646",
                              "5" = "1647–3538")) +
  labs(y = "Microclimate Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Global.Null.Quant.Plot

ggsave(Global.Null.Quant.Plot, file = "./Plots/Microclim.Global.Null.Quantile.pdf", height = 5, width = 5)

#### Topo ####

#### Geographic Distance Bins ####

range(df_plot$GeoDist)
# 0 to 3538.93

# Create quantile bins (5 bins) and convert to factor
df_plot$Geographic.Distance.Quantile <- cut(
  df_plot$GeoDist,
  breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(df_plot$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0  524.6038  908.4138  1273.9590  1646.7984 3538.9299

Global.Null.Quant.Plot = ggplot(df_plot, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–524",
                              "2" = "525–908",
                              "3" = "909–1273",
                              "4" = "1274–1646",
                              "5" = "1647–3538")) +
  labs(y = "Topographic Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Global.Null.Quant.Plot

ggsave(Global.Null.Quant.Plot, file = "./Plots/Topo.Global.Null.Quantile.pdf", height = 5, width = 5)

### Soil ###
# read in data for first plot
# just picked 500 because it has the largest range
dat = read.csv("./Results/Soil.Global.Null/soil.global.dist_ 500 .csv", row.names = 1)

range(dat$GeoDist)
# 0.002836702 to 346.7406

# Create quantile bins (5 bins) and convert to factor
dat$Geographic.Distance.Quantile <- cut(
  dat$GeoDist,
  breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(dat$GeoDist, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#       0%       20%       40%       60%       80%      100% 
# 0.002836702  52.84461  89.90785  124.9091  163.797 346.7406

Global.Null.Quant.Plot = ggplot(dat, aes(x = as.factor(Geographic.Distance.Quantile), y = SoilDist), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  scale_x_discrete(labels = c("1" = "0–53",
                              "2" = "54–90",
                              "3" = "91–124",
                              "4" = "125–163",
                              "5" = "164–346")) +
  labs(y = "Soil Distance",
       x = "Geographic Distance Quantile (hm)", fill = " ") +
  theme_classic()
Global.Null.Quant.Plot

ggsave(Global.Null.Quant.Plot, file = "./Plots/Soil.Global.Null.Quantile.png", height = 5, width = 5)


# Plot the data regular as well

plot_obj = ggplot(dat, aes(x = GeoDist, y = SoilDist)) +
  geom_point(shape = 21, fill = "grey", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  ylim(0,1)+
  xlab("Geographic Distance (hm)") + # should be hectometer
  ylab("Soil Distance") +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
plot_obj

# Save the plot
ggsave(plot_obj, file = "./Plots/Soil.Global.Null.png", height = 5, width = 5)









