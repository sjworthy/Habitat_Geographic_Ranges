# plotting

all.clim = read.csv("./Results/climate.data.output.1.csv")

# get columns we need, scale microenvironment columns

all.clim.2 = all.clim %>%
  dplyr::select(species,decimalLongitude,decimalLatitude,high_temp_C,low_temp_C,moisture_mm) %>%
  mutate(high.temp.scaled = scale(high_temp_C),
         low.temp.scaled = scale(low_temp_C),
         moisture.scaled = scale(moisture_mm))

# Focus on Acer_leucoderme: 265 individuals

A_leu = all.clim.2 %>%
  filter(species == "Acer leucoderme")

# creating spatial matrix
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])

# creating microclim data
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))

# calculate Haversine distance for spatial data
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

# distance decay plot

A.leu.microclim.vec = as.vector(A.leu.microclim.dat)
A.leu.geo.vec = as.vector(A.leu.geo.dist.2)

A.leu.df = data.frame(A.leu.microclim.vec,A.leu.geo.vec)

A.leu.high.temp = ggplot(A.leu.df, aes(y = high_temp_C, x = A.leu.geo.vec)) + 
  #geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(x = "Geographical Distance", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.high.temp

ggsave(A.leu.high.temp, file = "./Plots/A.leu.high.temp.diss.pdf", height = 5, width = 5)

A.leu.cold.temp = ggplot(A.leu.df, aes(y = low_temp_C, x = A.leu.geo.vec)) + 
  #geom_point(size = 3, alpha = 0.5) + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(x = "Geographical Distance", y = "Cold Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.cold.temp

ggsave(A.leu.cold.temp, file = "./Plots/A.leu.cold.temp.diss.pdf", height = 5, width = 5)

A.leu.moist.temp = ggplot(A.leu.df, aes(y = moisture_mm, x = A.leu.geo.vec)) + 
  #geom_point(size = 3, alpha = 0.5) +
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(x = "Geographical Distance", y = "Moisture Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.moist.temp

ggsave(A.leu.cold.temp, file = "./Plots/A.leu.ppt.diss.pdf", height = 5, width = 5)

# distance decay plots from manuscript

plot(A.leu.geo.dist.2,A.leu.microclim.dist,
     xlab="Geographic distance (degrees)",ylab="Microclimate distance",
     pch=21,bg="grey")

lm1<-lm(A.leu.microclim.dist~I(A.leu.geo.dist.2))
summary(lm1)
plot(function(x){cbind(1,x)%*%coef(lm1)},xlim=c(0,1500000),add=TRUE,lwd=2,col=2)


A.leu.high.temp.dist = gowdis(as.matrix(A.leu.microclim.dat$high_temp_C))
A.leu.low.temp.dist = gowdis(as.matrix(A.leu.microclim.dat$low_temp_C))
A.leu.ppt.dist = gowdis(as.matrix(A.leu.microclim.dat$moisture_mm))

plot(A.leu.geo.dist.2,A.leu.high.temp.dist,
     xlab="Geographic distance (degrees)",ylab="High temp",
     pch=21,bg="grey")

lm1<-lm(A.leu.high.temp.dist~I(A.leu.geo.dist.2))
summary(lm1)
plot(function(x){cbind(1,x)%*%coef(lm1)},xlim=c(0,1500000),add=TRUE,lwd=2,col=2)

plot(A.leu.geo.dist.2,A.leu.low.temp.dist,
     xlab="Geographic distance (degrees)",ylab="Low temp",
     pch=21,bg="grey")

lm1<-lm(A.leu.low.temp.dist~I(A.leu.geo.dist.2))
summary(lm1)
plot(function(x){cbind(1,x)%*%coef(lm1)},xlim=c(0,1500000),add=TRUE,lwd=2,col=2)


plot(A.leu.geo.dist.2,A.leu.ppt.dist,
     xlab="Geographic distance (degrees)",ylab="Moisture temp",
     pch=21,bg="grey")

lm1<-lm(A.leu.ppt.dist~I(A.leu.geo.dist.2))
summary(lm1)
plot(function(x){cbind(1,x)%*%coef(lm1)},xlim=c(0,1500000),add=TRUE,lwd=2,col=2)



### PCA of slopes, intercepts, R2 from MRM microclimate models

MRM.clim = read.csv("./Results/Q1.MRM.results.csv")

# distribution of R2 values and intercept values
ggplot(MRM.clim, aes(y = Intercept, x = R2))+
  geom_point()
# distribution of R2 values and slope values
ggplot(MRM.clim, aes(y = Slope, x = R2))+
  geom_point()


# PCA of slope and R2
pc = princomp(MRM.clim[,c(4,5,7)], cor = TRUE)
summary(pc)
pc$loadings
# PC1: positively associated with Intercept and negatively associated with R2
# PC2: positively associated with Slope

library(ggbiplot)
ggbiplot(pc, labels = MRM.clim$species)+
  theme_classic()
