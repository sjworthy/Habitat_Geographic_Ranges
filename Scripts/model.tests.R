# Example Analyses with select species

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)

# species 1: Quercus coccinea, R2 = 0.000289395, n = 641
# species 2: Acer leucoderme, R2 = 0.23, n = 265
# species 3: Fraxinus caroliniana, R2 = 0.50, n = 558
# species 4: Carya texana R2 = 0.76, n = 725

#### Get each species' data set ready ####

all.clim = read.csv("./Results/climate.data.output.1.csv")

Q.cocc = all.clim %>%
  filter(species == "Quercus coccinea")
A.leu = all.clim %>%
  filter(species == "Acer leucoderme")
F.car = all.clim %>% 
  filter(species == "Fraxinus caroliniana")
C.tex = all.clim %>%
  filter(species == "Carya texana")

# read in and merge with elevation data

Q.cocc.elev = read.csv("./Formatted.Data/elev/all_data_Quercus coccinea.csv")
A.leu.elev = read.csv("./Formatted.Data/elev/all_data_Acer leucoderme.csv")
F.car.elev = read.csv("./Formatted.Data/elev/all_data_Fraxinus caroliniana.csv")
C.text.elev = read.csv("./Formatted.Data/elev/all_data_Carya texana.csv")

Q.cocc.all = left_join(Q.cocc,Q.cocc.elev)
A.leu.all = left_join(A.leu,A.leu.elev)
F.car.all = left_join(F.car,F.car.elev)
C.tex.all = left_join(C.tex,C.text.elev)

#### Microclimate distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating microclimate data
Q.cocc.microclim.dat = Q.cocc.all %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A.leu.all %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F.car.all %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
C.tex.microclim.dat = C.tex.all %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
Q.cocc.microclim.dist = gowdis(as.matrix(Q.cocc.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))
C.tex.microclim.dist = gowdis(as.matrix(C.tex.microclim.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.clim.mod = MRM(Q.cocc.microclim.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.clim.mod, file = "./Results/test.results/Q.cocc.clim.mod.RDS")
A.leu.clim.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.clim.mod, file = "./Results/test.results/A.leu.clim.mod.RDS")
F.car.clim.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.3)
saveRDS(F.car.clim.mod, file = "./Results/test.results/F.car.clim.mod.RDS")
C.tex.clim.mod = MRM(C.tex.microclim.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.clim.mod, file = "./Results/test.results/C.tex.clim.mod.RDS")

#### Microclimate distance decay plots ####

Q.cocc.microclim.vec = as.vector(Q.cocc.microclim.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.microclim.vec,Q.cocc.geo.vec)

Q.cocc.microclim = ggplot(Q.cocc.df, aes(y = Q.cocc.microclim.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0", x = "Geographical Distance (km)", y = "Microclimate Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.microclim

ggsave(Q.cocc.microclim, file = "./Results/test.results/Plots/Q.cocc.microclim.pdf", height = 5, width = 5)

A.leu.microclim.vec = as.vector(A.leu.microclim.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.microclim.vec,A.leu.geo.vec)

A.leu.microclim = ggplot(A.leu.df, aes(y = A.leu.microclim.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.23", x = "Geographical Distance (km)", y = "Microclimate Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.microclim

ggsave(A.leu.microclim, file = "./Results/test.results/Plots/A.leu.microclim.pdf", height = 5, width = 5)

F.car.microclim.vec = as.vector(F.car.microclim.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.microclim.vec,F.car.geo.vec)

F.car.microclim = ggplot(F.car.df, aes(y = F.car.microclim.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.50", x = "Geographical Distance (km)", y = "Microclimate Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.microclim

ggsave(F.car.microclim, file = "./Results/test.results/Plots/F.car.microclim.pdf", height = 5, width = 5)

C.tex.microclim.vec = as.vector(C.tex.microclim.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.microclim.vec,C.tex.geo.vec)

C.tex.microclim = ggplot(C.tex.df, aes(y = C.tex.microclim.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.76", x = "Geographical Distance (km)", y = "Microclimate Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.microclim

ggsave(C.tex.microclim, file = "./Results/test.results/Plots/C.tex.microclim.pdf", height = 5, width = 5)




#### High temp distance ~ Geographic Distance ####
# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating microclimate data
Q.cocc.microclim.dat = Q.cocc.all %>%
  dplyr::select(high_temp_C)
A.leu.microclim.dat = A.leu.all %>%
  dplyr::select(high_temp_C)
F.car.microclim.dat = F.car.all %>%
  dplyr::select(high_temp_C)
C.tex.microclim.dat = C.tex.all %>%
  dplyr::select(high_temp_C)

# calculate gower distance for scaled microclimate data
Q.cocc.microclim.dist = gowdis(as.matrix(Q.cocc.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))
C.tex.microclim.dist = gowdis(as.matrix(C.tex.microclim.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.high.temp.mod = MRM(Q.cocc.microclim.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.high.temp.mod, file = "./Results/test.results/Q.cocc.high.temp.mod.RDS")
A.leu.high.temp.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.high.temp.mod, file = "./Results/test.results/A.leu.high.temp.mod.RDS")
F.car.high.temp.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.3)
saveRDS(F.car.high.temp.mod, file = "./Results/test.results/F.car.high.temp.mod.RDS")
C.tex.high.temp.mod = MRM(C.tex.microclim.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.high.temp.mod, file = "./Results/test.results/C.tex.high.temp.mod.RDS")

#### High temp distance decay plots ####

Q.cocc.microclim.vec = as.vector(Q.cocc.microclim.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.microclim.vec,Q.cocc.geo.vec)

Q.cocc.microclim = ggplot(Q.cocc.df, aes(y = Q.cocc.microclim.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0", x = "Geographical Distance (km)", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.microclim

ggsave(Q.cocc.microclim, file = "./Results/test.results/Plots/Q.cocc.high.temp.pdf", height = 5, width = 5)

A.leu.microclim.vec = as.vector(A.leu.microclim.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.microclim.vec,A.leu.geo.vec)

A.leu.microclim = ggplot(A.leu.df, aes(y = A.leu.microclim.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.08", x = "Geographical Distance (km)", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.microclim

ggsave(A.leu.microclim, file = "./Results/test.results/Plots/A.leu.high.temp.pdf", height = 5, width = 5)

F.car.microclim.vec = as.vector(F.car.microclim.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.microclim.vec,F.car.geo.vec)

F.car.microclim = ggplot(F.car.df, aes(y = F.car.microclim.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.15", x = "Geographical Distance (km)", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.microclim

ggsave(F.car.microclim, file = "./Results/test.results/Plots/F.car.high.temp.pdf", height = 5, width = 5)

C.tex.microclim.vec = as.vector(C.tex.microclim.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.microclim.vec,C.tex.geo.vec)

C.tex.microclim = ggplot(C.tex.df, aes(y = C.tex.microclim.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.39", x = "Geographical Distance (km)", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.microclim

ggsave(C.tex.microclim, file = "./Results/test.results/Plots/C.tex.high.temp.pdf", height = 5, width = 5)




#### Low temp distance ~ Geographic Distance ####
# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating microclimate data
Q.cocc.microclim.dat = Q.cocc.all %>%
  dplyr::select(low_temp_C)
A.leu.microclim.dat = A.leu.all %>%
  dplyr::select(low_temp_C)
F.car.microclim.dat = F.car.all %>%
  dplyr::select(low_temp_C)
C.tex.microclim.dat = C.tex.all %>%
  dplyr::select(low_temp_C)

# calculate gower distance for scaled microclimate data
Q.cocc.microclim.dist = gowdis(as.matrix(Q.cocc.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))
C.tex.microclim.dist = gowdis(as.matrix(C.tex.microclim.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.low.temp.mod = MRM(Q.cocc.microclim.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.low.temp.mod, file = "./Results/test.results/Q.cocc.low.temp.mod.RDS")
A.leu.low.temp.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.low.temp.mod, file = "./Results/test.results/A.leu.low.temp.mod.RDS")
F.car.low.temp.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.3)
saveRDS(F.car.low.temp.mod, file = "./Results/test.results/F.car.low.temp.mod.RDS")
C.tex.low.temp.mod = MRM(C.tex.microclim.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.low.temp.mod, file = "./Results/test.results/C.tex.low.temp.mod.RDS")


#### Low temp distance decay plots ####

Q.cocc.microclim.vec = as.vector(Q.cocc.microclim.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.microclim.vec,Q.cocc.geo.vec)

Q.cocc.microclim = ggplot(Q.cocc.df, aes(y = Q.cocc.microclim.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.03", x = "Geographical Distance (km)", y = "low Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.microclim

ggsave(Q.cocc.microclim, file = "./Results/test.results/Plots/Q.cocc.low.temp.pdf", height = 5, width = 5)

A.leu.microclim.vec = as.vector(A.leu.microclim.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.microclim.vec,A.leu.geo.vec)

A.leu.microclim = ggplot(A.leu.df, aes(y = A.leu.microclim.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.13", x = "Geographical Distance (km)", y = "low Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.microclim

ggsave(A.leu.microclim, file = "./Results/test.results/Plots/A.leu.low.temp.pdf", height = 5, width = 5)

F.car.microclim.vec = as.vector(F.car.microclim.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.microclim.vec,F.car.geo.vec)

F.car.microclim = ggplot(F.car.df, aes(y = F.car.microclim.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.36", x = "Geographical Distance (km)", y = "low Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.microclim

ggsave(F.car.microclim, file = "./Results/test.results/Plots/F.car.low.temp.pdf", height = 5, width = 5)

C.tex.microclim.vec = as.vector(C.tex.microclim.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.microclim.vec,C.tex.geo.vec)

C.tex.microclim = ggplot(C.tex.df, aes(y = C.tex.microclim.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.74", x = "Geographical Distance (km)", y = "low Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.microclim

ggsave(C.tex.microclim, file = "./Results/test.results/Plots/C.tex.low.temp.pdf", height = 5, width = 5)

#### PPT distance ~ Geographic Distance ####
# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating microclimate data
Q.cocc.microclim.dat = Q.cocc.all %>%
  dplyr::select(moisture_mm)
A.leu.microclim.dat = A.leu.all %>%
  dplyr::select(moisture_mm)
F.car.microclim.dat = F.car.all %>%
  dplyr::select(moisture_mm)
C.tex.microclim.dat = C.tex.all %>%
  dplyr::select(moisture_mm)

# calculate gower distance for scaled microclimate data
Q.cocc.microclim.dist = gowdis(as.matrix(Q.cocc.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))
C.tex.microclim.dist = gowdis(as.matrix(C.tex.microclim.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.ppt.mod = MRM(Q.cocc.microclim.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.ppt.mod, file = "./Results/test.results/Q.cocc.ppt.mod.RDS")
A.leu.ppt.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.ppt.mod, file = "./Results/test.results/A.leu.ppt.mod.RDS")
F.car.ppt.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.3)
saveRDS(F.car.ppt.mod, file = "./Results/test.results/F.car.ppt.mod.RDS")
C.tex.ppt.mod = MRM(C.tex.microclim.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.ppt.mod, file = "./Results/test.results/C.tex.ppt.mod.RDS")

#### PPT distance decay plots ####

Q.cocc.microclim.vec = as.vector(Q.cocc.microclim.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.microclim.vec,Q.cocc.geo.vec)

Q.cocc.microclim = ggplot(Q.cocc.df, aes(y = Q.cocc.microclim.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0", x = "Geographical Distance (km)", y = "ppt Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.microclim

ggsave(Q.cocc.microclim, file = "./Results/test.results/Plots/Q.cocc.ppt.pdf", height = 5, width = 5)

A.leu.microclim.vec = as.vector(A.leu.microclim.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.microclim.vec,A.leu.geo.vec)

A.leu.microclim = ggplot(A.leu.df, aes(y = A.leu.microclim.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.11", x = "Geographical Distance (km)", y = "ppt Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.microclim

ggsave(A.leu.microclim, file = "./Results/test.results/Plots/A.leu.ppt.pdf", height = 5, width = 5)

F.car.microclim.vec = as.vector(F.car.microclim.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.microclim.vec,F.car.geo.vec)

F.car.microclim = ggplot(F.car.df, aes(y = F.car.microclim.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.08", x = "Geographical Distance (km)", y = "PPT Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.microclim

ggsave(F.car.microclim, file = "./Results/test.results/Plots/F.car.ppt.pdf", height = 5, width = 5)

C.tex.microclim.vec = as.vector(C.tex.microclim.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.microclim.vec,C.tex.geo.vec)

C.tex.microclim = ggplot(C.tex.df, aes(y = C.tex.microclim.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.16", x = "Geographical Distance (km)", y = "PPT Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.microclim

ggsave(C.tex.microclim, file = "./Results/test.results/Plots/C.tex.ppt.pdf", height = 5, width = 5)


#### Topo distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(northness,eastness,mTPI,slope,elevation)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(northness,eastness,mTPI,slope,elevation)
F.car.topo.dat = F.car.all %>%
  dplyr::select(northness,eastness,mTPI,slope,elevation)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(northness,eastness,mTPI,slope,elevation)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.topo.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.topo.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.topo.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.topo.mod.RDS")

#### Topo distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.02", x = "Geographical Distance (km)", y = "Topo Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.topo.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.007", x = "Geographical Distance (km)", y = "Topo Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.topo.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.004", x = "Geographical Distance (km)", y = "Topo Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.topo.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.02", x = "Geographical Distance (km)", y = "Topo Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.topo.pdf", height = 5, width = 5)

#### Northness distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(northness)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(northness)
F.car.topo.dat = F.car.all %>%
  dplyr::select(northness)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(northness)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.northness.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.northness.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.northness.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.northness.mod.RDS")

#### Northness distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.005", x = "Geographical Distance (km)", y = "Northness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.northness.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.002", x = "Geographical Distance (km)", y = "Northness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.northness.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.0001", x = "Geographical Distance (km)", y = "Northness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.northness.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.00002", x = "Geographical Distance (km)", y = "Northness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.northness.pdf", height = 5, width = 5)

#### eastness distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(eastness)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(eastness)
F.car.topo.dat = F.car.all %>%
  dplyr::select(eastness)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(eastness)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.eastness.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.eastness.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.eastness.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.eastness.mod.RDS")

#### Eastness distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.006", x = "Geographical Distance (km)", y = "eastness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.eastness.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.00009", x = "Geographical Distance (km)", y = "eastness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.eastness.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.000009", x = "Geographical Distance (km)", y = "eastness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.eastness.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.003", x = "Geographical Distance (km)", y = "eastness Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.eastness.pdf", height = 5, width = 5)

#### mTPI distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(mTPI)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(mTPI)
F.car.topo.dat = F.car.all %>%
  dplyr::select(mTPI)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(mTPI)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.mTPI.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.mTPI.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.mTPI.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.mTPI.mod.RDS")


#### mTPI distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.02", x = "Geographical Distance (km)", y = "mTPI Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.mTPI.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.004", x = "Geographical Distance (km)", y = "mTPI Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.mTPI.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.001", x = "Geographical Distance (km)", y = "mTPI Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.mTPI.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.005", x = "Geographical Distance (km)", y = "mTPI Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.mTPI.pdf", height = 5, width = 5)

#### slope distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(slope)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(slope)
F.car.topo.dat = F.car.all %>%
  dplyr::select(slope)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(slope)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.slope.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.slope.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.slope.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.slope.mod.RDS")

#### slope distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.002", x = "Geographical Distance (km)", y = "slope Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.slope.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.002", x = "Geographical Distance (km)", y = "slope Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.slope.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.000001", x = "Geographical Distance (km)", y = "slope Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.slope.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.005", x = "Geographical Distance (km)", y = "slope Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.slope.pdf", height = 5, width = 5)

#### elevation distance ~ Geographic Distance ####

# creating spatial matrix
Q.cocc.spat.dat = as.matrix(Q.cocc.all[,c(3,4)])
A.leu.spat.dat = as.matrix(A.leu.all[,c(3,4)])
F.car.spat.dat = as.matrix(F.car.all[,c(3,4)])
C.tex.spat.dat = as.matrix(C.tex.all[,c(3,4)])

# creating topographic data
Q.cocc.topo.dat = Q.cocc.all %>%
  dplyr::select(elevation)
A.leu.topo.dat = A.leu.all %>%
  dplyr::select(elevation)
F.car.topo.dat = F.car.all %>%
  dplyr::select(elevation)
C.tex.topo.dat = C.tex.all %>%
  dplyr::select(elevation)

# calculate gower distance for scaled topoate data
Q.cocc.topo.dist = gowdis(as.matrix(Q.cocc.topo.dat))
A.leu.topo.dist = gowdis(as.matrix(A.leu.topo.dat))
F.car.topo.dist = gowdis(as.matrix(F.car.topo.dat))
C.tex.topo.dist = gowdis(as.matrix(C.tex.topo.dat))

# calculate Haversine distance for spatial data
Q.cocc.geo.dist = distm(Q.cocc.spat.dat, fun = distHaversine)
Q.cocc.geo.dist.2 = as.dist(Q.cocc.geo.dist) # convert to dist object
Q.cocc.geo.dist.3 = Q.cocc.geo.dist.2/1000 # convert to km
A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object
A.leu.geo.dist.3 = A.leu.geo.dist.2/1000 # convert to km
F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object
F.car.geo.dist.3 = F.car.geo.dist.2/1000 # convert to km
C.tex.geo.dist = distm(C.tex.spat.dat, fun = distHaversine)
C.tex.geo.dist.2 = as.dist(C.tex.geo.dist) # convert to dist object
C.tex.geo.dist.3 = C.tex.geo.dist.2/1000 # convert to km

# Perform MRM
Q.cocc.topo.mod = MRM(Q.cocc.topo.dist ~ Q.cocc.geo.dist.3)
saveRDS(Q.cocc.topo.mod, file = "./Results/test.results/Q.cocc.elevation.mod.RDS")
A.leu.topo.mod = MRM(A.leu.topo.dist ~ A.leu.geo.dist.3)
saveRDS(A.leu.topo.mod, file = "./Results/test.results/A.leu.elevation.mod.RDS")
F.car.topo.mod = MRM(F.car.topo.dist ~ F.car.geo.dist.3)
saveRDS(F.car.topo.mod, file = "./Results/test.results/F.car.elevation.mod.RDS")
C.tex.topo.mod = MRM(C.tex.topo.dist ~ C.tex.geo.dist.3)
saveRDS(C.tex.topo.mod, file = "./Results/test.results/C.tex.elevation.mod.RDS")

#### elevation distance decay plots ####

Q.cocc.topo.vec = as.vector(Q.cocc.topo.dist)
Q.cocc.geo.vec = as.vector(Q.cocc.geo.dist.3)

Q.cocc.df = data.frame(Q.cocc.topo.vec,Q.cocc.geo.vec)

Q.cocc.topo = ggplot(Q.cocc.df, aes(y = Q.cocc.topo.vec, x = Q.cocc.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Quercus coccinea R2 = 0.007", x = "Geographical Distance (km)", y = "elevation Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
Q.cocc.topo

ggsave(Q.cocc.topo, file = "./Results/test.results/Plots/Q.cocc.elevation.pdf", height = 5, width = 5)

A.leu.topo.vec = as.vector(A.leu.topo.dist)
A.leu.geo.vec = as.vector(A.leu.geo.dist.3)

A.leu.df = data.frame(A.leu.topo.vec,A.leu.geo.vec)

A.leu.topo = ggplot(A.leu.df, aes(y = A.leu.topo.vec, x = A.leu.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Acer leucoderme R2 = 0.02", x = "Geographical Distance (km)", y = "elevation Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.leu.topo

ggsave(A.leu.topo, file = "./Results/test.results/Plots/A.leu.elevation.pdf", height = 5, width = 5)

F.car.topo.vec = as.vector(F.car.topo.dist)
F.car.geo.vec = as.vector(F.car.geo.dist.3)

F.car.df = data.frame(F.car.topo.vec,F.car.geo.vec)

F.car.topo = ggplot(F.car.df, aes(y = F.car.topo.vec, x = F.car.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Fraxinus caroliniana R2 = 0.02", x = "Geographical Distance (km)", y = "elevation Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
F.car.topo

ggsave(F.car.topo, file = "./Results/test.results/Plots/F.car.elevation.pdf", height = 5, width = 5)

C.tex.topo.vec = as.vector(C.tex.topo.dist)
C.tex.geo.vec = as.vector(C.tex.geo.dist.3)

C.tex.df = data.frame(C.tex.topo.vec,C.tex.geo.vec)

C.tex.topo = ggplot(C.tex.df, aes(y = C.tex.topo.vec, x = C.tex.geo.vec)) + 
  geom_point(shape = 21, fill = "grey", color = "black") + 
  geom_smooth(color = "red")+
  geom_smooth(method = "lm", color = "blue")+
  labs(title = "Carya texana R2 = 0.016", x = "Geographical Distance (km)", y = "elevation Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
C.tex.topo

ggsave(C.tex.topo, file = "./Results/test.results/Plots/C.tex.elevation.pdf", height = 5, width = 5)
