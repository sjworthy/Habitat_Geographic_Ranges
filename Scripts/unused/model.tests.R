# Example Analyses with select species

library(tidyverse)
library(FD)
library(geosphere)
library(ecodist)
library(gghalves)

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

#### Bootstrapping Microclimate distance ~ Geographic Distance ####

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

# Function to perform bootstrapping and calculate standard effect size
bootstrap_mrm <- function(microclim_dist, geo_dist, n_bootstrap = 999) {
  
  # Convert geo_dist from class 'dist' to a matrix if needed
  geo_dist_matrix <- as.matrix(geo_dist)
  
  # Get the number of rows (occurrences) from the dist object
  n <- nrow(geo_dist_matrix)  # nrow will work after conversion to matrix
  
  # Initialize vector to store intercept and slope values
  slopes <- numeric(n_bootstrap)
  intercepts <- numeric(n_bootstrap)
  
  # Bootstrap loop
  for (i in 1:n_bootstrap) {
    
    # Sample with replacement from the geo distance matrix
    sampled_indices <- sample(1:n, n, replace = TRUE)
    
    # Create the bootstrapped geo distance matrix
    boot_geo_dist_matrix <- geo_dist_matrix[sampled_indices, sampled_indices]
    
    # Convert the bootstrapped matrix back to a 'dist' object
    boot_geo_dist <- as.dist(boot_geo_dist_matrix)

    # Run the MRM model with the bootstrapped geo distance matrix
    model <- MRM(microclim_dist ~ boot_geo_dist)
    
    # Extract the intercept and slope values
    intercept_value <- model$coef[1,1]
    slope_value <- model$coef[2,1]
    
    intercepts[i] <- intercept_value
    slopes[i] <- slope_value
    
  }
  
  # combine intercept and slope vectors into a dataframe
  null_output = as.data.frame(intercepts)
  null_output$slopes = slopes
  
  return(null_output)
}

### Q.cocc 
# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(Q.cocc.microclim.dist, Q.cocc.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.1927249, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Null Distribution of Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/Q.cocc.null.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.000007711331, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Null Distribution of Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/Q.cocc.null.slope.pdf", width = 5, height = 5)

clim.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
clim.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
clim.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
clim.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.1927249 - clim.intercept.null.means) / clim.intercept.nulls.sds
ses.slope <- (0.000007711331 - clim.slope.null.means) / clim.slope.nulls.sds

rank.intercept = rank(c(0.1927249,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.142

rank.slope= rank(c(0.000007711331,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 0.86

### A.leu

# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(A.leu.microclim.dist, A.leu.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.1341791297, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "A. leu Null Distribution of Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/A.leu.null.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.0001561188, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "A. leu Null Distribution of Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/A.leu.null.slope.pdf", width = 5, height = 5)

clim.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
clim.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
clim.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
clim.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.1341791297 - clim.intercept.null.means) / clim.intercept.nulls.sds
ses.slope <- (0.0001561188 - clim.slope.null.means) / clim.slope.nulls.sds

rank.intercept = rank(c(0.1341791297,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.0001561188,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

## F.car
# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(F.car.microclim.dist, F.car.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.0926204409 , color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F.car Null Distribution of Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/F.car.null.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.0001909017, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F. car Null Distribution of Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/F.car.null.slope.pdf", width = 5, height = 5)

clim.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
clim.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
clim.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
clim.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.0926204409 - clim.intercept.null.means) / clim.intercept.nulls.sds
ses.slope <- (0.0001909017 - clim.slope.null.means) / clim.slope.nulls.sds

rank.intercept = rank(c(0.0926204409,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.0001909017,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

## C.tex

# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(C.tex.microclim.dist, C.tex.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.0502335274 , color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "C.tex Null Distribution of Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/C.tex.null.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.0004084846, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F. car Null Distribution of Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/C.tex.null.slope.pdf", width = 5, height = 5)

clim.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
clim.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
clim.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
clim.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.0502335274 - clim.intercept.null.means) / clim.intercept.nulls.sds
ses.slope <- (0.0004084846 - clim.slope.null.means) / clim.slope.nulls.sds

rank.intercept = rank(c(0.0502335274,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.0004084846,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

#### Bootstrapping Topo distance ~ Geographic Distance ####

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

### Q.cocc 
# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(Q.cocc.topo.dist, Q.cocc.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.1605788312, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Null Distribution of Topo Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/Q.cocc.topo.null.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = -0.0000418118, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Null Distribution of Topo Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/Q.cocc.topo.null.slope.pdf", width = 5, height = 5)

topo.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
topo.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
topo.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
topo.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.1605788312 - topo.intercept.null.means) / topo.intercept.nulls.sds
ses.slope <- (-0.0000418118 - topo.slope.null.means) / topo.slope.nulls.sds

rank.intercept = rank(c(0.1605788312,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.142

rank.slope= rank(c(-0.0000418118,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 0.86

### A.leu

# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(A.leu.topo.dist, A.leu.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.1359138, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "A. leu Null Distribution of Topo Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/A.leu.null.topo.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.00001993009, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "A. leu Null Distribution of Topo Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/A.leu.null.topo.slope.pdf", width = 5, height = 5)

topo.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
topo.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
topo.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
topo.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.1359138 - topo.intercept.null.means) / topo.intercept.nulls.sds
ses.slope <- (0.00001993009 - topo.slope.null.means) / topo.slope.nulls.sds

rank.intercept = rank(c(0.1359138,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.00001993009,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

## F.car
# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(F.car.topo.dist, F.car.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.06911964 , color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F.car Null Distribution of Topo Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/F.car.null.topo.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.00001189631, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F. car Null Distribution of Topo Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/F.car.null.topo.slope.pdf", width = 5, height = 5)

topo.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
topo.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
topo.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
topo.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.06911964 - topo.intercept.null.means) / topo.intercept.nulls.sds
ses.slope <- (0.00001189631 - topo.slope.null.means) / topo.slope.nulls.sds

rank.intercept = rank(c(0.06911964,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.00001189631,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

## C.tex

# Run the bootstrap procedure
bootstrapped_effect_sizes <- bootstrap_mrm(C.tex.topo.dist, C.tex.geo.dist.3, n_bootstrap = 999)

null.intercept.plot = ggplot(bootstrapped_effect_sizes, aes(x = intercepts)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.09998326 , color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "C.tex Null Distribution of Topo Intercept with Observed Value",
    x = "Intercept Null Values",
    y = "Density") +
  theme_minimal()
null.intercept.plot

ggsave(null.intercept.plot, file = "./Results/test.results/Plots/C.tex.null.topo.intercept.pdf", width = 5, height = 5)

null.slope.plot = ggplot(bootstrapped_effect_sizes, aes(x = slopes)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = 0.00003588177, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "F. car Null Distribution of Topo Slope with Observed Value",
    x = "Slope Null Values",
    y = "Density") +
  theme_minimal()
null.slope.plot

ggsave(null.slope.plot, file = "./Results/test.results/Plots/C.tex.null.topo.slope.pdf", width = 5, height = 5)

topo.intercept.null.means <- mean(bootstrapped_effect_sizes$intercepts)
topo.intercept.nulls.sds <- sd(bootstrapped_effect_sizes$intercepts)
topo.slope.null.means <- mean(bootstrapped_effect_sizes$slopes)
topo.slope.nulls.sds <- sd(bootstrapped_effect_sizes$slopes)

ses.intercet <- (0.09998326 - topo.intercept.null.means) / topo.intercept.nulls.sds
ses.slope <- (0.00003588177 - topo.slope.null.means) / topo.slope.nulls.sds

rank.intercept = rank(c(0.09998326,bootstrapped_effect_sizes$intercepts))[1]
p.val.intercept = rank.intercept/1000 # 0.001

rank.slope= rank(c(0.00003588177,bootstrapped_effect_sizes$slopes))[1]
p.val.slope = rank.slope/1000 # 1

#### Density Plots ####

# combine species' data

Q.cocc.all = left_join(Q.cocc,Q.cocc.elev)
A.leu.all = left_join(A.leu,A.leu.elev)
F.car.all = left_join(F.car,F.car.elev)
C.tex.all = left_join(C.tex,C.text.elev)

all.dat = rbind(Q.cocc.all,A.leu.all,F.car.all,C.tex.all)

## elevation plot
elevation.plot = ggplot(all.dat, aes(x = species, y = elevation), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Elevation (m)",
       x = " ", fill = " ") +
  theme_classic()
elevation.plot

ggsave(elevation.plot, file = "./Results/test.results/Plots/elevation.pdf", height = 5, width = 5)

## high temp plot
high.temp.plot = ggplot(all.dat, aes(x = species, y = high_temp_C), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "High Temp (C)",
       x = " ", fill = " ") +
  theme_classic()
high.temp.plot

ggsave(high.temp.plot, file = "./Results/test.results/Plots/high_temp.pdf", height = 5, width = 5)

## low temp plot
low.temp.plot = ggplot(all.dat, aes(x = species, y = low_temp_C), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Low Temp (C)",
       x = " ", fill = " ") +
  theme_classic()
low.temp.plot

ggsave(low.temp.plot, file = "./Results/test.results/Plots/low_temp.pdf", height = 5, width = 5)

## microclimate distance plot

Q.cocc.df = as.data.frame(as.vector(Q.cocc.microclim.dist))
colnames(Q.cocc.df)[1] = "Microclimate.Distance"
Q.cocc.df$species = "Quercus coccinea"

A.leu.df = as.data.frame(as.vector(A.leu.microclim.dist))
colnames(A.leu.df)[1] = "Microclimate.Distance"
A.leu.df$species = "Acer leucoderme"

F.car.df = as.data.frame(as.vector(F.car.microclim.dist))
colnames(F.car.df)[1] = "Microclimate.Distance"
F.car.df$species = "Fraxinus caroliniana"

C.tex.df = as.data.frame(as.vector(C.tex.microclim.dist))
colnames(C.tex.df)[1] = "Microclimate.Distance"
C.tex.df$species = "Carya texana"

microclim.df = rbind(Q.cocc.df,A.leu.df,F.car.df,C.tex.df)

microclim.dist.plot = ggplot(microclim.df, aes(x = species, y = Microclimate.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Microclimate Distance",
       x = " ", fill = " ") +
  theme_classic()
microclim.dist.plot

ggsave(microclim.dist.plot, file = "./Results/test.results/Plots/microclim.dist.density.pdf", height = 5, width = 5)

## topographic distance plot

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

Q.cocc.df = as.data.frame(as.vector(Q.cocc.topo.dist))
colnames(Q.cocc.df)[1] = "Topo.Distance"
Q.cocc.df$species = "Quercus coccinea"

A.leu.df = as.data.frame(as.vector(A.leu.topo.dist))
colnames(A.leu.df)[1] = "Topo.Distance"
A.leu.df$species = "Acer leucoderme"

F.car.df = as.data.frame(as.vector(F.car.topo.dist))
colnames(F.car.df)[1] = "Topo.Distance"
F.car.df$species = "Fraxinus caroliniana"

C.tex.df = as.data.frame(as.vector(C.tex.topo.dist))
colnames(C.tex.df)[1] = "Topo.Distance"
C.tex.df$species = "Carya texana"

topo.df = rbind(Q.cocc.df,A.leu.df,F.car.df,C.tex.df)

topo.dist.plot = ggplot(topo.df, aes(x = species, y = Topo.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Topographic Distance",
       x = " ", fill = " ") +
  theme_classic()
topo.dist.plot

ggsave(topo.dist.plot, file = "./Results/test.results/Plots/topo.dist.density.pdf", height = 5, width = 5)

## microclimate dist plot by geographic dist

Q.cocc.df = as.data.frame(as.vector(Q.cocc.microclim.dist))
colnames(Q.cocc.df)[1] = "Microclimate.Distance"
Q.cocc.df$species = "Quercus coccinea"
Q.cocc.df$Geographic.Distance = as.vector(Q.cocc.geo.dist.3)

# Create quantile bins (5 bins) and convert to factor
Q.cocc.df$Geographic.Distance.Quantile <- cut(
  Q.cocc.df$Geographic.Distance,
  breaks = quantile(Q.cocc.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(Q.cocc.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0.0000  222.6753  371.6079  515.2171  710.5724 1424.8978 

Q.cocc.plot = ggplot(Q.cocc.df, aes(x = as.factor(Geographic.Distance.Quantile), y = Microclimate.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Microclimate Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
Q.cocc.plot

ggsave(Q.cocc.plot, file = "./Results/test.results/Plots/Q.cocc.bins.pdf", height = 5, width = 5)

A.leu.df = as.data.frame(as.vector(A.leu.microclim.dist))
colnames(A.leu.df)[1] = "Microclimate.Distance"
A.leu.df$species = "Acer leucoderme"
A.leu.df$Geographic.Distance = as.vector(A.leu.geo.dist.3)

# Create quantile bins (5 bins) and convert to factor
A.leu.df$Geographic.Distance.Quantile <- cut(
  A.leu.df$Geographic.Distance,
  breaks = quantile(A.leu.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(A.leu.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0.0000  191.6128  365.9663  562.1573  800.9971 1542.0604  

A.leu.plot = ggplot(A.leu.df, aes(x = as.factor(Geographic.Distance.Quantile), y = Microclimate.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Microclimate Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
A.leu.plot

ggsave(A.leu.plot, file = "./Results/test.results/Plots/A.leu.bins.pdf", height = 5, width = 5)

F.car.df = as.data.frame(as.vector(F.car.microclim.dist))
colnames(F.car.df)[1] = "Microclimate.Distance"
F.car.df$species = "Fraxinus caroliniana"
F.car.df$Geographic.Distance = as.vector(F.car.geo.dist.3)

# Create quantile bins (5 bins) and convert to factor
F.car.df$Geographic.Distance.Quantile <- cut(
  F.car.df$Geographic.Distance,
  breaks = quantile(F.car.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(F.car.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0.0000  228.7815  434.5618  688.4255  942.7673 1948.1028

F.car.plot = ggplot(F.car.df, aes(x = as.factor(Geographic.Distance.Quantile), y = Microclimate.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Microclimate Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
F.car.plot

ggsave(F.car.plot, file = "./Results/test.results/Plots/F.car.bins.pdf", height = 5, width = 5)


C.tex.df = as.data.frame(as.vector(C.tex.microclim.dist))
colnames(C.tex.df)[1] = "Microclimate.Distance"
C.tex.df$species = "Carya texana"
C.tex.df$Geographic.Distance = as.vector(C.tex.geo.dist.3)

# Create quantile bins (5 bins) and convert to factor
C.tex.df$Geographic.Distance.Quantile <- cut(
  C.tex.df$Geographic.Distance,
  breaks = quantile(C.tex.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

breaks = quantile(C.tex.df$Geographic.Distance, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
#        0%       20%       40%       60%       80%      100% 
# 0.0000  213.9236  354.7697  513.3723  745.1326 1374.4654

C.tex.plot = ggplot(C.tex.df, aes(x = as.factor(Geographic.Distance.Quantile), y = Microclimate.Distance), fill = "darkgray") +
  geom_half_point(side = "l", size = 0.1, color = "darkgray",
                  position = position_nudge(x=-.05), alpha = 0.3) +
  geom_half_boxplot(fill = NA, position = position_nudge(x=-.05)) +
  geom_half_violin(aes(fill = "darkgray"), side = "r", fill = "darkgray",
                   scale = "width") +
  labs(y = "Microclimate Distance",
       x = "Geographic Distance Quantile", fill = " ") +
  theme_classic()
C.tex.plot

ggsave(C.tex.plot, file = "./Results/test.results/Plots/C.text.bins.pdf", height = 5, width = 5)


