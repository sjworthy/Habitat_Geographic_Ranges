library(tidyverse)
library(FD)
library(geosphere)
library(ecodist) # get weird error sometimes for MRM, just reinstall package

# Two options for doing analysis: scale or not the distance matrices
# Q1.MRM results scaled within species and did not scale distance matrices

# first step for both workflows is to scale the microenvironmental data for all species at the same time

all.clim = read.csv("./Results/climate.data.output.1.csv")

# get columns we need, scale microenvironment columns

all.clim.2 = all.clim %>%
  dplyr::select(species,decimalLongitude,decimalLatitude,high_temp_C,low_temp_C,moisture_mm) %>%
  mutate(high.temp.scaled = scale(high_temp_C),
         low.temp.scaled = scale(low_temp_C),
         moisture.scaled = scale(moisture_mm))

# Test for three species: Acer_leucoderme, Aesculus_glabra, Fraxinus_caroliniana

A_leu = all.clim.2 %>%
  filter(species == "Acer leucoderme")
A_gla = all.clim.2 %>%
  filter(species == "Aesculus glabra")
F_car = all.clim.2 %>% 
  filter(species == "Fraxinus caroliniana")

#### Method 1: scale, center response, no scaling of distance matrix ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
F.car.microclim.dat = F_car %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

# distance decay plot

A.gla.microclim.vec = as.vector(A.gla.microclim.dat)
A.gla.geo.vec = as.vector(A.gla.geo.dist.2)

A.gla.df = data.frame(A.gla.microclim.vec,A.gla.geo.vec)

A.gla.high.temp = ggplot(A.gla.df, aes(y = high.temp.scaled, x = A.gla.geo.vec)) + 
  #geom_point(size = 3, alpha = 0.5) + 
  geom_smooth()+
  geom_smooth(method = "lm")+
  labs(x = "Geographical Distance", y = "High Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.gla.high.temp

A.gla.cold.temp = ggplot(A.gla.df, aes(y = low.temp.scaled, x = A.gla.geo.vec)) + 
  #geom_point(size = 3, alpha = 0.5) + 
  geom_smooth()+
  geom_smooth(method = "lm")+
  labs(x = "Geographical Distance", y = "Cold Temp Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.gla.cold.temp

A.gla.moist.temp = ggplot(A.gla.df, aes(y = moisture.scaled, x = A.gla.geo.vec)) + 
 #geom_point(size = 3, alpha = 0.5) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  labs(x = "Geographical Distance", y = "Moisture Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
A.gla.moist.temp


# Perform MRM
A.gla.MRM.mod = MRM(A.gla.microclim.dist ~ A.gla.geo.dist.2)
# intercept = 5.792871e-02, slope = 2.551610e-07, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.2)
# intercept = 1.341791e-01, slope = 1.561188e-07, p = 0.001, R2 = 0.2339119
F.car.MRM.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.2)
# intercept = 9.262044e-02, slope = 1.909017e-07, p = 0.001, R2 = 0.5023438

#### Method 2: scale, center response, scale, center of distance matrix #####
# Function to standardize variables in distance format
decostandDist<-function(dist,method="standardize",...){
  
  require(vegan)
  
  distCopy<-dist
  #  lab<-labels(dist)
  
  #  Mat<-as.matrix(dist)
  #  nMat<-nrow(Mat)
  #  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  distStd<-decostand(as.vector(dist),"standardize",...)
  distCopy[1:length(distCopy)]<-distStd
  #  rownames(resuMat)<-lab
  #  colnames(resuMat)<-lab
  
  #  resu<-as.dist(resuMat)
  return(distCopy)
  
}

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
F.car.microclim.dat = F_car %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

A.gla.microclim.dist.Std = decostandDist(A.gla.microclim.dist)
A.leu.microclim.dist.Std = decostandDist(A.leu.microclim.dist)
F.car.microclim.dist.Std = decostandDist(F.car.microclim.dist)

A.gla.geo.dist.Std = decostandDist(A.gla.geo.dist.2)
A.leu.geo.dist.Std = decostandDist(A.leu.geo.dist.2)
F.car.geo.dist.Std = decostandDist(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Std ~ A.gla.geo.dist.Std)
# intercept = 1.450922e-14, p = 0.001, slope = 6.432393e-01, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Std ~ A.leu.geo.dist.Std)
# intercept = -1.889351e-16, p = 0.072, slope = 4.836444e-01, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Std ~ F.car.geo.dist.Std)
# intercept = -1.890622e-15, p = 0.001, slope = 7.087622e-01, p = 0.001, R2 = 0.5023438

#### Method 3 and 4 ####

# when scaling microenvironment, only scaling, no centering.

all.clim.var = all.clim %>%
  dplyr::select(species,decimalLongitude,decimalLatitude,high_temp_C,low_temp_C,moisture_mm) %>%
  mutate(across(c(4:6), ~ . / sd(.)))

# Test for three species: Acer_leucoderme, Aesculus_glabra, Fraxinus_caroliniana

A_leu = all.clim.var %>%
  filter(species == "Acer leucoderme")
A_gla = all.clim.var %>%
  filter(species == "Aesculus glabra")
F_car = all.clim.var %>% 
  filter(species == "Fraxinus caroliniana")

#### Method 3: scale response, no scaling of distance matrix ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

# Perform MRM
A.gla.MRM.mod = MRM(A.gla.microclim.dist ~ A.gla.geo.dist.2)
# intercept = 5.792871e-02, slope = 2.551610e-07, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.2)
# intercept = 1.341791e-01, slope = 1.561188e-07, p = 0.001, R2 = 0.2339119
F.car.MRM.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.2)
# intercept = 9.262044e-02, slope = 1.909017e-07, p = 0.001, R2 = 0.5023438

#### Method 4 scale response, scale, center of distance matrix ####
decostandDist<-function(dist,method="standardize",...){
  
  require(vegan)
  
  distCopy<-dist
  #  lab<-labels(dist)
  
  #  Mat<-as.matrix(dist)
  #  nMat<-nrow(Mat)
  #  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  distStd<-decostand(as.vector(dist),"standardize",...)
  distCopy[1:length(distCopy)]<-distStd
  #  rownames(resuMat)<-lab
  #  colnames(resuMat)<-lab
  
  #  resu<-as.dist(resuMat)
  return(distCopy)
  
}

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

A.gla.microclim.dist.Std = decostandDist(A.gla.microclim.dist)
A.leu.microclim.dist.Std = decostandDist(A.leu.microclim.dist)
F.car.microclim.dist.Std = decostandDist(F.car.microclim.dist)

A.gla.geo.dist.Std = decostandDist(A.gla.geo.dist.2)
A.leu.geo.dist.Std = decostandDist(A.leu.geo.dist.2)
F.car.geo.dist.Std = decostandDist(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Std ~ A.gla.geo.dist.Std)
# intercept = 1.450922e-14, p = 0.001, slope = 6.432393e-01, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Std ~ A.leu.geo.dist.Std)
# intercept = -3.554178e-16, p = 0.013, slope = 4.836444e-01, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Std ~ F.car.geo.dist.Std)
# intercept = -1.926621e-15, p = 0.001, slope = 7.087622e-01, p = 0.001, R2 = 0.5023438


#### Method 5: no scaling response, no scaling distance ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

# Perform MRM
A.gla.MRM.mod = MRM(A.gla.microclim.dist ~ A.gla.geo.dist.2)
# intercept = 5.792871e-02, p = 1.000 slope = 2.551610e-07, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod = MRM(A.leu.microclim.dist ~ A.leu.geo.dist.2)
# intercept = 1.341791e-01, slope = 1.561188e-07, p = 0.001, R2 = 0.2339119
F.car.MRM.mod = MRM(F.car.microclim.dist ~ F.car.geo.dist.2)
# intercept = 9.262044e-02, slope = 1.909017e-07, p = 0.001, R2 = 0.5023438

#### Method 6: no scaling response, scale, center distance ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

A.gla.microclim.dist.Std = decostandDist(A.gla.microclim.dist)
A.leu.microclim.dist.Std = decostandDist(A.leu.microclim.dist)
F.car.microclim.dist.Std = decostandDist(F.car.microclim.dist)

A.gla.geo.dist.Std = decostandDist(A.gla.geo.dist.2)
A.leu.geo.dist.Std = decostandDist(A.leu.geo.dist.2)
F.car.geo.dist.Std = decostandDist(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Std ~ A.gla.geo.dist.Std)
# intercept = 1.450922e-14, p = 0.001, slope = 6.432393e-01, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Std ~ A.leu.geo.dist.Std)
# intercept = -2.332170e-16, p = 0.023, slope = 4.836444e-01, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Std ~ F.car.geo.dist.Std)
# intercept = -1.991440e-15, p = 0.001, slope = 7.087622e-01, p = 0.001, R2 = 0.5023438



#### Method 7: scale, center response, scale distance matrix ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)
F.car.microclim.dat = F_car %>%
  dplyr::select(high.temp.scaled,low.temp.scaled,moisture.scaled)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

decostandDist.Var <-function(dist){
  
  require(vegan)
  
  distCopy<-dist
  #  lab<-labels(dist)
  
  #  Mat<-as.matrix(dist)
  #  nMat<-nrow(Mat)
  #  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  dist.vec<-as.vector(dist)
  distVar<-dist.vec/sd(dist.vec)
  distCopy[1:length(distCopy)]<-distVar
  #  rownames(resuMat)<-lab
  #  colnames(resuMat)<-lab
  
  #  resu<-as.dist(resuMat)
  return(distCopy)
  
}

A.gla.microclim.dist.Var = decostandDist.Var(A.gla.microclim.dist)
A.leu.microclim.dist.Var = decostandDist.Var(A.leu.microclim.dist)
F.car.microclim.dist.Var = decostandDist.Var(F.car.microclim.dist)

A.gla.geo.dist.Var = decostandDist.Var(A.gla.geo.dist.2)
A.leu.geo.dist.Var = decostandDist.Var(A.leu.geo.dist.2)
F.car.geo.dist.Var = decostandDist.Var(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Var ~ A.gla.geo.dist.Var)
# intercept = 0.4215578, p = 1.000, slope = 0.6432393, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Var ~ A.leu.geo.dist.Var)
# intercept = 1.1797645, p = 1.000, slope = 0.4836444, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Var ~ F.car.geo.dist.Var)
# intercept = 0.9259137, p = 1.000, slope = 0.7087622, p = 0.001, R2 = 0.5023438


#### Method 8: scale response, scale distance matrix ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

A.gla.microclim.dist.Var = decostandDist.Var(A.gla.microclim.dist)
A.leu.microclim.dist.Var = decostandDist.Var(A.leu.microclim.dist)
F.car.microclim.dist.Var = decostandDist.Var(F.car.microclim.dist)

A.gla.geo.dist.Var = decostandDist.Var(A.gla.geo.dist.2)
A.leu.geo.dist.Var = decostandDist.Var(A.leu.geo.dist.2)
F.car.geo.dist.Var = decostandDist.Var(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Var ~ A.gla.geo.dist.Var)
# intercept = 0.4215578, p = 1, slope = 0.6432393, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Var ~ A.leu.geo.dist.Var)
# intercept = 1.1797645, p = 1.00, slope = 0.4836444, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Var ~ F.car.geo.dist.Var)
# intercept = 0.9259137, p = 1.00, slope = 0.7087622, p = 0.001, R2 = 0.5023438

#### Method 9: no scaling response, scale distance matrix ####

# creating spatial matrix
A.gla.spat.dat = as.matrix(A_gla[,c(2,3)])
A.leu.spat.dat = as.matrix(A_leu[,c(2,3)])
F.car.spat.dat = as.matrix(F_car[,c(2,3)])

# creating microclim data, scaled
A.gla.microclim.dat = A_gla %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
A.leu.microclim.dat = A_leu %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)
F.car.microclim.dat = F_car %>%
  dplyr::select(high_temp_C,low_temp_C,moisture_mm)

# calculate gower distance for scaled microclimate data
A.gla.microclim.dist = gowdis(as.matrix(A.gla.microclim.dat))
A.leu.microclim.dist = gowdis(as.matrix(A.leu.microclim.dat))
F.car.microclim.dist = gowdis(as.matrix(F.car.microclim.dat))

# calculate Haversine distance for spatial data
A.gla.geo.dist = distm(A.gla.spat.dat, fun = distHaversine)
A.gla.geo.dist.2 = as.dist(A.gla.geo.dist) # convert to dist object

A.leu.geo.dist = distm(A.leu.spat.dat, fun = distHaversine)
A.leu.geo.dist.2 = as.dist(A.leu.geo.dist) # convert to dist object

F.car.geo.dist = distm(F.car.spat.dat, fun = distHaversine)
F.car.geo.dist.2 = as.dist(F.car.geo.dist) # convert to dist object

decostandDist.Var <-function(dist){
  
  require(vegan)
  
  distCopy<-dist
  #  lab<-labels(dist)
  
  #  Mat<-as.matrix(dist)
  #  nMat<-nrow(Mat)
  #  resuMat<-matrix(decostand(as.vector(Mat),method=method,...),nMat,nMat)
  dist.vec<-as.vector(dist)
  distVar<-dist.vec/sd(dist.vec)
  distCopy[1:length(distCopy)]<-distVar
  #  rownames(resuMat)<-lab
  #  colnames(resuMat)<-lab
  
  #  resu<-as.dist(resuMat)
  return(distCopy)
  
}

A.gla.microclim.dist.Var = decostandDist.Var(A.gla.microclim.dist)
A.leu.microclim.dist.Var = decostandDist.Var(A.leu.microclim.dist)
F.car.microclim.dist.Var = decostandDist.Var(F.car.microclim.dist)

A.gla.geo.dist.Var = decostandDist.Var(A.gla.geo.dist.2)
A.leu.geo.dist.Var = decostandDist.Var(A.leu.geo.dist.2)
F.car.geo.dist.Var = decostandDist.Var(F.car.geo.dist.2)

# Perform MRM
A.gla.MRM.mod.std = MRM(A.gla.microclim.dist.Var ~ A.gla.geo.dist.Var)
# intercept = 0.4215578, p = 1.000, slope = 0.6432393, p = 0.001, R2 = 0.4137568
A.leu.MRM.mod.std = MRM(A.leu.microclim.dist.Var ~ A.leu.geo.dist.Var)
# intercept = 1.1797645, p = 1.000, slope = 0.4836444, p = 0.001, R2 = 0.2339119
F.car.MRM.mod.std = MRM(F.car.microclim.dist.Var ~ F.car.geo.dist.Var)
# intercept = 0.9259137, p = 1.000, slope = 0.7087622, p = 0.001, R2 = 0.5023438