#### Trying ConR for EOO estimation ####
# EOO: Extent of Occurrence
# https://gdauby.github.io/ConR/index.html
# https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.11230

#install.packages("devtools")
#devtools::install_github("gdauby/ConR")
print(library(ConR))

gbif = read.csv("gbif.final.csv", row.names = 1)

# subset and rearrange the columns
occ.range = gbif[,c(2,3,1)]
EOO.norm = EOO.computing(occ.range, file.name = "EOO")
write.csv(EOO.norm, file = "EOO.norm.csv")
EOO.alpha = EOO.computing(occ.range, method.range = "alpha.hull", file.name = "EOO.alpha")
write.csv(EOO.alpha, file = "EOO.alpha.csv")