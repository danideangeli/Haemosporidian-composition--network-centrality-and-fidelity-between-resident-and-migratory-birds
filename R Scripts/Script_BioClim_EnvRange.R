#### BAIXANDO DADOS DO BIOCLIM PARA AS COORDENADAS DA DANIELA ####
#### DeAngeli et al - Tese de Doutorado - Capitulo 3
#### Objetivos: Baixar dados climaticos a serem usados para determinar a amplitude ambiental dos parasitos

# Obs: Amostras 7077 e 7088 estão com coordenadas identicas mas referenciadas em locais diferentes. Pode ter sido erro na hora de entrar as coordenadas

#### PACOTES ####

library(reshape2)
#Spatial Analyzes - BioClim
library(raster)
library(sp)
library(dismo)
#PCA Analyzes
library(factoextra)
library(corrplot)
#Cluster Analyzes
library(cluster)
library(vegan)
library(ape)
#Hill analyzes
library("hillR")
#Plot
library(ggfortify)

#### CORES ####

## GGPLOT2 colors
# Recovering the four colors used for GGPLOT2 to plot the map 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(6)


#### DADOS ####

Table <- read.csv(file = "Sites.csv", sep = ";", na.strings = "", dec = ",") # Loading the data frame
Table$Latitude <- as.numeric(Table$Latitude) # Converting Lat to numeric
Table$Latitude[Table$Latitude == 0] <- NA # Converting zero to NA
Table$Longitude <- as.numeric(Table$Longitude) # Converting Long to numeric
Table$Longitude[Table$Longitude == 0] <- NA # converting zero to NA
Table <- subset(Table, !is.na(Table$Longitude) | !is.na(Table$Latitude2))

Table <- data.frame(Table)
Table <- filter(Table, Host_Environment == "Wild")
Table <- filter(Table, Latitude != 0)
Table <- filter(Table, parasiteGenus != "Leucocytozoon")

# Creating a ID variable to lat_long combinations 
Table$Loc <- paste(Table$Latitude,Table$Longitude, sep = "_") # Criando uma coluna com ID para os sites
Loc.ID <- cbind.data.frame(Loc = unique(Table$Loc), ID = 1:length(unique(Table$Loc)))
Table <- merge.data.frame(Table, Loc.ID, by = "Loc")
rm(Loc.ID)

# Creating an occurrence matrix
Occ.Mat <- dcast(data = Table, formula = Lineage ~ ID, fun.aggregate = length)
rownames(Occ.Mat) <- Occ.Mat[,1]
Occ.Mat <- Occ.Mat[,-1]
sort(rowSums(Occ.Mat))
colnames(Occ.Mat)

# Crating a df with selected columns to analyzes (Long, Lat, Continent and Country)
coords <- cbind.data.frame(x = Table$Longitude, y = Table$Latitude, Continent = Table$continent, Country = Table$country, ID = Table$ID)
coords <- coords[!duplicated(coords[,c(1,2,5)]),]
coords <- coords[complete.cases(coords),]
coords <- subset(coords, abs(coords$x) <= 180 & abs(coords$y) <= 180)

#### IMPORTING DATA ####
#dir.create(path = "Inputs/BioClim") # creating a directory for saving BioClim data
bioclim.data <- getData(name = "worldclim",
                        var = c("bio"),
                        res = 10) # Loading BioClim Data from WorldClim

coords.range <- c(min(coords$x),max(coords$x),min(coords$y),max(coords$y)) # Creating a vector with coords range
points <- SpatialPoints(coords[,c(1,2)], proj4string = bioclim.data@crs) # Creating a spatial object from coords
values <- extract(bioclim.data,points) # Extracting the BioClim values to our Sites coords

# Creating a Dataframe with all variables
df <- cbind.data.frame(coordinates(points),coords[,-c(1,2)],values)
df <- df[complete.cases(df),] # removing sites with NA values for BioClim
rownames(df) <- df$ID

df.st <- cbind.data.frame(df[,c(1:5)],scale(df[,-c(1:5)]))
colnames(df.st)[c(1,2)] <- c("long","lat")
rm(values, coords)

#### PLOTANDO OS PONTOS NO MAPA
png("Sites.png", width = 10, height = 6, units = "in", res = 300)
plot(points, add = T)
dev.off()

#### PCA of localities ####

PCA.Locs <- prcomp(x = df.st[,-c(1:5)])

png("Figures/PCA_Locs1.png", width = 10, height = 10, units = "in",res = 300)
autoplot(PCA.Locs, data = df.st, colour = "Continent")
dev.off()

png("Figures/PCA_Locs.png", width = 10, height = 10, units = "in",res = 300)
fviz_pca_biplot(PCA.Locs,
                col.var = "gray30", # Color by contributions to the PC
                col.ind = df.st$Continent, palette = cols,
                repel = T, addEllipses = F, geom.ind = c("point"),legend = "bottom",
                title = "", pointsize = 2.5,alpha = 0.75,
                invisible=c("quali"), ellipse.level=0.95, pointshape = 19
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
          legend.text = element_text(size = 15))
dev.off()

#### Cluster of Localities ####

Dist.Locs <- vegdist(df.st[,-c(1:5)], method = "euclidian")
Clust.Locs <- as.phylo(hclust(Dist.Locs))

Occ.Mat1 <- Occ.Mat[,Clust.Locs$tip.label]
Occ.Mat1 <- Occ.Mat1[rowSums(Occ.Mat1) != 0,]

HillLineages <- hill_phylo(comm = Occ.Mat1, tree = Clust.Locs)
EnvRange <- as.data.frame(HillLineages)
write.csv(EnvRange, file = "Results/EnvRange.csv", row.names = T)

#### Occurrence network ####

library(igraph)
G1 <- graph_from_incidence_matrix(Occ.Mat1, directed = F)
Vertex.color <- ifelse(V(G1)$name %in% rownames(Occ.Mat1),"skyblue", "orange")
plot(G1, vertex.label = "", vertex.size = 5, vertex.color = Vertex.color)
