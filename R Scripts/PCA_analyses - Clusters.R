library(tidyverse)
library(readxl)
library(data.table)
library(vegan)
library(gdata)
library(effects)
library(ggplot2)
library(pca3d)
library(reshape2)
library(factoextra)
library(corrplot)
library(cluster)
library(plyr)
library(sf)
library(sp)
library(raster)
library(dismo)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3/Tabelas")

##### Loading and preparing data #####

Table <- read_excel("Cap3.xlsx")
Table$Latitude <- as.numeric(Table$Latitude) # Converting Lat to numeric
Table$Latitude[Table$Latitude == 0] <- NA # Converting zero to NA
Table$Longitude <- as.numeric(Table$Longitude) # Converting Long to numeric
Table$Longitude[Table$Longitude == 0] <- NA # converting zero to NA
Table <- subset(Table, !is.na(Table$Longitude) | !is.na(Table$Latitude)) # Excluding possible NAs
Table <- data.frame(Table)
Table <- filter(Table, Latitude != "NA")
Table <- filter(Table, Longitude != "NA")

# Creating a ID variable to lat_long combinations 
Table$Loc <- paste(Table$Latitude,Table$Longitude, sep = "__") # Criando uma coluna com ID para os sites
Loc.ID <- cbind.data.frame(Loc = unique(Table$Loc), ID = 1:length(unique(Table$Loc)))
Table <- merge.data.frame(Table, Loc.ID, by = "Loc")
rm(Loc.ID)

Table1 <- with(Table, names(table(Loc)[table(Loc) > 9]))
Table1 <- Table[Table$Loc%in% Table1, ]
Table1$Loc <- factor(Table1$Loc)

##### Clustaring data by distance #####

#my.sf.point <- st_as_sf(x = Table, 
                        #coords = c("Longitude", "Latitude"),
                        #crs = "+proj=longlat +datum=WGS84")

#distance matrix in feet
#st_distance(my.sf.point)

#which poiint are within 50 miles (~80467.2 meters)
#l <- st_is_within_distance(my.sf.point, dist = 5000000)

#Table$within_1000km <- rowSums(as.matrix(l))-1

#m <- as.matrix(l)
#colnames(m) <- c(1:nrow(Table))
#rownames(m) <- c(1:nrow(Table))
#Table$points_within_1000km <- apply( m, 1, function(u) paste( names(which(u)), collapse="," ) )
#Table$clusterid <- dplyr::group_indices(Table, Table$points_within_1000km) 

#length(unique(Table$clusterid))

#Table$points_within_1000km <- NULL
#Table$within_1000km <- NULL

#pts <- st_as_sf(Table, coords = c("Latitude", "Longitude"),remove=FALSE)
#pts_agg <- aggregate(pts,
                     #pts$geometry,
                     #FUN = mean, 
                     #join = function(Latitude, Longitude) st_is_within_distance(Latitude, Longitude, dist = 50))


#Table <- plyr::arrange(Table, clusterid)
#data.table::setcolorder(Table, c("clusterid"))

##### Clustering data by climate regions #####

sites <- Table1[(9:10)]
sites <- distinct(sites)

coords <- cbind.data.frame(x = Table1$Longitude, y = Table1$Latitude, ID = Table1$ID)
coords <- coords[!duplicated(coords[,c(1,2,3)]),]
coords <- coords[complete.cases(coords),]
coords <- subset(coords, abs(coords$x) <= 180 & abs(coords$y) <= 180)

bioclim.data <- getData(name = "worldclim",
                        var = c("bio"),
                        res = 10) # Loading BioClim Data from WorldClim

coords.range <- c(min(coords$x),max(coords$x),min(coords$y),max(coords$y)) # Creating a vector with coords range
points <- SpatialPoints(coords[,c(1,2)], proj4string = bioclim.data@crs) # Creating a spatial object from coords
values <- extract(bioclim.data,points) # Extracting the BioClim values to our Sites coords

sites <- cbind.data.frame(sites[(1:2)], values)
sites <- sites[-c(49),]
sites.scale <- scale(sites[3:21])
sites.scale <- as.data.frame(sites.scale)
sites.scale <- cbind.data.frame(sites[,1:2], sites.scale)

dis.sites <- dist(sites.scale) 
dist.sites2 <- as.dist(1-cor(t(sites.scale)))

hier.clust <- hclust(dis.sites,method="complete")
plot(hier.clust)
rect.hclust(hier.clust, h = 10, border = "red")

pamclu <- cluster::pam((sites.scale),k=10)
plot(silhouette(pamclu),main=NULL)
cluster.data <- as.data.frame(pamclu$silinfo$widths)
setDT(cluster.data, keep.rownames = TRUE)[]

setDT(sites, keep.rownames = TRUE)[]

clusters <- inner_join(sites, cluster.data, by = "rn")
clusters <- clusters[,2:23]
clusters <- clusters %>%
  arrange(cluster)
setcolorder(clusters, c("cluster"))

clusters$Loc <- paste(clusters$Latitude,clusters$Longitude, sep = "__") 
clustersID <- clusters[,c(1,23)]

Table2 <- inner_join(Table1, clustersID, by = "Loc")
setcolorder(Table2, c("cluster"))

##### Clustered data manipulation #####

dados = Table2 %>%
  group_by(cluster, Host_Status, add = TRUE) %>%
  dplyr::summarise(n_lineages = n(),
            par_richness = n_distinct(Lineage_Name))%>%
  arrange(cluster)

dados1 <- filter(dados, par_richness > 2)

dados2 <- with(dados1, names(table(cluster)[table(cluster) > 1]))
dados2 <- dados1[dados1$cluster%in% dados2, ]
dados2$cluster <- factor(dados2$cluster)

Table2$CluSta <- paste(Table2$cluster, " ", Table2$Host_Status)

Table2$Host_Status=as.factor(paste(Table2$Host_Status))
Table2$Host_Status<- relevel(Table2$Host_Status, ref="R")

levels(Table2$Host_Status)

Occ.Mat <- dcast(data = Table2, formula =  CluSta ~ Lineage_Name, fun.aggregate = length)
rownames(Occ.Mat) <- Occ.Mat[,1]
Occ.Mat = Occ.Mat %>%
  arrange(CluSta)

species_category <- Table2[, c(14, 8)]
species_category <- distinct(species_category)
species_category = species_category %>%
  arrange(CluSta)
data.table::setcolorder(Occ.Mat, c("CluSta"))
rownames(Occ.Mat) <- Occ.Mat[,1]
Occ.Mat <- Occ.Mat[,-c(1,2)]
sort(rowSums(Occ.Mat))
colnames(Occ.Mat)

attach(species_category)

##### Composition and dissimilarity analyses #####

dis <- vegdist(Occ.Mat,method="bray")
mod.coleta <- betadisper(dis, Host_Status ,type="centroid")
anova(mod.coleta)

permutest(mod.coleta, pairwise = T, permutations = 999)


mod.coleta.HSD <- TukeyHSD(mod.coleta)
plot(mod.coleta.HSD)
plot(mod.coleta)

#comp <- anosim(mod.coleta, Host_Status, permutations = 999, distance = "bray", strata = NULL,
       #parallel = getOption("mc.cores"))
#summary(comp)
#plot(comp)

##### Composition and dissimilarity analyses resident vs non-resident species #####

Table3 <- Table2
Table3$Host_Status <- gsub("PM", "M", Table3$Host_Status)

Table3$CluSta <- paste(Table3$cluster, " ", Table3$Host_Status)

Table3$Host_Status=as.factor(paste(Table3$Host_Status))
Table3$Host_Status<- relevel(Table3$Host_Status, ref="R")

levels(Table3$Host_Status)

Occ.Mat1 <- dcast(data = Table3, formula =  CluSta ~ Lineage_Name, fun.aggregate = length)
rownames(Occ.Mat1) <- Occ.Mat1[,1]
Occ.Mat1 = Occ.Mat1 %>%
  arrange(CluSta)

species_category1 <- Table3[, c(14, 8)]
species_category1 <- distinct(species_category1)
species_category1 = species_category1 %>%
  arrange(CluSta)
data.table::setcolorder(Occ.Mat1, c("CluSta"))
rownames(Occ.Mat1) <- Occ.Mat1[,1]
Occ.Mat1 <- Occ.Mat1[,-c(1,2)]
sort(rowSums(Occ.Mat1))
colnames(Occ.Mat1)

attach(species_category1)

dis1 <- vegdist(Occ.Mat1,method="bray")
mod.coleta1 <- betadisper(dis1, Host_Status ,type="centroid")
anova(mod.coleta1)

permutest(mod.coleta1, pairwise = T, permutations = 999)

mod.coleta.HSD1 <- TukeyHSD(mod.coleta1)
plot(mod.coleta.HSD1)
plot(mod.coleta1)

##### MAP #####
library(rworldmap)

newmap <- getMap(resolution = "high")
plot(newmap, asp = 1, col = "white", xlim = c(-60, -45), ylim = c(-40, 10))

points(clusters$Longitude, clusters$Latitude, col= clusters$cluster, cex = 2,
       pch = c(15,16,17,18,19,20,8,9,10,11)[(pch=c(clusters$cluster))])

#Defining specific colors and types of points for each cluster

cluster1 <- filter(clusters, cluster == 1)
points(cluster1$Longitude, cluster1$Latitude, col= "red", cex = 2,pch = 15)

cluster2 <- filter(clusters, cluster == 2)
points(cluster2$Longitude, cluster2$Latitude, col= "blue", cex = 2,pch = 16)

cluster3 <- filter(clusters, cluster == 3)
points(cluster3$Longitude, cluster3$Latitude, col= "green", cex = 2,pch = 17)

cluster4 <- filter(clusters, cluster == 4)
points(cluster4$Longitude, cluster4$Latitude, col= "purple", cex = 2,pch = 18)

cluster5 <- filter(clusters, cluster == 5)
points(cluster5$Longitude, cluster5$Latitude, col= "orange", cex = 2,pch = 19)

cluster6 <- filter(clusters, cluster == 6)
points(cluster6$Longitude, cluster6$Latitude, col= "yellow", cex = 2,pch = 20)

cluster7 <- filter(clusters, cluster == 7)
points(cluster7$Longitude, cluster7$Latitude, col= "grey", cex = 2,pch = 21)

cluster8 <- filter(clusters, cluster == 8)
points(cluster8$Longitude, cluster8$Latitude, col= "black", cex = 2,pch = 22)

cluster9 <- filter(clusters, cluster == 9)
points(cluster9$Longitude, cluster9$Latitude, col= "brown", cex = 2,pch = 23)

cluster10 <- filter(clusters, cluster == 10)
points(cluster10$Longitude, cluster10$Latitude, col= "black", cex = 2,pch = 24)



