library(tidyverse)
library(readxl)
library(data.table)
library(vegan)
library(gdata)
library(effects)
library(ggplot2)
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

##### Data manipulation #####

Table1$biomes <- as.factor(Table1$biomes)

dados = Table1 %>%
  group_by(biomes, Host_Status, add = TRUE) %>%
  dplyr::summarise(n_lineages = n(),
                   par_richness = n_distinct(Lineage_Name))%>%
  arrange(biomes)

dados1 <- filter(dados, par_richness > 2)

dados2 <- with(dados1, names(table(biomes)[table(biomes) > 1]))
dados2 <- dados1[dados1$biomes%in% dados2, ]
dados2$biomes <- factor(dados2$biomes)

Table1$BioSta <- paste(Table1$biomes, " ", Table1$Host_Status)

Table1$Host_Status=as.factor(paste(Table1$Host_Status))
Table1$Host_Status<- relevel(Table1$Host_Status, ref="R")

levels(Table1$Host_Status)

Occ.Mat <- dcast(data = Table1, formula =  BioSta ~ Lineage_Name, fun.aggregate = length)
rownames(Occ.Mat) <- Occ.Mat[,1]
Occ.Mat = Occ.Mat %>%
  arrange(BioSta)

species_category <- Table1[, c(13, 7)]
species_category <- distinct(species_category)
species_category = species_category %>%
  arrange(BioSta)
data.table::setcolorder(Occ.Mat, c("BioSta"))
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
plot(mod.coleta, main = NULL, sub = NULL)

##### Composition and dissimilarity analyses resident vs non-resident species #####

Table3 <- Table1
Table3$Host_Status <- gsub("PM", "M", Table3$Host_Status)

Table3$BioSta <- paste(Table3$biomes, " ", Table3$Host_Status)

Table3$Host_Status=as.factor(paste(Table3$Host_Status))
Table3$Host_Status<- relevel(Table3$Host_Status, ref="R")

levels(Table3$Host_Status)

Occ.Mat1 <- dcast(data = Table3, formula =  BioSta ~ Lineage_Name, fun.aggregate = length)
rownames(Occ.Mat1) <- Occ.Mat1[,1]
Occ.Mat1 = Occ.Mat1 %>%
  arrange(BioSta)

species_category1 <- Table3[, c(13, 7)]
species_category1 <- distinct(species_category1)
species_category1 = species_category1 %>%
  arrange(BioSta)
data.table::setcolorder(Occ.Mat1, c("BioSta"))
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
plot(mod.coleta1, main = NULL, sub = NULL)

##### MAP #####

newmap <- rworldmap::getMap(resolution = "high")
pdf(file="Fig 1.pdf", width = 8, height = 6)
plot(newmap, asp = 1, col = "white", xlim = c(-60, -45), ylim = c(-40, 10))
points(Table1$Longitude, Table1$Latitude, col= Table1$biomes, cex = 1.8, pch = 18)
legend(x = "bottomleft", title = "Biomes", pch = 18, cex = 1.0, pt.cex = 1.6, bty = "n",
       legend = c("Amazonia", "Andean Forest", "Caatinga", "Cerrado",
                                    "Grassland", "Atlantic Rain Forest", "Pantanal"),
       col = 1:7)
dev.off()
       #col = c("Black", "Red", "green3", "Blue", "cyan", 
ggsave(filename = "Fig1.pdf", scale = 2, dpi = 600)                                                #"magenta", "yellow"))

print(unique(Table1$biomes))

 
