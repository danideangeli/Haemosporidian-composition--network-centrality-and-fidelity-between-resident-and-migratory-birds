library(tidyverse)
library(readxl)
library(vegan)
library(gdata)
library(effects)
library(ggplot2)
library(pca3d)
library(reshape2)
library(factoextra)
library(corrplot)
library(cluster)
library(vegan)
library(plyr)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3")

Table <- read_excel("Cap3.xlsx")
Table$Latitude <- as.numeric(Table$Latitude) # Converting Lat to numeric
Table$Latitude[Table$Latitude == 0] <- NA # Converting zero to NA
Table$Longitude <- as.numeric(Table$Longitude) # Converting Long to numeric
Table$Longitude[Table$Longitude == 0] <- NA # converting zero to NA
Table <- subset(Table, !is.na(Table$Longitude) | !is.na(Table$Latitude)) # Excluding possible NAs
Table <- data.frame(Table)

# Creating a ID variable to lat_long combinations 
Table$Loc <- paste(Table$Latitude,Table$Longitude, sep = "_") # Criando uma coluna com ID para os sites
Loc.ID <- cbind.data.frame(Loc = unique(Table$Loc), ID = 1:length(unique(Table$Loc)))
Table <- merge.data.frame(Table, Loc.ID, by = "Loc")
rm(Loc.ID)

Table <- plyr::arrange(Table, ID)
data.table::setcolorder(Table, c("ID"))

dados = Table %>%
  group_by(Loc, Host_Status, add = TRUE) %>%
  dplyr::summarise(n_lineages = n(),
            par_richness = n_distinct(Lineage_Name))%>%
  arrange(Loc)

dados1 <- filter(dados, par_richness > 2)

dados2 <- with(dados1, names(table(Loc)[table(Loc) > 2]))
dados2 <- dados1[dados1$Loc%in% dados2, ]
dados2$Loc <- factor(dados2$Loc)

Table1 <- with(Table, names(table(Loc)[table(Loc) > 9]))
Table1 <- Table[Table$Loc%in% Table1, ]
Table1$Loc <- factor(Table1$Loc)

Occ.Mat <- dcast(data = Table, formula = species ~ Lineage_Name, fun.aggregate = length)
rownames(Occ.Mat) <- Occ.Mat[,1]

species_category <- Table[, c(7, 8)]
species_category <- distinct(species_category)
Occ.Mat <- inner_join(Occ.Mat, species_category, by = "species")
data.table::setcolorder(Occ.Mat, c("species", "Host_Status"))
rownames(Occ.Mat) <- Occ.Mat[,1]
Occ.Mat <- Occ.Mat[,-c(1,2)]
sort(rowSums(Occ.Mat))
colnames(Occ.Mat)

attach(species_category)

dis <- vegdist(Occ.Mat,method="bray")
mod.coleta <- betadisper(dis, Host_Status ,type="centroid")
anova(mod.coleta)

permutest(mod.coleta, pairwise = T, permutations = 999)

mod.coleta.HSD <- TukeyHSD(mod.coleta)
plot(mod.coleta.HSD)
plot(mod.coleta)

