R Scripts: Haemosporidian taxonomic composition, network centrality and partner fidelity between resident and migratory avian hosts 
Daniela de Angeli Dutra¹*, Alan Fecchio², Érika Martins Braga³, Robert Poulin¹

#Network Analyses

library(tidyverse)
library(readxl)
library(data.table)
library(bipartite)
library(picante)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(treeman)
library(ape)
library(devtools)
library(phytools)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3/Tabelas")

dados <- read_excel("cap3.xlsx") # loading data

##### Data Manipulation #####

dados1 = dados %>%
  group_by(Lineage_Name, species, biomes, add = TRUE) %>%
  summarise(n_interactions = n()) %>%
  arrange(Lineage_Name)

dados2 <- dados %>%
  group_by(Lineage_Name, biomes, add = TRUE) %>%
  summarise(n_local_partners = n(),
            n_local_species = n_distinct(species)) %>%
  arrange(Lineage_Name)

dados3<- dados %>%
  group_by(biomes, add = TRUE) %>%
  summarise(n_putative_species = n_distinct(species),
            n_putative_infections = n()) %>%
  arrange(biomes)

dados2 <- filter(dados2, n_local_partners > 9)

dados2$biomes <- NULL

dados4 <- inner_join(dados1, dados2, by = "Lineage_Name") 
dados4 <- as.data.table(dados4)

dados4a <- filter(dados4, biomes == "Amazonia")
dados4b <- filter(dados4, biomes == "Andinean_Forest")
dados4c <- filter(dados4, biomes == "Caatinga")
dados4d <- filter(dados4, biomes == "Cerrado")
dados4f <- filter(dados4, biomes == "Mata_Atlantica")

##### Creating occurrence matrices #####

Occ.Mata <- dcast.data.table(dados4a, Lineage_Name ~ species, fun = sum, value.var = "n_interactions")
Occ.Mata <- column_to_rownames(Occ.Mata, var = "Lineage_Name")

nullmatrixA <- vaznull(9999, Occ.Mata)

Occ.Matb <- dcast.data.table(dados4b, Lineage_Name ~ species, fun = sum, value.var = "n_interactions")
Occ.Matb <- column_to_rownames(Occ.Matb, var = "Lineage_Name")

nullmatrixB <- vaznull(9999, Occ.Matb)

Occ.Matc <- dcast.data.table(dados4c, Lineage_Name ~ species, fun = sum, value.var = "n_interactions")
Occ.Matc <- column_to_rownames(Occ.Matc, var = "Lineage_Name")

nullmatrixC <- vaznull(9999, Occ.Matc)

Occ.Matd <- dcast.data.table(dados4d, Lineage_Name ~ species, fun = sum, value.var = "n_interactions")
Occ.Matd <- column_to_rownames(Occ.Matd, var = "Lineage_Name")

nullmatrixD <- vaznull(9999, Occ.Matd)

Occ.Matf <- dcast.data.table(dados4f, Lineage_Name ~ species, fun = sum, value.var = "n_interactions")
Occ.Matf <- column_to_rownames(Occ.Matf, var = "Lineage_Name")

nullmatrixF <- vaznull(9999, Occ.Matf)

##### Calculating normalized degree and other measures ####

HostPropertiesa <- specieslevel(Occ.Mata, level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
HostPropertiesa <- rownames_to_column(HostPropertiesa, var = "species")
HostPropertiesa$biomes <- "Amazonia"

HostPropertiesb <- specieslevel(Occ.Matb, level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
HostPropertiesb <- rownames_to_column(HostPropertiesb, var = "species")
HostPropertiesb$biomes <- "Andinean_Forest"

HostPropertiesc <- specieslevel(Occ.Matc, level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
HostPropertiesc <- rownames_to_column(HostPropertiesc, var = "species")
HostPropertiesc$biomes <- "Caatinga"

HostPropertiesd <- specieslevel(Occ.Matd, level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
HostPropertiesd <- rownames_to_column(HostPropertiesd, var = "species")
HostPropertiesd$biomes <- "Cerrado"

HostPropertiesf <- specieslevel(Occ.Matf, level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
HostPropertiesf <- rownames_to_column(HostPropertiesf, var = "species")
HostPropertiesf$biomes <- "Mata_Atlantica"

## Combaning all biomes into one data ##

HostProperties <- rbind.data.frame(HostPropertiesa,HostPropertiesb, HostPropertiesc, HostPropertiesd,
                                   HostPropertiesf) %>%
  arrange(species)

## Adding host migratory status information ##

dados5 <- cbind.data.frame(dados$species, dados$Host_Status)
dados5 <- distinct(dados5)

dados6 <- inner_join(HostProperties, dados5, by = c("species" = "dados$species"))
names(dados6)[names(dados6) == "dados$Host_Status"] <- "Host_Status"

##### Loop Vectors #####

ZNullPropertiesaND <- matrix(NA, 81, 9999)
ZNullPropertiesbND <- matrix(NA, 89, 9999)
ZNullPropertiescND <- matrix(NA, 34, 9999)
ZNullPropertiesdND <- matrix(NA, 73, 9999)
ZNullPropertiesfND <- matrix(NA, 68, 9999)

ZNullPropertiesaWC <- matrix(NA, 81, 9999)
ZNullPropertiesbWC <- matrix(NA, 89, 9999)
ZNullPropertiescWC <- matrix(NA, 34, 9999)
ZNullPropertiesdWC <- matrix(NA, 73, 9999)
ZNullPropertiesfWC <- matrix(NA, 68, 9999)

ZNullPropertiesaWB <- matrix(NA, 81, 9999)
ZNullPropertiesbWB <- matrix(NA, 89, 9999)
ZNullPropertiescWB <- matrix(NA, 34, 9999)
ZNullPropertiesdWB <- matrix(NA, 73, 9999)
ZNullPropertiesfWB <- matrix(NA, 68, 9999)

##### LOOP #####

for(i in 1:9999){
  print(c("iteration",i))
  
  ## Calculating normalised degree and other indixes for ranzomizations##
  
  NullPropertiesa <- specieslevel(nullmatrixA[[i]], level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
  
  NullPropertiesb <- specieslevel(nullmatrixB[[i]], level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
  
  NullPropertiesc <- specieslevel(nullmatrixC[[i]], level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
  
  NullPropertiesd <- specieslevel(nullmatrixD[[i]], level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
  
  NullPropertiesf <- specieslevel(nullmatrixF[[i]], level = "higher", index=c("normalised degree", "betweenness","closeness"), PDI.normalise=F)
  
  ZNullPropertiesaND [,i] <- NullPropertiesa$normalised.degree
  ZNullPropertiesbND [,i] <- NullPropertiesb$normalised.degree
  ZNullPropertiescND [,i] <- NullPropertiesc$normalised.degree
  ZNullPropertiesdND [,i] <- NullPropertiesd$normalised.degree
  ZNullPropertiesfND [,i] <- NullPropertiesf$normalised.degree
  
  ZNullPropertiesaWC [,i] <- NullPropertiesa$weighted.closeness
  ZNullPropertiesbWC [,i] <- NullPropertiesb$weighted.closeness
  ZNullPropertiescWC [,i] <- NullPropertiesc$weighted.closeness
  ZNullPropertiesdWC [,i] <- NullPropertiesd$weighted.closeness
  ZNullPropertiesfWC [,i] <- NullPropertiesf$weighted.closeness
  
  ZNullPropertiesaWB [,i] <- NullPropertiesa$weighted.betweenness
  ZNullPropertiesbWB [,i] <- NullPropertiesb$weighted.betweenness
  ZNullPropertiescWB [,i] <- NullPropertiesc$weighted.betweenness
  ZNullPropertiesdWB [,i] <- NullPropertiesd$weighted.betweenness
  ZNullPropertiesfWB [,i] <- NullPropertiesf$weighted.betweenness
  
}

NullPropertiesa <- cbind.data.frame(select(HostPropertiesa,species, normalised.degree, weighted.betweenness, weighted.closeness),
                                    normalised.degree = ZNullPropertiesaND, weighted.betweenness = ZNullPropertiesaWB,weighted.closeness = ZNullPropertiesaWC)
NullPropertiesb <- cbind.data.frame(select(HostPropertiesb,species, normalised.degree, weighted.betweenness, weighted.closeness),
                                    normalised.degree = ZNullPropertiesbND, weighted.betweenness = ZNullPropertiesbWB,weighted.closeness = ZNullPropertiesbWC)
NullPropertiesc <- cbind.data.frame(select(HostPropertiesc,species, normalised.degree, weighted.betweenness, weighted.closeness),
                                    normalised.degree = ZNullPropertiescND, weighted.betweenness = ZNullPropertiescWB,weighted.closeness = ZNullPropertiescWC)
NullPropertiesd <- cbind.data.frame(select(HostPropertiesd,species, normalised.degree, weighted.betweenness, weighted.closeness),
                                    normalised.degree = ZNullPropertiesdND, weighted.betweenness = ZNullPropertiesdWB,weighted.closeness = ZNullPropertiesdWC)
NullPropertiesf <- cbind.data.frame(select(HostPropertiesf,species, normalised.degree, weighted.betweenness, weighted.closeness),
                                    normalised.degree = ZNullPropertiesfND,weighted.betweenness = ZNullPropertiesfWB,weighted.closeness = ZNullPropertiesfWC)


NullPropertiesa$biomes <- "Amazonia"
NullPropertiesb$biomes <- "Andinean_Forest"
NullPropertiesc$biomes <- "Caatinga"
NullPropertiesd$biomes <- "Cerrado"
NullPropertiesf$biomes <- "Mata_Atlantica"

## Combaning all biomes into one data ##

NullProperties <- rbind.data.frame(NullPropertiesa,NullPropertiesb, NullPropertiesc, NullPropertiesd,
                                   NullPropertiesf) %>%
  arrange(species)

## Adding host migratory status information ##

Nulldados6 <- inner_join(NullProperties, dados5, by = c("species" = "dados$species"))
names(Nulldados6)[names(Nulldados6) == "dados$Host_Status"] <- "Host_Status"


##### Adding phylogenetic information #####

load("random_trees2.RData")
tree <- random_trees2[[1]]
tree <- as.phylo(tree)
fullbirds <- as.data.frame(tree$tip.label)
mybirds <- as.data.frame(unique(dados6$species))
names(mybirds)[names(mybirds) == "unique(dados6$species)"] <- "species"

length(which(tree$tip.label%in%as.character(mybirds[,1])))
todrop<-tree$tip.label[which(tree$tip.label%in%as.character(mybirds[,1])==FALSE)]
mybirds_tree<-drop.tip(tree,todrop)

keep.tip(tree,as.character(mybirds[,1]))

Nulldados7 <-
  Nulldados6 %>%
  filter(species %in% mybirds_tree$tip.label)

length(unique(Nulldados7$species))
length(mybirds_tree$tip.label)

##### Calculating Z-Score #####

Nulldados7$MeanND <- rowMeans(Nulldados7[5:10003])
Nulldados7$MeanWB <- rowMeans(Nulldados7[10004:20002])
Nulldados7$MeanWC <- rowMeans(Nulldados7[20003:30001])
Nulldados7$sdND <- matrixStats::rowSds(as.matrix(Nulldados7[5:10003]))
Nulldados7$sdWB <- matrixStats::rowSds(as.matrix(Nulldados7[10004:20002]))
Nulldados7$sdWC <- matrixStats::rowSds(as.matrix(Nulldados7[20003:30001]))

Nulldados7$ZND <- (Nulldados7$normalised.degree-Nulldados7$MeanND)/Nulldados7$sdND
Nulldados7$ZWB <- (Nulldados7$weighted.betweenness-Nulldados7$MeanWB)/Nulldados7$sdWB
Nulldados7$ZWC <- (Nulldados7$weighted.closeness-Nulldados7$MeanWC)/Nulldados7$sdWC

save.image("NullTeste2.RData")

##### Bayesian Models #####

dados7 <- select(Nulldados7, "species", "Host_Status", "biomes","ZND", "ZWB", "ZWC")

library(brms)

dados7$Host_Status <- as.factor(paste(dados7$Host_Status)) 
dados7$Host_Status <- relevel(dados7$Host_Status, ref = "R") 
levels(dados7$Host_Status)

hist(dados7$ZND, breaks = 100) # testing data distribution 
hist(dados7$ZWB, breaks = 100) # testing data distribution 
hist(dados7$ZWC, breaks = 100) # testing data distribution 

## Creating phylogenetic matrix ##

inv.phylo <- MCMCglmm::inverseA(mybirds_tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
class(A)

variaveis <- brms::bf(ZND~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = skew_normal())

prior <- get_prior(variaveis, data = dados7) #getting priors
prior

model <- brm(ZND~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7,
             family = skew_normal(), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(student_t(3, -0.3, 2.5), "Intercept"),                                             
               prior(student_t(3, 0, 2.5), "sd"),
               prior(student_t(3, 0, 2.5), "sigma")
             ))
summary(model)

a <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot1ND <- plot(conditional_effects(model), points = FALSE, theme = a)

plot2ND <- plot1ND$Host_Status + labs(x = "Migratory host category", y = "Normalised degree") #caption = "Plasmodium")
plot2ND

variaveis <- brms::bf(ZWB~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = skew_normal())

prior <- get_prior(variaveis, data = dados7) #getting priors
prior

modelWB <- brm(ZWB~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7,
               family = skew_normal(), chains = 4,
               iter = 4000,
               data2 = list(A = A),
               prior = c(
                 prior(student_t(3, -0.1, 2.5), "Intercept"),                                             
                 prior(student_t(3, 0, 2.5), "sd"),
                 prior(student_t(3, 0, 2.5), "sigma")
               ))
summary(modelWB)

plot1WB <- plot(conditional_effects(modelWB), points = FALSE, theme = a)

plot2WB <- plot1WB$Host_Status + labs(x = "Migratory host category", y = "Weoghted betweenness") #caption = "Plasmodium")
plot2WB

variaveis <- brms::bf(ZWC~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = gaussian())

prior <- get_prior(variaveis, data = dados7) #getting priors
prior

modelWC <- brm(ZWC~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7,
               family = gaussian(), chains = 4,
               iter = 4000,
               data2 = list(A = A),
               prior = c(
                 prior(student_t(3, 0.1, 2.5), "Intercept"),                                             
                 prior(student_t(3, 0, 2.5), "sd"),
                 prior(student_t(3, 0, 2.5), "sigma")
               ))
summary(modelWC)

plot1WC <- plot(conditional_effects(modelWC), points = FALSE, theme = a)

plot2WC <- plot1WC$Host_Status + labs(x = "Migratory host category", y = "Weighted closeness") #caption = "Plasmodium")
plot2WC

##### Analysing using resident vs non-resident data only #####

dados8 <- dados7
dados8$Host_Status <- gsub("PM", "M", dados8$Host_Status)
print(unique(dados8$Host_Status))

dados8$Host_Status <- as.factor(dados8$Host_Status)
dados8$Host_Status <- relevel(dados8$Host_Status, ref = "R")

variaveis1 <- brms::bf(ZND~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family =  skew_normal()) #check

prior1 <- get_prior(variaveis1, data = dados8) #getting priors
prior1

model1 <- brm(ZND~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8,
              family = skew_normal(), chains = 4,
              iter = 4000,
              data2 = list(A = A),
              prior = c(
                prior(student_t(3, -0.3, 2.5), "Intercept"),                                             
                prior(student_t(3, 0, 2.5), "sd"),
                prior(student_t(3, 0, 2.5), "sigma")
              ))
summary(model1)

plot3ND <- plot(conditional_effects(model1), points = FALSE, theme = a)

plot4ND <- plot3ND$Host_Status + labs(x = "Migratory host category", y = "Normalised degree") #caption = "Plasmodium")
plot4ND

variaveis1 <- brms::bf(ZWB~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family =  skew_normal()) #check

prior1 <- get_prior(variaveis1, data = dados8) #getting priors
prior1

model1WB <- brm(ZWB~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8,
                family = skew_normal(), chains = 4,
                iter = 4000,
                data2 = list(A = A),
                prior = c(
                  prior(student_t(3, -0.1, 2.5), "Intercept"),                                             
                  prior(student_t(3, 0, 2.5), "sd"),
                  prior(student_t(3, 0, 2.5), "sigma")
                ))
summary(model1WB)

plot3WB <- plot(conditional_effects(model1WB), points = FALSE, theme = a)

plot4WB <- plot3WB$Host_Status + labs(x = "Migratory host category", y = "weighted betweenness") #caption = "Plasmodium")
plot4WB

variaveis1 <- brms::bf(ZWC~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family =  gaussian()) #check

prior1 <- get_prior(variaveis1, data = dados8) #getting priors
prior1

model1WC <- brm(ZWC~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8,
                family = gaussian(), chains = 4,
                iter = 4000,
                data2 = list(A = A),
                prior = c(
                  prior(student_t(3, 0.1, 2.5), "Intercept"),                                             
                  prior(student_t(3, 0, 2.5), "sd"),
                  prior(student_t(3, 0, 2.5), "sigma")
                ))
summary(model1WC)

plot3WC <- plot(conditional_effects(model1WC), points = FALSE, theme = a)

plot4WC <- plot3WC$Host_Status + labs(x = "Migratory host category", y = "weighted closeness") #caption = "Plasmodium")
plot4WC

save.image("FinalZscoredModels.RData")

#Haemosporidian Taxonomic Composition Analyses

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
