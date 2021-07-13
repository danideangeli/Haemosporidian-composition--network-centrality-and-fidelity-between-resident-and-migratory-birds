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

dados4 <- inner_join(dados1, dados2, by = "Lineage_Name") %>%
        mutate(ND_Host_Parasite = n_interactions/n_local_partners)
dados4 <- as.data.table(dados4)

dados4a <- filter(dados4, biomes == "Amazonia")
dados4b <- filter(dados4, biomes == "Andinean_Forest")
dados4c <- filter(dados4, biomes == "Caatinga")
dados4d <- filter(dados4, biomes == "Cerrado")
#dados4e <- filter(dados4, biomes == "Grassland")
dados4f <- filter(dados4, biomes == "Mata_Atlantica")
#dados4g <- filter(dados4, biomes == "Pantanal")
#dados4h <- filter(dados4, biomes == "Patagonian_Forest")

## Creating occurrence matrixes ##

Occ.Mata <- dcast.data.table(dados4a, Lineage_Name ~ species, fun = sum)
Occ.Mata <- column_to_rownames(Occ.Mata, var = "Lineage_Name")

Occ.Matb <- dcast.data.table(dados4b, Lineage_Name ~ species, fun = sum)
Occ.Matb <- column_to_rownames(Occ.Matb, var = "Lineage_Name")

Occ.Matc <- dcast.data.table(dados4c, Lineage_Name ~ species, fun = sum)
Occ.Matc <- column_to_rownames(Occ.Matc, var = "Lineage_Name")

Occ.Matd <- dcast.data.table(dados4d, Lineage_Name ~ species, fun = sum)
Occ.Matd <- column_to_rownames(Occ.Matd, var = "Lineage_Name")

Occ.Matf <- dcast.data.table(dados4f, Lineage_Name ~ species, fun = sum)
Occ.Matf <- column_to_rownames(Occ.Matf, var = "Lineage_Name")

## Calculating normalised degree and other indixes##

HostPropertiesa <- specieslevel(Occ.Mata, level = "higher", index=c("normalised degree", "PDI","effective partners", 
                                                                "species specificity"), PDI.normalise=F)
HostPropertiesa <- rownames_to_column(HostPropertiesa, var = "species")
HostPropertiesa$biomes <- "Amazonia"

HostPropertiesb <- specieslevel(Occ.Matb, level = "higher", index=c("normalised degree", "PDI","effective partners", 
                                                                  "species specificity"), PDI.normalise=F)
HostPropertiesb <- rownames_to_column(HostPropertiesb, var = "species")
HostPropertiesb$biomes <- "Andinean_Forest"

HostPropertiesc <- specieslevel(Occ.Matc, level = "higher", index=c("normalised degree", "PDI","effective partners", 
                                                                  "species specificity"), PDI.normalise=F)
HostPropertiesc <- rownames_to_column(HostPropertiesc, var = "species")
HostPropertiesc$biomes <- "Caatinga"

HostPropertiesd <- specieslevel(Occ.Matd, level = "higher", index=c("normalised degree", "PDI","effective partners", 
                                                                  "species specificity"), PDI.normalise=F)
HostPropertiesd <- rownames_to_column(HostPropertiesd, var = "species")
HostPropertiesd$biomes <- "Cerrado"

HostPropertiesf <- specieslevel(Occ.Matf, level = "higher", index=c("normalised degree", "PDI","effective partners", 
                                                                  "species specificity"), PDI.normalise=F)
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

write.csv2(dados6, "Fidelitytable2.csv")

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

dados7 <-
  dados6 %>%
  filter(species %in% mybirds_tree$tip.label)

length(unique(dados7$species))
length(mybirds_tree$tip.label)

Netdata <- left_join(dados7, dados4, by = "species")

##### Bayesian Models #####

library(brms)

dados7$Host_Status <- as.factor(paste(dados7$Host_Status)) 
dados7$Host_Status <- relevel(dados7$Host_Status, ref = "R") 
levels(dados7$Host_Status)

hist(dados7$normalised.degree) # testing data distribution 
ggpubr::ggdensity(dados7$normalised.degree)# testing data distribution

## Creating phylogenetic matrix ##

inv.phylo <- MCMCglmm::inverseA(mybirds_tree, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
class(A)
#A <- ape::vcv.phylo(mybirds_tree)

variaveis <- brms::bf(normalised.degree~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = zero_one_inflated_beta)

prior <- get_prior(variaveis, data = dados7) #getting priors
prior

model <- brm(normalised.degree~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7,
             family = zero_one_inflated_beta(), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(beta(1, 1), "coi"),
               prior(student_t(3, 0, 2.5), "Intercept"),                                             
               prior(gamma(0.01, 0.01),"phi"),
               prior(student_t(3, 0, 2.5), "sd"),
               prior(beta(1, 1),"zoi")
             ))
summary(model)

a <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot1 <- plot(conditional_effects(model), points = FALSE, theme = a)

plot2 <- plot1$Host_Status + labs(x = "Migratory host category", y = "Normalised degree") #caption = "Plasmodium")
plot2


##### Analysing using resident vs non-resident data only #####

dados8 <- dados7
dados8$Host_Status <- gsub("PM", "M", dados8$Host_Status)
print(unique(dados8$Host_Status))

dados8$Host_Status <- as.factor(dados8$Host_Status)
dados8$Host_Status <- relevel(dados8$Host_Status, ref = "R")

variaveis1 <- brms::bf(normalised.degree~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family =  zero_one_inflated_beta()) #check

prior1 <- get_prior(variaveis1, data = dados8) #getting priors
prior1

model1 <- brm(normalised.degree~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8,
             family = zero_one_inflated_beta(), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(beta(1, 1), "coi"),
               prior(student_t(3, 0, 2.5), "Intercept"),                                             
               prior(gamma(0.01, 0.01),"phi"),
               prior(student_t(3, 0, 2.5), "sd"),
               prior(beta(1, 1),"zoi")
             ))
summary(model1)

plot3 <- plot(conditional_effects(model1), points = FALSE, theme = a)

plot4 <- plot3$Host_Status + labs(x = "Migratory host category", y = "Normalised degree") #caption = "Plasmodium")
plot4

#save.image("Fidelity Data.RData") # Saving R Data

##### Computing closeness and betwenness #####

#Occ.Mat <- dcast.data.table(dados4, Lineage_Name ~ species, fun = sum)
#Occ.Mat <- column_to_rownames(Occ.Mat, var = "Lineage_Name")

#NullMatrix <- vaznull(10000, Occ.Mat)

#NullProp <- specieslevel(NullMatrix[[1]], level = "higher", index=c("betweenness", "closeness"), PDI.normalise=F)

#list_of_matrix<-list(matrix(1:1000,nrow=100),matrix(1:1000,nrow=100))
#for (i in 1:10){
  #print (matrix)
  #NullProp[,i] <- specieslevel(NullMatrix[[i]], level = "higher", index=c("betweenness", "closeness"), PDI.normalise=F)
  #do some other task
  #print(matrix*10)}

##### Analysing closeness and betweenness #####

CloseBetweena <- specieslevel(Occ.Mata, level = "higher", index=c("closeness", "betweenness"), PDI.normalise=F)
CloseBetweena <- rownames_to_column(CloseBetweena, var = "species")
CloseBetweena$biomes <- "Amazonia"

CloseBetweenb <- specieslevel(Occ.Matb, level = "higher", index=c("closeness", "betweenness"), PDI.normalise=F)
CloseBetweenb <- rownames_to_column(CloseBetweenb, var = "species")
CloseBetweenb$biomes <- "Andinean_Forest"

CloseBetweenc <- specieslevel(Occ.Matc, level = "higher", index=c("closeness", "betweenness"), PDI.normalise=F)
CloseBetweenc <- rownames_to_column(CloseBetweenc, var = "species")
CloseBetweenc$biomes <- "Caatinga"

CloseBetweend <- specieslevel(Occ.Matd, level = "higher", index=c("closeness", "betweenness"), PDI.normalise=F)
CloseBetweend <- rownames_to_column(CloseBetweend, var = "species")
CloseBetweend$biomes <- "Cerrado"

CloseBetweenf <- specieslevel(Occ.Matf, level = "higher", index=c("closeness", "betweenness"), PDI.normalise=F)
CloseBetweenf <- rownames_to_column(CloseBetweenf, var = "species")
CloseBetweenf$biomes <- "Mata_Atlantica"

## Combaning all biomes into one data ##

CloseBetween <- rbind.data.frame(CloseBetweena,CloseBetweenb, CloseBetweenc, CloseBetweend,
                                   CloseBetweenf) %>%
  arrange(species)

## Adding host migratory status information ##

dados5 <- cbind.data.frame(dados$species, dados$Host_Status)
dados5 <- distinct(dados5)

dados6B <- inner_join(CloseBetween, dados5, by = c("species" = "dados$species"))
names(dados6B)[names(dados6B) == "dados$Host_Status"] <- "Host_Status"

write.csv2(dados6, "Centrality2.csv")

dados7B <-
  dados6B %>%
  filter(species %in% mybirds_tree$tip.label)

length(unique(dados7B$species))
length(mybirds_tree$tip.label)

NetdataB <- left_join(dados7B, dados4, by = "species")

##### Bayesian Models #####

library(brms)

dados7B$Host_Status <- as.factor(paste(dados7B$Host_Status)) 
dados7B$Host_Status <- relevel(dados7B$Host_Status, ref = "R") 
levels(dados7B$Host_Status)

hist(dados7B$weighted.closeness) # testing data distribution 
ggpubr::ggdensity(dados7B$weighted.closeness)# testing data distribution

##### Closeness and betweenness Bayesian models #####

variaveis3 <- brms::bf(weighted.closeness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = zero_one_inflated_beta)

prior <- get_prior(variaveis3, data = dados7B) #getting priors
prior

modelB <- brm(weighted.closeness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7B,
             family = zero_one_inflated_beta(), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(beta(1, 1), "coi"),
               prior(student_t(3, 0, 2.5), "Intercept"),                                             
               prior(gamma(0.01, 0.01),"phi"),
               prior(student_t(3, 0, 2.5), "sd"),
               prior(beta(1, 1),"zoi")
             ))
summary(modelB)

plot1B <- plot(conditional_effects(model), points = FALSE, theme = a)

plot2B <- plot1$Host_Status + labs(x = "Migratory host category", y = "Weighted Closeness") #caption = "Plasmodium")
plot2B

plot2$data

##### Analysing using resident vs non-resident data only #####

dados8B <- dados7B
dados8B$Host_Status <- gsub("PM", "M", dados8B$Host_Status)
print(unique(dados8B$Host_Status))

dados8B$Host_Status <- as.factor(dados8B$Host_Status)
dados8B$Host_Status <- relevel(dados8B$Host_Status, ref = "R")

variaveis1 <- brms::bf(weighted.closeness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family =  zero_one_inflated_beta()) #check

prior1 <- get_prior(variaveis1, data = dados8B) #getting priors
prior1

model1B <- brm(weighted.closeness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8B,
              family = zero_one_inflated_beta(), chains = 4,
              iter = 4000,
              data2 = list(A = A),
              prior = c(
                prior(beta(1, 1), "coi"),
                prior(student_t(3, 0, 2.5), "Intercept"),                                             
                prior(gamma(0.01, 0.01),"phi"),
                prior(student_t(3, 0, 2.5), "sd"),
                prior(beta(1, 1),"zoi")
              ))
summary(model1B)

plot3B <- plot(conditional_effects(model1B), points = FALSE, theme = a)

Plot4B <- plot3B$Host_Status + labs(x = "Migratory host category", y = "Weighted Closeness") #caption = "Plasmodium")
Plot4B

##### Betweenness models ######

hist(dados7B$weighted.betweenness) # testing data distribution 
ggpubr::ggdensity(dados7B$weighted.betweenness)# testing data distribution

variaveis1 <- brms::bf(weighted.betweenness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), family = zero_one_inflated_beta)

prior <- get_prior(variaveis, data = dados7B) #getting priors
prior

modelC <- brm(weighted.betweenness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados7B,
             family = zero_one_inflated_beta(), chains = 4,
             iter = 4000,
             data2 = list(A = A),
             prior = c(
               prior(beta(1, 1), "coi"),
               prior(student_t(3, 0, 2.5), "Intercept"),                                             
               prior(gamma(0.01, 0.01),"phi"),
               prior(student_t(3, 0, 2.5), "sd"),
               prior(beta(1, 1),"zoi")
             ))
summary(modelC)

plot5 <- plot(conditional_effects(modelC), points = FALSE, theme = a)

plot6 <- plot5$Host_Status + labs(x = "Migratory host category", y = "Weighted Betweenness") #caption = "Plasmodium")
plot6


##### Resident versus non-resident analyses #####

model1C <- brm(weighted.betweenness~Host_Status + (1|biomes) + (1|gr(species, cov = A)), data = dados8B,
              family = zero_one_inflated_beta(), chains = 4,
              iter = 4000,
              data2 = list(A = A),
              prior = c(
                prior(beta(1, 1), "coi"),
                prior(student_t(3, 0, 2.5), "Intercept"),                                             
                prior(gamma(0.01, 0.01),"phi"),
                prior(student_t(3, 0, 2.5), "sd"),
                prior(beta(1, 1),"zoi")
              ))
summary(model1C)

plot7 <- plot(conditional_effects(model1C), points = FALSE, theme = a)

plot8 <- plot7$Host_Status + labs(x = "Migratory host category", y = "Weighted Betweenness") #caption = "Plasmodium")
plot8


save.image("Fidelity Data Final.RData") # Saving R Data


##### END #####
