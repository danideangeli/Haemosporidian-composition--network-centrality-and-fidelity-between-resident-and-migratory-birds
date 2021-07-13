library(bipartite)
#library(biGraph)
library(bipartiteD3)
#devtools::install_github("pedroj/bipartite_plots")
library(ggbipart)
library(igraph)
library(statnet)
library(vegan)
library(readxl)
library(tidyverse)
library(data.table)
library(network)
library(ggplot2)
library(sna)
library(intergraph)
#devtools::install_github("briatte/ggnet")
library(ggnet)
library(ergm)
library(intergraph)
library(RColorBrewer)
library(tidygraph)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3/Tabelas")

##### Data Manipulation - not in use anymore #####

dados <- read_excel("cap3.xlsx") 

dados$Loc <- paste(dados$Latitude,dados$Longitude, sep = "__") #creating coloumn with an ID for all sites
Loc.ID <- cbind.data.frame(Loc = unique(dados$Loc), ID = 1:length(unique(dados$Loc)))
dados <- merge.data.frame(dados, Loc.ID, by = "Loc")
rm(Loc.ID)

dados1 <- with(dados, names(table(Lineage_Name)[table(Lineage_Name) > 4])) 
dados1 <- dados[dados$Lineage_Name%in% dados1, ]
dados1$Lineage_Name <- factor(dados1$Lineage_Name) 

dados2 <- with(dados1, names(table(species)[table(species) > 4])) 
dados2 <- dados1[dados1$species%in% dados2, ]
dados2$species <- factor(dados2$species) 

dados3 <- dados2 %>%
  group_by(Lineage_Name, species, add = TRUE) %>%
  summarise(Abundance = n()) %>%
  arrange(Lineage_Name)
dados3 <- as.data.table(dados3)

##### Loadind Network Data #####

load("Netdata.RData")
dados3 <- as.data.table(Netdata)
rm(Netdata)
#write.csv2(dados3, "Network.csv")

##### Creating Occurence Matrix #####

Occ.Mat <- dcast(dados3, species ~ Lineage_Name, fun=sum, value.var = "n_interactions")
Occ.Mat <- column_to_rownames(Occ.Mat, var = "species")

Occ.Mat$Sum <- rowSums(Occ.Mat) # Creating a coloumn with the number of occurrences for each lineage
Occ.Mat <- arrange(Occ.Mat, desc(Sum)) # Ordering lineages from more sampled to less sampled
Occ.Mat <- filter(Occ.Mat, Sum >= 1 )
Occ.Mat$Sum <- NULL 

##### Plotting Matrixes and Nets #####

plotweb(Occ.Mat)
visweb(Occ.Mat, labsize = 3)

Occ.Mat1 <- as.matrix(Occ.Mat)
class(Occ.Mat1)
plotPAC(Occ.Mat1, plot.scale=1, fill.col="Blue", arrow.col="cyan", 
        circles=FALSE, radius=1, scaling = 0.2 )

mod <- computeModules(Occ.Mat)
plotModuleWeb(mod)                   

plotweb(Occ.Mat, method="cca", text.rot="90",labsize=0.5,col.low="green4", col.high="gold",
        col.interaction="skyblue3", bor.col.interaction ="blue")

gplot(as.one.mode(Occ.Mat), label=colnames(Occ.Mat), gmode="graph", 
      label.cex=0.6, vertex.cex=2)

Net <- bip_railway(Occ.Mat, label = T, nodesize = 2)
Net+ coord_flip()

Host_Status <- rownames_to_column(Occ.Mat)
Host_Status <- as.data.frame(Host_Status$rowname)
names(Host_Status)[names(Host_Status) == "Host_Status$rowname"] <- "species"
dados4<- unique(dados3[ , c("species", "Host_Status")])
Host_Status <- merge(Host_Status, dados4 , by = "species")

rede <- graph_from_incidence_matrix(Occ.Mat1, mode = "all", weighted = TRUE)
rede$Host_Status <- table(V(rede)$Host_Status)

coul  <- brewer.pal(3, "Set1") 
my_color <- coul[as.numeric(as.factor(V(rede1$Host_Status)))]

Occ.Mat2 <- Occ.Mat
Occ.Mat2 <- rownames_to_column(Occ.Mat2)
Occ.Mat2 <- inner_join(Occ.Mat2, Host_Status, by = c("rowname" = "species"))
Occ.Mat2 <- arrange(Occ.Mat2, desc(Host_Status))
Occ.Mat2 <- column_to_rownames(Occ.Mat2, var = "rowname")
Occ.Mat3 <- as.matrix(Occ.Mat2[,1:40])

bip.net<- bip_init_network(t(Occ.Mat3))

bip_ggnet(bip.net, Occ.Mat3, label = FALSE, shape= "mode",
          color= "mode", palette = "Set1", layout.exp = 0,
          edge.label = NULL, label.size = 3, size = 3)

x2 = data.frame(Type = network.vertex.names(bip.net))
x2 = factor(c(rep("Parasite", 40), rep("Resident", 227), rep("Partial Migrant",16), rep("Full Migrant",6)))

bip.net %v% "color" <- as.character(x2)
y2 <- RColorBrewer::brewer.pal(12, "Paired")[c(6,2,4,10)]
names(y2) <- levels(x2)
bip_ggnet(bip.net, Occ.Mat1, label = FALSE, shape= "mode",
          color= "color", palette = y2, layout.exp = 0,
          edge.label = NULL, label.size = 3, size = 3)

save.image("RedeImage.RData")

##### Calculating Indixes #####

networklevel(Occ.Mat, index=c("ISA", "weighted NODF", "Fisher alpha", "links per species", 
                              "connectance"), SAmethod="log")

grouplevel(Occ.Mat, level="both", 
           index=c("mean number of links", "weighted cluster coefficient", 
                   "effective partners", "niche overlap"), dist="bray")

SpProperties <- specieslevel(Occ.Mat, index=c("normalised degree", "PDI","effective partners", 
                              "species specificity"), PDI.normalise=F)

HostProperties <- SpProperties$`higher level`
