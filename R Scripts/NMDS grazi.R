library(vegan)
library(gdata)
library(effects)
library(ggplot2)
library(pca3d)
table <-read.table("Comp.txt", h=T, check.names=FALSE, row.names = 1)
table<-t(table)
meta_table <- read.table("Comp.txt", h=T, check.names=FALSE)

sol <- metaMDS(table)
par (mfrow = c(1,2))
ordiplot (sol, cex = 1.0, type = 't')
stressplot (sol, cex = 1.0)

table <-read.table("Comp.txt", h=T, check.names=FALSE, row.names = 1)
table2 <-read.table("Comp2.txt", h=T, check.names=FALSE, row.names = 1)

attach(table2)

dis <- vegdist(table,method="bray")
mod.coleta <- betadisper(dis, Coleta ,type="centroid")
anova(mod.coleta)

permutest(mod.coleta, pairwise = T, permutations = 999)

mod.coleta.HSD <- TukeyHSD(mod.coleta)
plot(mod.coleta.HSD)
plot(mod.coleta)
Stangle("Comp.txt")
