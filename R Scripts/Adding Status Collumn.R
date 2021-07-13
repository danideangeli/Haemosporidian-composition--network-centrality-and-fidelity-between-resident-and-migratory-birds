library(tidyverse)
library(readxl)

setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3")

dados <- read_excel("PhD2021.xlsx")

unique(dados$Status)

dados1 <- filter(dados, Status != "Missing")

dados2 = dados1 %>%
  group_by(Especie)%>%
  group_by(Status, add = TRUE) %>%
  summarise(n_individuals = n()) %>%
  arrange(Especie)

dados2$n_individuals <- NULL

dados3 <- filter(dados, Status == "Missing")

dados4 <- left_join(dados3, dados2, by = c("Especie" = "Especie"))

length(unique(dados4$Status.y))
unique(dados4$Status.y)

write.csv2(dados4, "PhD2021A.csv")

dados <- read_excel("PhD2021.xlsx")

dados4 <- filter(dados, Temp == "0")

library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=10)

r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

lats <- dados$Latitude
lons <- dados$Longitude
lats <- as.numeric(lats)
lons <- as.numeric(lons)

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)

write.csv2(df, "missingworldclim.csv")
