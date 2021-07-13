library(igraph) 
library(network) 
library(sna)
library(ggraph)
library(visNetwork)
library(threejs)
library(networkD3)
library(ndtv)
setwd("~/Lab Poulin/PhD/PhD/Dados/Capítulo 3")

##### Exemplo 1 #####

nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
head(nodes)
head(links)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net
?graph_from_data_frame

E(net)       # The edges of the "net" object
V(net)       # The vertices of the "net" object
E(net)$type  # Edge attribute "type"
V(net)$media # Vertex attribute "media"

# Find nodes and edges by attribute:
# (that returns oblects of type vertex sequence/edge sequence)
V(net)[media=="BBC"]
E(net)[type=="mention"]

# You can also examine the network matrix directly:
net[1,]
net[5,7]

# Get an edge list or a matrix:
as_edgelist(net, names=T)
as_adjacency_matrix(net, attr="weight")

# Or data frames describing nodes and edges:
as_data_frame(net, what="edges")
as_data_frame(net, what="vertices")

plot(net)
net <- simplify(net, remove.multiple = F, remove.loops = T) 
plot(net)
plot(net, edge.arrow.size=.2,vertex.label=NA)

##### Exemplo 2 #####

nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)
  
head(nodes2)
head(links2)

links2 <- as.matrix(links2)
dim(links2)
dim(nodes2)

head(nodes2)
head(links2)

net2 <- graph_from_incidence_matrix(links2)
table(V(net2)$type)

plot(net2, edge.arrow.size=.2,vertex.label=NA)

plot(net, edge.arrow.size=.2, edge.color="turquoise2", edge.curved=.1,
     vertex.color="turquoise2", vertex.frame.color="#ffffff",
     vertex.label=V(net)$media, vertex.label.color="black") 

# Compute node degrees (#links) and use that to set node size:
deg <- igraph::degree(net, mode="all")
V(net)$size <- deg*3

# We could also use the audience size value:
V(net)$size <- V(net)$audience.size*0.6

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label <- NA

# Set edge width based on weight:
E(net)$width <- E(net)$weight/6

#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "red"

# Generate colors based on media type:
colrs <- c("turquoise3", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]

# We can even set the network layout:
graph_attr(net, "layout") <- layout_with_lgl
plot(net, vertex.label=V(net)$media, vertex.label.color="black", edge.color="purple") 

legend(x=-2.2, y=-0.6, c("Newspaper","Television", "Online News"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=1, bty="n", ncol=1)

plot(net, vertex.shape="none", vertex.label=V(net)$media, 
     vertex.label.font=2.5, vertex.label.color="black",
     vertex.label.cex=1, edge.color="gray85")

# Add color to the edges based on their sources
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- V(net)$color[edge.start]

plot(net, edge.color=edge.col, edge.curved=.2)  

##### Layouts examples #####

net.bg <- sample_pa(100) 
V(net.bg)$size <- 8
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- "" 
E(net.bg)$arrow.mode <- 0
plot(net.bg)

plot(net.bg, layout=layout_randomly)

l <- layout_in_circle(net.bg)
plot(net.bg, layout=l)

l2 <- layout_randomly(net.bg)
plot(net.bg, layout=l2)

l3 <- layout_on_sphere(net.bg)
plot(net.bg, layout=l3)

l4 <- layout_with_fr(net.bg, niter=100)
plot(net.bg, layout=l4)

ws  <-  c(1, rep(100, ecount(net.bg)-1))
lw <- layout_with_fr(net.bg, weights=ws)
plot(net.bg, layout=lw) 

par(mfrow=c(2,2), mar=c(0,0,0,0))   # plot four figures - 2 rows, 2 columns
plot(net.bg, layout=layout_with_fr)
plot(net.bg, layout=layout_with_fr)
plot(net.bg, layout=l4)
plot(net.bg, layout=l4)

l <- layout_with_fr(net.bg)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(2,2), mar=c(0,0,0,0))
plot(net.bg, rescale=F, layout=l*0.4)
plot(net.bg, rescale=F, layout=l*0.6)
plot(net.bg, rescale=F, layout=l*0.8)
plot(net.bg, rescale=F, layout=l*1.0)

dev.off()

l <- layout_with_fr(net.bg, dim=3)
plot(net.bg, layout=l)

l <- layout_with_kk(net.bg)
plot(net.bg, layout=l)

l <- layout_with_graphopt(net.bg)
plot(net.bg, layout=l)

# The charge parameter below changes node repulsion:
l1 <- layout_with_graphopt(net.bg, charge=0.02)
l2 <- layout_with_graphopt(net.bg, charge=0.00000001)

par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(net.bg, layout=l1)
plot(net.bg, layout=l2)

dev.off()

plot(net.bg, layout=layout_with_lgl)

plot(net.bg, layout=layout_with_mds)

# All layouts in one graph

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow=c(4,4), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net)) 
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }

dev.off()

# Selecting only the most frequent connections

hist(links$weight)
mean(links$weight)
sd(links$weight)

cut.off <- mean(links$weight) 
net.sp <- delete_edges(net, E(net)[weight<cut.off])
plot(net.sp, layout=layout_with_kk) 

par(mfrow=c(1,2))

# Community detection (by optimizing modularity over partitions):
clp <- cluster_optimal(net)
class(clp)

# Community detection returns an object of class "communities" 
# which igraph knows how to plot: 
plot(clp, net)

# We can also plot the communities without relying on their built-in plot:
V(net)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])
dev.off()
