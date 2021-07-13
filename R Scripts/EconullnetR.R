devtools::install_github("ivaughan/econullnetr")
library(econullnetr)

set.seed(1234)
view(Silene)
sil.null <- generate_null_net(Silene[, 2:7], Silene.plants[, 2:6], sims = 10,
                              c.samples = Silene[, 1],
                              r.samples = Silene.plants[, 1], prog.count = FALSE)


# Network-level statistics
bipartite_stats(sil.null, index.type = "networklevel",
                indices = c("linkage density", "weighted connectance", 
                            "interaction evenness"), intereven = "sum", 
                prog.count = FALSE)
