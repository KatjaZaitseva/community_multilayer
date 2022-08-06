Community algorithms comparison
================
Ekaterina Zaitseva
7/3/2022

# Data Wrangling

Before starting the data cleaning process, we load all required packages
in RStudio.

<!-- Tips to make rmarkdown nice: https://github.com/rstudio/rmarkdown 
https://github.com/bbest/rmarkdown-example-->

``` r
library(multinet)
library(bipartite)
library(dplyr)
library(reshape2)
library(tidyverse)
library(igraph)
#the library for multilevel
library(devtools, quietly = TRUE)
library(Matrix)
library(BBmisc)
library(aricode)

library(muxViz)

library(networkdata)

require(gridExtra)
library(scales)

library(blockmodeling)
library(R.utils)
```

# Functions

## Community detection comparison output metrics

``` r
num_communities <- function(x) max(x$cid) + 1
grouped_actors <- function(x) x %>% group_by(cid) %>% summarise(count = n_distinct(actor))
grouped_layers <- function(x) x %>% group_by(cid) %>% summarise(count = n_distinct(layer))
avg_community_size <- function(x) length(x[,1])/num_communities(x)
avg_actors_per_community <- function(x) mean(grouped_actors(x)$count)
avg_layers_per_community <- function(x) mean(grouped_layers(x)$count)
num_clustered_vertices <- function(x) length(unique(x[, c(1, 2)])[,1])
num_clustered_actors <- function(x) length(unique(x$actor))
actor_overlapping <- function(x) length(unique(x[, c(1, 3)])[, 1]) / length(unique(x[, 1]))
```

## Direct blockmodeling for multilevel structure

``` r
multiBlockmodeling <- function(net, n_k) {
  net_edges <- data.frame(nodeID_from = edges_ml(net)$from_actor, 
                          nodeID_to = edges_ml(net)$to_actor,
                          layerLabel = edges_ml(net)$from_layer, weight = 1)
  net_layers <- data.frame(layerID = seq(1:length(layers_ml(net))), 
                           layerLabel = layers_ml(net))
  net_edges_l <- merge(net_edges, net_layers, by = 'layerLabel')
  node_to <- c(unlist(sort(unique(net_edges_l['nodeID_to']))))
  node_from <- c(unlist(sort(unique(net_edges_l['nodeID_from']))))
  
  #extend dataset to all nodes to make symmetric matrix with nodes for every layer
  node_diff <- c(setdiff(node_from,node_to),setdiff(node_to,node_from))
  if (length(node_diff) > 0){
    net_edges_extended <- net_edges_l
    for(i in 1:length(node_diff)){
      for(j in 1:3){
        net_edges_extended <- net_edges_extended %>% 
          add_row(layerID = j, nodeID_from = node_diff[i],
                  nodeID_to = node_diff[i], weight=0, layerLabel = net_layers[j,2])
      }
    }
  }
  #make 3D-array
  if (length(node_diff) > 0){
    phys_3d_array <- xtabs(weight ~ nodeID_from + nodeID_to + layerID, net_edges_extended)
  }else{
    phys_3d_array <- xtabs(weight ~ nodeID_from + nodeID_to + layerID, net_edges_l)
  }
  #blockmodeling
  #blocks = com - complete block matrix
  #hom - homogeneity blockmodeling 
  #and then for hom homFun = ss -  sum of squares homogeneity blockmodeling
  res <- optRandomParC(M = phys_3d_array, k = n_k, rep = 30, printRep = FALSE,
                       approaches = "hom", homFun = "ss", blocks = "com",
                       seed = 42)
  plot(res)
  #clustering results
  cid <- res$best$best1$clu - 1
  actor <- sort(actors_ml(net)$actor)
  block_comm <- data.frame(cid=cid, actor=actor)
  net_vertices <- data.frame(actor=vertices_ml(net)$actor, layer=vertices_ml(net)$layer)
  
  block_comm <- merge(x = block_comm, y = net_vertices[ , c("actor", "layer")], by = "actor", all.x=TRUE)
  block_comm <- block_comm[,c('actor', 'layer', 'cid')]
  return(block_comm)
}
```

## Spectral clustering on Multi-Layer graphs with k-means

The code adapted with the help of this
[resource](https://github.com/justin830827/Community-Detection-in-Multilayer-Graph.git).

Create list of Laplacian matrices

``` r
getLapList <- function(net, norm = TRUE, sparse = FALSE) {
  gen_lap_list <- list()
  for (l in layers_ml(net)){
    gen_l_igraph <- as.igraph(net, layers = c(l), 
                              merge.actors = TRUE, 
                              all.actors = FALSE)
    gen_l_simple <- igraph::simplify(gen_l_igraph, remove.loops = TRUE)
    gen_l_simple <- set_graph_attr(gen_l_simple, "layout", layout_with_kk(gen_l_simple))
    lap_matrix <- laplacian_matrix(gen_l_simple, norm=TRUE, sparse=FALSE)
    gen_lap_list <- append(gen_lap_list, list(lap_matrix))
  }
  return(gen_lap_list)
}
```

Create a special matrix U

``` r
getU <- function(lap, k) {
  rowlap <- rownames(lap)
  ei <- eigen(Matrix(lap))
  ei_val <- ei$values
  ei_vec <- ei$vectors
  par(mfrow=c(1,1))
  plot(ei_val, type="b", main="Eigenvalues")
  
  ind_smallest <- tail(order(ei_val, decreasing = TRUE),k)
  vec_smallest <- ei_vec[,c(ind_smallest)]
  rownames(vec_smallest) <- rowlap
  
  #Matrix(Re(vec_smallest), sparse=FALSE)
  as.matrix(Re(vec_smallest), sparse=FALSE)
}
```

The function for the algorithm

``` r
SCML <- function(net, k, alpha) {
  lap_list <- getLapList(net)
  U_list <- list()
  uu_mul <- list()
  for (i in 1:length(lap_list)){
    l <- lap_list[[i]]
    U <- getU(l, k)
    U_list[[i]] <- as.matrix(U)
  }
  #print(lap_list)
  #print(U_list)
  
  net_all_rows <- unique(unlist(lapply(lap_list, rownames)))
  net_all_cols <- unique(unlist(lapply(lap_list, colnames)))
  
  n <- length(net_all_cols)
  #print(n)
  
  for (i in 1:length(U_list)){
    u <- U_list[[i]]
    uu <- u %*% t(u)
    uu_mul[[i]] <- uu
  }
  lap_sum <- matrix(0, nrow = length(net_all_rows), ncol = length(net_all_cols))
  rownames(lap_sum) <- net_all_rows
  colnames(lap_sum) <- net_all_cols
  
  uu_sum <- matrix(0, nrow = length(net_all_rows), ncol = length(net_all_cols))
  rownames(uu_sum) <- net_all_rows
  colnames(uu_sum) <- net_all_cols
  
  for (l in lap_list){
    #matrices can be diff dim
    lap_sum[rownames(l), colnames(l)] = lap_sum[rownames(l), colnames(l)] + l
  }
  for (u in uu_mul){
    #matrices can be diff dim
    uu_sum[rownames(u), colnames(u)] <- uu_sum[rownames(u), colnames(u)] + u
  }
  
  #print(lap_sum)
  #print(uu_sum)
  
  lap_mod <- lap_sum - (alpha * uu_sum)
  lap_mod <- as.matrix(lap_mod)
  
  #print(lap_mod)
  
  U <- getU(lap_mod, k)
  #print(U)
  U_norm <- BBmisc::normalize(U, method = "range", range = c(0, 1))
  #print(U_norm)
  
  clusters_net <- kmeans(U_norm, k, nstart=50, iter.max=500)
  k_means_com <- data.frame(actor=row.names(lap_mod),cid=as.numeric(clusters_net$cluster)-1)
  net_df <- data.frame(actor=vertices_ml(net)$actor,layer=vertices_ml(net)$layer)
  k_means_com <- merge(x = k_means_com, y = net_df[ , c("actor", "layer")], by = "actor", all.x=TRUE)
  k_means_com <- k_means_com[,c('actor', 'layer', 'cid')]
  
  return(k_means_com)
}
```

You can also embed plots, for example:

![](community_algorithms_comparison_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
