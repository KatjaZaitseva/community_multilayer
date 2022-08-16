Community algorithms comparison
================
Ekaterina Zaitseva
7/3/2022

# Data Wrangling

Before starting the data cleaning process, we load all required packages
in RStudio.

<!-- Tips to make rmarkdown nice: https://github.com/rstudio/rmarkdown 
https://github.com/bbest/rmarkdown-example-->

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

1.  Create list of Laplacian matrices

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

2.  Create a special matrix U

``` r
getU <- function(lap, k, plot=TRUE) {
  rowlap <- rownames(lap)
  ei <- eigen(Matrix(lap))
  ei_val <- ei$values
  ei_vec <- ei$vectors
  if (plot){
    par(mfrow=c(1,1))
    plot(ei_val, type="b", main="Eigenvalues")
  }
  
  ind_smallest <- tail(order(ei_val, decreasing = TRUE),k)
  vec_smallest <- ei_vec[,c(ind_smallest)]
  rownames(vec_smallest) <- rowlap
  
  #Matrix(Re(vec_smallest), sparse=FALSE)
  as.matrix(Re(vec_smallest), sparse=FALSE)
}
```

3.  The function for the Spectral Clustering algorithm on Multi-Layer
    graphs (SC-ML)

``` r
SCML <- function(net, k, alpha, plotShow=TRUE) {
  lap_list <- getLapList(net)
  U_list <- list()
  uu_mul <- list()
  for (i in 1:length(lap_list)){
    l <- lap_list[[i]]
    U <- getU(l, k, plot=plotShow)
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
  
  U <- getU(lap_mod, k, plot=plotShow)
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

## Tunning parameters for SCML (k, alpha)

``` r
tuning_k <- function(net, alpha, true_labels=FALSE, labels=NULL) {
  kl <- seq(2, 15, 1)
  nmil <- modl <- list()
  for (k in 1:length(kl)){
    alg <- SCML(net, kl[k], alpha, plotShow=FALSE)
    if (true_labels){
      metrics <- evaluation_SCML(alg, net, labels_col = labels)
      nmil <- append(nmil, metrics[1])
    }else{
      modl <- append(modl, modularity_ml(net, alg, 
                                         gamma = 1, omega = 1))
    }
  }
  if (true_labels){
    k_df <- data.frame(unlist(kl), unlist(nmil))
    names(k_df) <- c("k", "nmi")
    par(mfrow=c(1,1))
    #plot nmi
    ggplot(k_df, aes(x = k, y = nmi)) +
      geom_line(linetype = "dashed") +
      geom_point(size=2) + 
      labs(title = "SCML k tuning",
           y = "NMI", x = "# clusters")
  }else{
    k_df <- data.frame(unlist(kl), unlist(modl))
    names(k_df) <- c("k", 'mod')
    par(mfrow=c(1,1))
    #plot nmi
    ggplot(k_df, aes(x = k, y = mod)) +
      geom_line(linetype = "dashed") +
      geom_point(size=2) + 
      labs(title = "SCML k tuning",
           y = "Modularity", x = "# clusters")
  }
}

tuning_alpha <- function(net, k, true_labels=FALSE, labels=NULL) {
  alpha <- seq(0.1, 1.2, 0.1)
  nmil <- modl <- list()
  for (a in 1:length(alpha)){
    alg <- SCML(net, k, alpha[a], plotShow=FALSE)
    if (true_labels){
      metrics <- evaluation_SCML(alg, net, labels_col = labels)
      nmil <- append(nmil, metrics[1])
    }else{
      modl <- append(modl, modularity_ml(net, alg, 
                                         gamma = 1, omega = 1))
    }
  }
  if (true_labels){
    alpha_df <- data.frame(unlist(alpha), unlist(nmil))
    names(alpha_df) <- c("alpha", "nmi")
    par(mfrow=c(1,1))
    #plot nmi
    ggplot(alpha_df, aes(x = alpha, y = nmi)) +
      geom_line(linetype = "dashed") +
      geom_point(size=2) + 
      labs(title = "SCML alpha tuning",
           y = "NMI", x = "alpha")
  }else{
    alpha_df <- data.frame(unlist(alpha), unlist(modl))
    names(alpha_df) <- c("alpha", "mod")
    par(mfrow=c(1,1))
    #plot nmi
    ggplot(alpha_df, aes(x = alpha, y = mod)) +
      geom_line(linetype = "dashed") +
      geom_point(size=2) + 
      labs(title = "SCML alpha tuning",
           y = "Modularity", x = "alpha")
  }
}
```

## Accuracy evaluation metrics

``` r
purity_score <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

evaluation <- function(comm_str, net, labels_col, gen=FALSE, gen_true=NULL, 
                            NA_flag=FALSE, onelayer=FALSE) {
  a <- actors_ml(net)$actor
  v <- vertices_ml(net)
  
  if (gen){
    net_df_com <- gen_true
  }else{
    true_labels <- get_values_ml(net, attribute=labels_col, actors=a)
    
    if (NA_flag){
      net_df <- data.frame(actor=a, group=true_labels[,1])
      net_df <- net_df %>% subset(group!="NA") %>%
        #remain 1 last character -> -1
        mutate(cid = substr(group, nchar(group)-1+1, nchar(group))) %>%
        mutate_at(c('cid'), as.numeric) %>% 
        select(c('actor','cid'))
    }
    else{
      net_df <- data.frame(actor=a, cid=true_labels[,1])
      net_df <- net_df %>%
        mutate_at(c('cid'), as.numeric)
    }
    net_layer_df <- data.frame(actor=v$actor,layer=v$layer)
    net_df_com <- merge(x = net_df, y = net_layer_df[ , c("actor", "layer")], by = "actor", all.x=TRUE)
    net_df_com <- net_df_com[,c('actor', 'layer', 'cid')]
  }
  if (NA_flag){
    comm_str <- comm_str[comm_str$actor %in% net_df_com$actor, , drop = FALSE]
    comm_str <- comm_str[order(comm_str$actor),]
  }else{
    comm_str <- comm_str[order(comm_str$actor),]
  }
  net_df_com <- net_df_com[order(net_df_com$actor),]
  
  #nmi <- NMI(comm_str$cid, net_df_com$cid, variant='sum')
  if (dim(comm_str)[1]==dim(net_df_com)[1]){
    if (onelayer) {
      NMI(comm_str$cid, net_df_com$cid, variant='sum')
    }else{
      nmi <- nmi_ml(net, comm_str, net_df_com)
    }
    comm_str_agg <- comm_str %>% group_by(actor, cid) %>% arrange(actor)
    net_df_com_agg <- comm_str %>% group_by(actor, cid) %>% arrange(actor)
    pur <- purity_score(comm_str$cid, net_df_com$cid)
    randi <- rand.index(comm_str$cid, net_df_com$ci)
  } else {
    pur <- NA
    nmi <- NA
    randi <- NA
  }
  
  return(c(nmi,pur))
  # cat(" NMI", nmi_ml(net, kmeans_com, net_df_com), "\n",
  #     "Purity score", purity_score(kmeans_com$cid, net_df_com$cid), "\n",
  #     "Modularity", modularity_ml(net, kmeans_com, gamma = 1, omega = 1))
}
```

## Multinet built-in algorithms + 2 new (SC-ML and Blockmodeling)

``` r
multilayer_algorithms <- function(net, k, alpha, p.min.actors=3, p.min.layers=2,
                                  p.gamma=0.8, p.omega=1){
  sys_time_a <- list()
  #flattening
  #weighted
  start_time <- Sys.time()
  c1 <- flat_ec_ml(net)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #unweighted
  start_time <- Sys.time()
  c2 <- flat_nw_ml(net)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #layer by layer
  start_time <- Sys.time()
  c3 <- abacus_ml(net, min.actors = p.min.actors,  min.layers = p.min.layers)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #multilayer
  #density
  #k - minimum number of actors in a clique
  #m - minimum number of common layers in a clique
  start_time <- Sys.time()
  c4 <- clique_percolation_ml(net, k=p.min.actors, m=p.min.layers)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #optimization
  #gamma - the Kronecker delta
  #omega - inter-layer weight parameter
  #higher values in omega  will result in communities spanning multiple layers, 
  #because inlcuding the same actor on different layers in the same community 
  #increases the value of modularity
  start_time <- Sys.time()
  c5 <- glouvain_ml(net, gamma=p.gamma, omega=p.omega)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #random walks
  start_time <- Sys.time()
  c6 <- infomap_ml(net, overlapping=FALSE, directed=FALSE, self.links=TRUE)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  #label propagation
  start_time <- Sys.time()
  c7 <- mdlp_ml(net)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  
  start_time <- Sys.time()
  c8 <- SCML(net, k, alpha)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  
  start_time <- Sys.time()
  withTimeout(Sys.sleep(10), timeout = 30)
  c9 <- multiBlockmodeling(net, k)
  end_time <- Sys.time()
  sys_time_a <- c(sys_time_a, round(end_time - start_time, 2))
  
  alg_list <- list()
  for (i in 1:9){
    alg_list[[i]] <- eval(parse(text = paste0('c',i)))
  }
  
  return(list(alg_list, sys_time_a))
}
```

## Algorithms comparison

``` r
algorithms_comparison <- function(net, alg_list, sys_time_l, mod.gamma=1, mod.omega=1,
                                  true_labels=FALSE, labels=NULL, genBool=FALSE,
                                  gen_truel=NULL, flag_NA=FALSE){
  #calculate
  #(1) the number of communities generated
  #(2) the average community size
  #(3) the proportion of vertices included in at least one cluster 
  #(which is 1 for complete community detection methods)
  #(4) the proportion of actors included in at least one cluster 
  #(which is 1 for complete community detection methods)
  #(5) the ratio between the number of actor-community pairs and the number of clustered actors, 
  #indicating the level of overlapping (which is 1 for partitioning 
  #community detection methods and higher for overlapping methods)
  
  #the function is adapted from https://bitbucket.org/uuinfolab/20csur/src/master/R/summary.R
  num_com <- avg_act <- avg_lay <- modul <- nmil <- purl <- sysl <-  clust_actorsl <- list()
  #randil <- list()
  for (i in 1:length(alg_list)){
    alg <- alg_list[[i]]
    num_com <- append(num_com, num_communities(alg))
    avg_act <- append(avg_act, avg_actors_per_community(alg))
    avg_lay <- append(avg_lay, avg_layers_per_community(alg))
    clust_actorsl <-  append(clust_actorsl, round(num_clustered_actors(alg) / num_actors_ml(net),2))
    modul <- append(modul, modularity_ml(net, alg, 
                                         gamma = mod.gamma, omega = mod.omega))
    if (true_labels || genBool){
      metrics <- evaluation(alg, net, labels_col = labels, gen = genBool, gen_true = gen_truel, NA_flag = flag_NA)
      nmil <- append(nmil, metrics[1])
      purl <- append(purl, metrics[2])
      #randil <- append(randil, metrics[3])
    }
  }
  
  if (true_labels || genBool){
    com_stats <- data.frame(num_communities = unlist(num_com), avg_actors_per_community = unlist(avg_act),
                            clustered_actors =  unlist(clust_actorsl), modularity =  unlist(modul),
                            nmi = unlist(nmil), accuracy = unlist(purl), exec_time = unlist(sys_time_l))
    row.names(com_stats) <- c("flatten_weighted","flatten_unweighted",
                              "abacus", "clique p.", "louvain", "infomap","mdlp","SCML", "Blockm")
  } else {
    com_stats <- data.frame(num_communities = unlist(num_com), avg_actors_per_community = unlist(avg_act),
                            clustered_actors =  unlist(clust_actorsl), modularity =  unlist(modul), 
                            exec_time = unlist(sys_time_l))
    row.names(com_stats) <- c("flatten_weighted","flatten_unweighted",
                              "abacus", "clique p.", "louvain", "infomap","mdlp","SCML","Blockm")
  }
  print(paste("The algorithms comparison for the network:",deparse(substitute(net))))
  print(round(com_stats,2))
}
```

## Algorithm SC-ML for one layer

``` r
onelayerSCML <- function(net, k, layer, drawGraph = FALSE) {
  gen_l_igraph <- as.igraph(net, layers = c(layer), 
                            merge.actors = TRUE, 
                            all.actors = FALSE)
  gen_l_simple <- igraph::simplify(gen_l_igraph, remove.loops = TRUE)
  gen_l_simple <- set_graph_attr(gen_l_simple, "layout", layout_with_kk(gen_l_simple))
  lap_matrix <- laplacian_matrix(gen_l_simple, norm=TRUE, sparse=FALSE)
  U <- getU(lap_matrix, k, plot=FALSE)
  U_norm <- BBmisc::normalize(U, method = "range", range = c(0, 1))
  
  clusters_net <- kmeans(U_norm, k, nstart=50, iter.max=500)
  k_means_com <- data.frame(actor=row.names(lap_matrix),cid=as.factor(clusters_net$cluster))
  if(drawGraph){
    par(mfrow = c(1, 1))
    color_map <- brewer.pal(max(k_means_com$cid), "Set3") 
    plot(k_means_com, gen_l_igraph, vertex.label=NA, vertex.size=8,
         vertex.color=color_map, edge.color = 'grey',
         mark.col=NA, mark.border=NA) #remove polygon
  }else{
    return(k_means_com)
  }
}
```

# Real network data sets

## Data set 1: Krackhardt High-Tech Managers

Load the network Krackhardt-High-Tech

``` r
tech_nodes <- read.table("Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_nodes.txt", header=TRUE)
tech_layers <- read.table("Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_layers.txt", header=TRUE)
tech_layers$layerLabel[tech_layers$layerLabel == "Reports_to"] <- "reportsto"
tech_links <- read.table("Krackhardt-High-Tech_Multiplex_Social/Dataset/Krackhardt-High-Tech_multiplex.edges", 
                         col.names = c('layerID', 'nodeID_from', 'nodeID_to', 'weight'))
tech_links <- merge(x = tech_links, y = tech_layers, by = 'layerID', all.x = TRUE)
tech_links_2 <- tech_links[, c(5, 2, 3, 4)]
colnames(tech_links_2) <- c("layer_from","nodeID_from","nodeID_to","weight")
tech_links_2['layer_to'] <- tech_links_2['layer_from']

krack_net <- ml_empty()

#add layers
add_layers_ml(krack_net, unlist(tech_layers[2]), directed=FALSE)

#add vertices
by_nodefrom_layer <- tech_links_2 %>% group_by(nodeID_from, layer_from) %>% summarise(n = n())
by_nodeto_layer <- tech_links_2 %>% group_by(nodeID_to, layer_from) %>% summarise(n = n())
colnames(by_nodefrom_layer) <- c("nodeID", "layer", "n")
colnames(by_nodeto_layer) <- c("nodeID", "layer", "n")
by_node_layer <- rbind(by_nodefrom_layer, by_nodeto_layer)
by_node_layer_agg <- by_node_layer %>% group_by(nodeID, layer) %>% summarise(n = n())
vertices <- data.frame(actors = by_node_layer_agg['nodeID'],
                       layers = by_node_layer_agg['layer'])
add_vertices_ml(krack_net, vertices)

#add edges
intra_layer_edges <- data.frame(actors_from = tech_links_2['nodeID_from'],
                                layers_from = tech_links_2['layer_from'],
                                actors_to = tech_links_2['nodeID_to'],
                                layers_to = tech_links_2['layer_to'])
colnames(intra_layer_edges) <- c('actor1','layer1','actor2','layer2')
add_edges_ml(krack_net, intra_layer_edges)

#attributes
add_attributes_ml(krack_net, colnames(tech_nodes[c(2:5)]), type='numeric')
#add age
set_values_ml(krack_net, "nodeAge", unlist(tech_nodes[1]),
              values = unlist(tech_nodes["nodeAge"]))
#add tenure
set_values_ml(krack_net, "nodeTenure", unlist(tech_nodes[1]), values = unlist(tech_nodes["nodeTenure"]))
#add level
set_values_ml(krack_net, "nodeLevel", unlist(tech_nodes[1]), values = unlist(tech_nodes["nodeLevel"]))
#add department
set_values_ml(krack_net, "nodeDepartment", unlist(tech_nodes[1]), values = unlist(tech_nodes["nodeDepartment"]))
```

## Data set 2: CS-AARHUS

``` r
aucs_net <- ml_aucs()
```

# Plot networks

<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->

A multilayer network Krackhardt High-Tech Managers representing 3 types
of interactions between 21 managers

``` r
par(mfrow = c(1, 1))

#plot(krack_net)
#plot aligned network between layers with layer 'Advice'
#palette(brewer.pal(n = 8, name = "Accent"))
palette('default')
#align all the actors (w_inter=1) with respect to their layout in the 1 layer
l_krack <- layout_multiforce_ml(krack_net, w_inter = 1, w_in = c(1, 0, 0),
                          gravity = c(1, 0, 0))

par(mfrow = c(1, 1)) 

#plot with level
level <- get_values_ml(krack_net, actors = vertices_ml(krack_net)[[1]],
                            attribute = "nodeLevel")
gr <- values2graphics(level)
#remove the labels (vertex.labels = "")
plot(krack_net, layout = l_krack, grid = c(1, 3), vertex.labels = "",
     vertex.color = gr$color, edge.col = "gray40")
#add a legend with the names of the layers
legend("bottom", legend = c("CEO","Vice President",'Manager'), title = 'Level',
       col = gr$legend.col, pt.bg = gr$legend.col, pch = gr$legend.pch, bty = "n", pt.cex = 1,
       cex = 0.7, inset = c(0, -0.15), xpd=TRUE, horiz = TRUE)
```

![](community_algorithms_comparison_files/figure-gfm/plot%20krack-1.png)<!-- -->

A multilayer network CS-AARHUS representing 6 types of interactions
between 61 employees at a university research department

``` r
#plot aligned network between layers with layer 'work'
l_aucs <- layout_multiforce_ml(aucs_net, w_inter = 1, w_in = c(0, 0, 0, 0, 1),
                                gravity = c(0, 0, 0, 0, 1))
#plot with role
role <- get_values_ml(aucs_net, actors = vertices_ml(aucs_net)[[1]],
                       attribute = "role")
role <- as.factor(role[[1]])
gr <- values2graphics(role)
#remove the labels (vertex.labels = "")
plot(aucs_net, layout = l_aucs, grid = c(1, 5), vertex.labels = "",
     vertex.color = gr$color, edge.col = "gray40")
#add a legend with the names of the layers
legend("bottom", legend = c("Admin","Assistant","Associate","Emeritus","NA","PhD", 
                            "Phd (visiting)","Postdoc","Professor"), title = 'Level',
       col = gr$legend.col, pt.bg = gr$legend.col, pch = gr$legend.pch, bty = "n", pt.cex = 1,
       cex = 0.6, inset = c(0, -0.15), xpd=TRUE, horiz = TRUE)
```

![](community_algorithms_comparison_files/figure-gfm/plot%20aucs-1.png)<!-- -->

# Measuring networks

-   dens - the ratio between the number of edges and the number of
    possible edges
-   cc - clustering coefficient. The ratio between the triangles and the
    connected triples in the layer
-   apl - average path length. The average graph-distance between all
    pairs of vertices in the layer

``` r
round(summary(krack_net),2) 
```

    ##             n   m dir nc slc dens   cc  apl dia
    ## _flat_     21 244   0  1  21 1.16 0.79 1.24   2
    ## advice     21 145   0  1  21 0.69 0.73 1.31   2
    ## friendship 21  79   0  1  21 0.38 0.47 1.65   3
    ## reportsto  21  20   0  1  21 0.10 0.00 2.98   4

``` r
round(summary(aucs_net),2) 
```

    ##           n   m dir nc slc dens   cc  apl dia
    ## _flat_   61 620   0  1  61 0.34 0.48 2.06   4
    ## coauthor 25  21   0  8   6 0.07 0.43 1.50   3
    ## facebook 32 124   0  1  32 0.25 0.48 1.96   4
    ## leisure  47  88   0  2  44 0.08 0.34 3.12   8
    ## lunch    60 193   0  1  60 0.11 0.57 3.19   7
    ## work     60 194   0  1  60 0.11 0.34 2.39   4

# Layers analysis

-   jaccard.actors helps to answer the question if all actors are
    present on all layers, having a value equaled to 1 when the answer
    is yes
-   pearson.degree helps to answer the question whether actors having a
    high (or low) degree on one layer and behave similarly in other
    layers. The smallest value (−1) indicates that high-degree actors in
    one layer are low-degree in the other and vice versa, while the
    largest value (1) is returned if high-degree (respectively,
    low-degree) actors in one layer are high-degree (respectively,
    low-degree) actors in the other

``` r
round(layer_comparison_ml(krack_net, method="jaccard.actors"),2)
```

    ##            advice friendship reportsto
    ## advice          1          1         1
    ## friendship      1          1         1
    ## reportsto       1          1         1

``` r
round(layer_comparison_ml(aucs_net, method="jaccard.actors"),2)
```

    ##          lunch facebook coauthor leisure work
    ## lunch     1.00     0.53     0.42    0.78 0.97
    ## facebook  0.53     1.00     0.30    0.52 0.53
    ## coauthor  0.42     0.30     1.00    0.41 0.42
    ## leisure   0.78     0.52     0.41    1.00 0.78
    ## work      0.97     0.53     0.42    0.78 1.00

``` r
round(layer_comparison_ml(aucs_net, method="jaccard.actors"),2)
```

    ##          lunch facebook coauthor leisure work
    ## lunch     1.00     0.53     0.42    0.78 0.97
    ## facebook  0.53     1.00     0.30    0.52 0.53
    ## coauthor  0.42     0.30     1.00    0.41 0.42
    ## leisure   0.78     0.52     0.41    1.00 0.78
    ## work      0.97     0.53     0.42    0.78 1.00

``` r
round(layer_comparison_ml(aucs_net, method = "pearson.degree"),2)
```

    ##          lunch facebook coauthor leisure work
    ## lunch     1.00     0.31     0.15    0.28 0.25
    ## facebook  0.31     1.00     0.55    0.38 0.54
    ## coauthor  0.15     0.55     1.00    0.48 0.43
    ## leisure   0.28     0.38     0.48    1.00 0.07
    ## work      0.25     0.54     0.43    0.07 1.00

Actors analysis for krack network

``` r
#degree of a vertex - number of adjacent edges
par(mfrow = c(2, 2))
hist(degree_ml(krack_net), 
     breaks = max(degree_ml(krack_net)),
     main = "flattened",
     xlab = "degree")

color_map <- brewer.pal(num_layers_ml(krack_net), "Accent")
for (i in 1:num_layers_ml(krack_net)) {
  d <- degree_ml(krack_net, layers = layers_ml(krack_net)[[i]])
  hist(d, breaks = max(d, na.rm = TRUE), main = layers_ml(krack_net)[[i]],
       xlab = "degree", col = color_map[i])
}
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
par(mfrow = c(1, 1)) 
```

Actors analysis for aucs network

``` r
par(mfrow = c(2, 3))
hist(degree_ml(aucs_net), 
     breaks = max(degree_ml(aucs_net)),
     main = "flattened",
     xlab = "degree")

color_map <- brewer.pal(num_layers_ml(aucs_net), "Accent")
for (i in 1:num_layers_ml(aucs_net)) {
  d <- degree_ml(aucs_net, layers = layers_ml(aucs_net)[[i]])
  hist(d, breaks = max(d, na.rm = TRUE), main = layers_ml(aucs_net)[[i]],
       xlab = "degree", col = color_map[i])
}
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
par(mfrow = c(1, 1)) 
```

# Generated network data sets

Definitions \* in pillar community structures each actor belongs to the
same community on all layers \* in semipillar community structures the
communities in one layer are different from the other layers \* in
partitioning community structures each vertex belongs to one community
\* in overlapping community structures some vertices belong to multiple
communities

Types \* PEP (pillar partitioning) \* PEO (pillar overlapping) \* SEP
(semipillar partitioning) \* SEO (semipillar overlapping)

``` r
gen_net33_3_3 <- generate_communities_ml("pep", num.actors=33, 
                                   num.layers=3, num.communities=3,
                                   overlap=0)
gen_net40_3_8 <- generate_communities_ml("pep", num.actors=40, 
                                      num.layers=3, num.communities=8,
                                      overlap=0)
gen_net99_3_3 <- generate_communities_ml("pep", num.actors=99, 
                                    num.layers=3, num.communities=3,
                                    overlap=0)
gen_net104_3_8 <- generate_communities_ml("pep", num.actors=104, 
                                      num.layers=3, num.communities=8,
                                      overlap=0)
gen_net104_6_8 <- generate_communities_ml("pep", num.actors=104, 
                                      num.layers=6, num.communities=8,
                                      overlap=0)
```

Plot one of the generated networks

``` r
par(mfrow = c(1, 1))

lg <- layout_multiforce_ml(gen_net104_3_8$net, w_inter = 1, w_in = c(1, 0, 0),
                           gravity = c(1, 0, 0))
plot(gen_net104_3_8$net, layout = lg, grid = c(1, 3), vertex.labels = "")
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Algorithms comparison for every network

## Tuning parameters for SC-ML

``` r
tuning_k(krack_net, 0.5)
```

![](community_algorithms_comparison_files/figure-gfm/tuning%20parameters%20for%20SC-ML-1.png)<!-- -->

``` r
tuning_alpha(krack_net, 4)
```

![](community_algorithms_comparison_files/figure-gfm/tuning%20parameters%20for%20SC-ML-2.png)<!-- -->

``` r
tuning_k(aucs_net, 0.5)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
tuning_alpha(aucs_net, 6)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

## List of algorithms for real network data sets

``` r
krack_alg <- multilayer_algorithms(krack_net, 4, 0.5, p.min.layers = 3, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 127.1944 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
aucs_alg <- multilayer_algorithms(aucs_net, 6, 0.5, p.min.layers = 5, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-10.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-11.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 491.5809 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-8-12.png)<!-- -->

## List of algorithms for artificially generated network data sets

``` r
g33_3_3_alg <- multilayer_algorithms(gen_net33_3_3$net, 3, 0.5, p.min.layers = 3, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 185.6002 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->

``` r
g40_3_8_alg <- multilayer_algorithms(gen_net40_3_8$net, 8, 0.5, p.min.layers = 3, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 173.7522 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-10.png)<!-- -->

``` r
g99_3_3_alg <- multilayer_algorithms(gen_net99_3_3$net, 3, 0.5, p.min.layers = 3, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-11.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-12.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-13.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-14.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 1738.401 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-15.png)<!-- -->

``` r
g104_3_8_alg <- multilayer_algorithms(gen_net104_3_8$net, 8, 0.5, p.min.layers = 3, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-16.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-17.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-18.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-19.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 854.4442 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-20.png)<!-- -->

``` r
g104_6_8_alg <- multilayer_algorithms(gen_net104_6_8$net, 8, 0.5, p.min.layers = 6, p.min.actors = 2)
```

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-21.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-22.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-23.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-24.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-25.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-26.png)<!-- -->![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-27.png)<!-- -->

    ## 
    ## 
    ## Optimization of all partitions completed
    ## 1 solution(s) with minimal error = 1676.489 found.

![](community_algorithms_comparison_files/figure-gfm/unnamed-chunk-9-28.png)<!-- -->
## Comparison of algorithms for real network data-sets w/out
grount-truth labels

``` r
algorithms_comparison(krack_net, krack_alg[[1]], krack_alg[[2]])
```

    ## [1] "The algorithms comparison for the network: krack_net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 2                    10.50             1.00
    ## flatten_unweighted               2                    10.50             1.00
    ## abacus                           5                     2.80             0.67
    ## clique p.                        4                     4.00             0.76
    ## louvain                          3                     7.00             1.00
    ## infomap                          1                    21.00             1.00
    ## mdlp                             1                    21.00             1.00
    ## SCML                             4                     5.25             1.00
    ## Blockm                           4                     5.25             1.00
    ##                    modularity exec_time
    ## flatten_weighted         0.30      0.01
    ## flatten_unweighted       0.31      0.01
    ## abacus                   0.19      0.02
    ## clique p.                0.24      0.00
    ## louvain                  0.33      0.02
    ## infomap                  0.25      0.00
    ## mdlp                     0.25      0.01
    ## SCML                     0.32      0.17
    ## Blockm                   0.19     11.08

``` r
algorithms_comparison(aucs_net, aucs_alg[[1]], aucs_alg[[2]])
```

    ## [1] "The algorithms comparison for the network: aucs_net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 5                    12.20             1.00
    ## flatten_unweighted               5                    12.20             1.00
    ## abacus                           3                     2.33             0.11
    ## clique p.                        3                     2.00             0.10
    ## louvain                          5                    12.20             1.00
    ## infomap                          5                    12.20             1.00
    ## mdlp                             5                    12.20             1.00
    ## SCML                             6                    10.17             1.00
    ## Blockm                           6                    10.17             1.00
    ##                    modularity exec_time
    ## flatten_weighted         0.52      0.02
    ## flatten_unweighted       0.52      0.02
    ## abacus                   0.07      0.05
    ## clique p.                0.06      0.01
    ## louvain                  0.52      0.07
    ## infomap                  0.51      0.00
    ## mdlp                     0.45      0.03
    ## SCML                     0.51      0.41
    ## Blockm                   0.39     26.02

## Comparison of algorithms with grount-truth labels

``` r
algorithms_comparison(krack_net, krack_alg[[1]], krack_alg[[2]], 
                      true_labels=TRUE, labels='nodeDepartment')
```

    ## [1] "The algorithms comparison for the network: krack_net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 2                    10.50             1.00
    ## flatten_unweighted               2                    10.50             1.00
    ## abacus                           5                     2.80             0.67
    ## clique p.                        4                     4.00             0.76
    ## louvain                          3                     7.00             1.00
    ## infomap                          1                    21.00             1.00
    ## mdlp                             1                    21.00             1.00
    ## SCML                             4                     5.25             1.00
    ## Blockm                           4                     5.25             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.30 0.20     0.38      0.01
    ## flatten_unweighted       0.31 0.40     0.43      0.01
    ## abacus                   0.19   NA       NA      0.02
    ## clique p.                0.24   NA       NA      0.00
    ## louvain                  0.33 0.54     0.67      0.02
    ## infomap                  0.25 0.00     0.38      0.00
    ## mdlp                     0.25 0.00     0.38      0.01
    ## SCML                     0.32 0.86     0.90      0.17
    ## Blockm                   0.19 0.34     0.52     11.08

``` r
algorithms_comparison(aucs_net, aucs_alg[[1]], aucs_alg[[2]], 
                      true_labels=TRUE, labels='group', flag_NA=TRUE)
```

    ## [1] "The algorithms comparison for the network: aucs_net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 5                    12.20             1.00
    ## flatten_unweighted               5                    12.20             1.00
    ## abacus                           3                     2.33             0.11
    ## clique p.                        3                     2.00             0.10
    ## louvain                          5                    12.20             1.00
    ## infomap                          5                    12.20             1.00
    ## mdlp                             5                    12.20             1.00
    ## SCML                             6                    10.17             1.00
    ## Blockm                           6                    10.17             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.52 0.87     0.78      0.02
    ## flatten_unweighted       0.52 0.87     0.78      0.02
    ## abacus                   0.07   NA       NA      0.05
    ## clique p.                0.06   NA       NA      0.01
    ## louvain                  0.52 0.87     0.78      0.07
    ## infomap                  0.51 0.81     0.65      0.00
    ## mdlp                     0.45 0.70     0.64      0.03
    ## SCML                     0.51 0.82     0.81      0.41
    ## Blockm                   0.39 0.54     0.53     26.02

## Comparison of algorithms for generated network with grount-truth labels

``` r
algorithms_comparison(gen_net33_3_3$net, g33_3_3_alg[[1]], g33_3_3_alg[[2]], 
                      genBool=TRUE, gen_truel = gen_net33_3_3$com)
```

    ## [1] "The algorithms comparison for the network: gen_net33_3_3$net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 3                    11.00             1.00
    ## flatten_unweighted               3                    11.00             1.00
    ## abacus                           5                     6.20             0.94
    ## clique p.                        4                     2.50             0.30
    ## louvain                          3                    11.00             1.00
    ## infomap                          3                    11.00             1.00
    ## mdlp                             8                     4.12             1.00
    ## SCML                             3                    11.00             1.00
    ## Blockm                           3                    11.00             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.64 1.00     1.00      0.01
    ## flatten_unweighted       0.64 1.00     1.00      0.01
    ## abacus                   0.50   NA       NA      0.03
    ## clique p.                0.14   NA       NA      0.00
    ## louvain                  0.64 1.00     1.00      0.02
    ## infomap                  0.64 1.00     1.00      0.00
    ## mdlp                     0.43 0.57     0.91      0.01
    ## SCML                     0.51 0.75     0.85      0.17
    ## Blockm                   0.46 0.65     0.67     10.98

``` r
algorithms_comparison(gen_net40_3_8$net, g40_3_8_alg[[1]], g40_3_8_alg[[2]], 
                      genBool=TRUE, gen_truel = gen_net40_3_8$com)
```

    ## [1] "The algorithms comparison for the network: gen_net40_3_8$net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 8                     5.00             1.00
    ## flatten_unweighted               8                     5.00             1.00
    ## abacus                           9                     2.22             0.50
    ## clique p.                        5                     2.20             0.28
    ## louvain                          8                     5.00             1.00
    ## infomap                          8                     5.00             1.00
    ## mdlp                            15                     2.67             1.00
    ## SCML                             8                     5.00             1.00
    ## Blockm                           8                     5.00             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.58 1.00     1.00      0.01
    ## flatten_unweighted       0.58 1.00     1.00      0.01
    ## abacus                   0.22   NA       NA      0.02
    ## clique p.                0.13   NA       NA      0.01
    ## louvain                  0.58 1.00     1.00      0.02
    ## infomap                  0.58 1.00     1.00      0.00
    ## mdlp                     0.41 0.66     0.72      0.01
    ## SCML                     0.57 0.97     0.98      0.22
    ## Blockm                   0.38 0.52     0.50     13.07

``` r
algorithms_comparison(gen_net99_3_3$net, g99_3_3_alg[[1]], g99_3_3_alg[[2]], 
                      genBool=TRUE, gen_truel = gen_net99_3_3$com)
```

    ## [1] "The algorithms comparison for the network: gen_net99_3_3$net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 3                       33             1.00
    ## flatten_unweighted               3                       33             1.00
    ## abacus                           3                       33             1.00
    ## clique p.                        5                       17             0.86
    ## louvain                          3                       33             1.00
    ## infomap                          3                       33             1.00
    ## mdlp                             3                       33             1.00
    ## SCML                             3                       33             1.00
    ## Blockm                           3                       33             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.64 1.00     1.00      0.05
    ## flatten_unweighted       0.64 1.00     1.00      0.04
    ## abacus                   0.64 1.00     1.00      0.04
    ## clique p.                0.46   NA       NA      0.03
    ## louvain                  0.64 1.00     1.00      0.12
    ## infomap                  0.64 1.00     1.00      0.00
    ## mdlp                     0.64 1.00     1.00      0.11
    ## SCML                     0.64 1.00     1.00      0.47
    ## Blockm                   0.38 0.65     0.67     47.02

``` r
algorithms_comparison(gen_net104_3_8$net, g104_3_8_alg[[1]], g104_3_8_alg[[2]], 
                      genBool=TRUE, gen_truel = gen_net104_3_8$com)
```

    ## [1] "The algorithms comparison for the network: gen_net104_3_8$net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 8                    13.00             1.00
    ## flatten_unweighted               8                    13.00             1.00
    ## abacus                           8                    12.62             0.97
    ## clique p.                       20                     3.05             0.59
    ## louvain                          8                    13.00             1.00
    ## infomap                          8                    13.00             1.00
    ## mdlp                            20                     5.20             1.00
    ## SCML                             8                    13.00             1.00
    ## Blockm                           8                    13.00             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.71 1.00     1.00      0.03
    ## flatten_unweighted       0.71 1.00     1.00      0.03
    ## abacus                   0.68   NA       NA      0.04
    ## clique p.                0.25   NA       NA      0.03
    ## louvain                  0.71 1.00     1.00      0.08
    ## infomap                  0.71 1.00     1.00      0.00
    ## mdlp                     0.49 0.80     0.95      0.05
    ## SCML                     0.71 1.00     1.00      0.33
    ## Blockm                   0.49 0.74     0.66      1.26

``` r
algorithms_comparison(gen_net104_6_8$net, g104_6_8_alg[[1]], g104_6_8_alg[[2]], 
                      genBool=TRUE, gen_truel = gen_net104_6_8$com)
```

    ## [1] "The algorithms comparison for the network: gen_net104_6_8$net"
    ##                    num_communities avg_actors_per_community clustered_actors
    ## flatten_weighted                 8                    13.00             1.00
    ## flatten_unweighted               8                    13.00             1.00
    ## abacus                           8                    12.38             0.95
    ## clique p.                        3                     2.00             0.06
    ## louvain                          8                    13.00             1.00
    ## infomap                          8                    13.00             1.00
    ## mdlp                            19                     5.47             1.00
    ## SCML                             8                    13.00             1.00
    ## Blockm                           8                    13.00             1.00
    ##                    modularity  nmi accuracy exec_time
    ## flatten_weighted         0.78 1.00     1.00      0.03
    ## flatten_unweighted       0.78 1.00     1.00      0.03
    ## abacus                   0.74   NA       NA      0.06
    ## clique p.                0.03   NA       NA      0.04
    ## louvain                  0.78 1.00     1.00      0.17
    ## infomap                  0.78 1.00     1.00      0.01
    ## mdlp                     0.57 0.71     0.84      0.11
    ## SCML                     0.78 1.00     1.00      0.49
    ## Blockm                   0.69 0.85     0.77      2.53

## SCML testing on one layer

``` r
for (l in layers_ml(krack_net)){
  k_means_com <- onelayerSCML(krack_net, 4, l)
  k_means_com <- k_means_com[order(k_means_com$actor),]
  a <- actors_ml(krack_net, layers = l)$actor
  v <- vertices_ml(krack_net, layers = l)
  true_labels <- get_values_ml(krack_net, attribute='nodeDepartment', actors=a)
  net_df <- data.frame(actor=a, cid=true_labels[,1])
  net_df <- net_df %>%
    mutate_at(c('cid'), as.numeric) 
  net_df_l <- net_df[net_df$actor %in% k_means_com$actor, , drop = FALSE]
  net_df_l <- net_df_l[order(net_df_l$actor),]
  head(net_df_l)
  cat("NMI for layer", l, round(NMI(k_means_com$cid, net_df_l$cid, variant='sum'),2),
      "; purity score", round(purity_score(k_means_com$cid, net_df_l$cid), 2), 
      "\n")
} 
```

    ## NMI for layer advice 0.33 ; purity score 0.43 
    ## NMI for layer friendship 0.57 ; purity score 0.71 
    ## NMI for layer reportsto 0.96 ; purity score 0.95

``` r
for (l in layers_ml(aucs_net)){
  k_means_com <- onelayerSCML(aucs_net, 6, l)
  k_means_com <- k_means_com[order(k_means_com$actor),]
  a <- actors_ml(aucs_net)$actor
  true_labels <- get_values_ml(aucs_net, attribute="group",actors=a)
  net_df <- data.frame(actor=a, group=true_labels[,1])
  net_df_na <- net_df %>% subset(group=="NA")
  net_df <- net_df %>% subset(group!="NA") %>%
  #remain 1 last character -> -1
    mutate(cid = substr(group, nchar(group)-1+1, nchar(group))) %>% 
    mutate_at(c('cid'), as.numeric) %>% 
    select(c('actor','cid'))
  k_means_com <- k_means_com[!k_means_com$actor %in% net_df_na$actor, , drop = FALSE]
  
  net_df_l <- net_df[net_df$actor %in% k_means_com$actor, , drop = FALSE]
  net_df_l <- net_df_l[order(net_df_l$actor),]
  head(net_df_l)
  cat("NMI for layer", l, round(NMI(k_means_com$cid, net_df_l$cid, variant='sum'),2),
      "; purity score", round(purity_score(k_means_com$cid, net_df_l$cid),2), "\n")
} 
```

    ## NMI for layer lunch 0.83 ; purity score 0.76 
    ## NMI for layer facebook 0.7 ; purity score 0.74 
    ## NMI for layer coauthor 0.77 ; purity score 0.76 
    ## NMI for layer leisure 0.62 ; purity score 0.68 
    ## NMI for layer work 0.73 ; purity score 0.73
