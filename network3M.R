# function - network3M #
########################################################################
# File: network3M.R
# Aim : Draw a network diagram of the correlation between metabolites,microbes and phenotypes 
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen, Tao Sun, Dandan Liang, Mengci Li
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
########################################################################################
## Input:                                                                             ##
##    individualCorData --- A list included coefficient, p value and p.adjust value   ##
##                 of the correlation between metabolome data and microbiome data     ##  
##                 from individualCor.R                                               ##
##    rThreshold --- The correlation coefficient threshold to filter nodes            ##
##                 Default: 0.3)                                                      ##           
##    fdrThreshold --- The FDR threshold to filter nodes                              ##
##                Default: 0.05                                                       ##
##    centrality --- Set nodes size according to the centrality                       ##
##                Default: degree                                                     ##
##    phenoData --- The dataframe for confounder                                      ##
##                       Default: NA(rows:samples,columns:confounder)                 ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
##    network_edge_df ---  A list included coefficient, p value and p.adjust value    ##
##                 of the correlation between metabolome data and microbiome          ##
##                 data                                                               ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
network3M <- function(individualCorData,rThreshold = 0.3,fdrThreshold = 0.05,centrality = "degree",phenoData = NA)
{
  
  if(!file.exists("./results/Inter-Cor/pairwise")){
    dir.create("./results/Inter-Cor/pairwise",recursive = T)
  }

  ############ Exist continuous phenotype variables ############
  if(names(individualCorData)[1] == "mic2metaCor"){
    mic_meta_cor <- individualCorData[["mic2metaCor"]][["r"]]
    mic_meta_p <- individualCorData[["mic2metaCor"]][["p"]]
    mic_meta_p.adjust <- individualCorData[["mic2metaCor"]][["p.adjust"]]
    
    pheno_meta_cor <- individualCorData[["pheno2metaCor"]][["r"]]
    pheno_meta_p <- individualCorData[["pheno2metaCor"]][["p"]]
    pheno_meta_p.adjust <- individualCorData[["pheno2metaCor"]][["p.adjust"]]
    
    pheno_mic_cor <- individualCorData[["pheno2micCor"]][["r"]]
    pheno_mic_p <- individualCorData[["pheno2micCor"]][["p"]]
    pheno_mic_p.adjust <- individualCorData[["pheno2micCor"]][["p.adjust"]]
    
    ### Generate edge table
    CorrDF <- function(cormat, pmat) {
      ut = matrix (TRUE, nrow = nrow (cormat), ncol = ncol (pmat))
      data.frame(
        from = rownames(cormat)[row(cormat)[ut]],
        to = colnames(cormat)[col(cormat)[ut]],
        r =(cormat)[ut],
        p.adjust = pmat[ut]
      )
    }
    
    mic_meta_cor_df <- CorrDF(mic_meta_cor,mic_meta_p.adjust)
    pheno_meta_cor_df <- CorrDF(pheno_meta_cor,pheno_meta_p.adjust)
    pheno_mic_cor_df <- CorrDF(pheno_mic_cor,pheno_mic_p.adjust)
    
    ### The relationship of r > rThreshold and p-value < fdrThreshold is taken as the edge 
    ### of the selected co-occurrence network
    mic_meta_cor_df <- mic_meta_cor_df[which(abs(mic_meta_cor_df$r) > rThreshold),]
    mic_meta_cor_df <- mic_meta_cor_df[which(mic_meta_cor_df$p.adjust < fdrThreshold),]
    pheno_meta_cor_df <- pheno_meta_cor_df[which(abs(pheno_meta_cor_df$r) > rThreshold),]
    pheno_meta_cor_df <- pheno_meta_cor_df[which(pheno_meta_cor_df$p.adjust < fdrThreshold),]
    pheno_mic_cor_df <- pheno_mic_cor_df[which(abs(pheno_mic_cor_df$r) > rThreshold),]
    pheno_mic_cor_df <- pheno_mic_cor_df[which(pheno_mic_cor_df$p.adjust < fdrThreshold),]
    edge_df <- rbind(mic_meta_cor_df,pheno_meta_cor_df,pheno_mic_cor_df)
    if(nrow(edge_df) == 0 ){
      stop("There are no edges that satisfy these thresholds!")
    }
    
    
    ### Generate nodes attribute table 
    ### The two edge columns merge into nodes
    nodeattrib <- data.frame(node = union(edge_df$from,edge_df$to))

    
    ### Categorize nodes
    nodeattrib$group <- 0
    for (i in as.character(nodeattrib$node))
    {
      if (i %in% rownames(mic_meta_cor) == TRUE){
        nodeattrib[nodeattrib$node == i,"group"] <- "Microbe"
      }else if(i %in% colnames(mic_meta_cor) == TRUE){ 
        nodeattrib[nodeattrib$node == i,"group"] <- "Metabolite"
      }else{
        nodeattrib[nodeattrib$node == i,"group"] <- "Phenotype"
      }
    }
    rownames(nodeattrib) <- nodeattrib$node
    
    ### Draw a co-occurrence network diagram with igraph
    co_net <- graph_from_data_frame(edge_df,direct = F, vertices = nodeattrib) 
    
    ### The above coloring information is mapped into the nodes attribute table
    net_nodes <- rownames(nodeattrib[nodeattrib$group %in% co_net,])
    
    ### Add a Hexagon shape 
    Hexagon <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- "#6300ea"
      }
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      symbols(x = coords[,1], y = coords[,2], bg = vertex.color, fg = vertex.color,
              stars = cbind(vertex.size, vertex.size, vertex.size,vertex.size,vertex.size,vertex.size),
              add = TRUE, inches = FALSE)
    }
    add_shape("Hexagon",plot = Hexagon)
    
    triangle <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- "#6300ea"
      }
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      symbols(x = coords[,1], y = coords[,2], bg = vertex.color, fg = vertex.color,
              stars = cbind(vertex.size, vertex.size, vertex.size),
              add = TRUE, inches = FALSE)
    }
    add_shape("triangle",plot = triangle)
    
    ### Set the node shape
    V(co_net)$shape <- V(co_net)$group
    V(co_net)$shape[V(co_net)$shape == "Microbe"] <- "circle"
    V(co_net)$shape[V(co_net)$shape == "Metabolite"] <- "triangle"
    V(co_net)$shape[V(co_net)$shape == "Phenotype"] <- "square"
    
    ### Color different categories of nodes
    nodes_cat <- c("Microbe","Metabolite","Phenotype")
    unique(V(co_net)$group)
    V(co_net)$color <- V(co_net)$group
    V(co_net)$color[!V(co_net)$color %in% nodes_cat] <- "gray30"
    V(co_net)$color[V(co_net)$color == "Microbe"] <- "#00b692"
    V(co_net)$color[V(co_net)$color == "Metabolite"] <- "#6300ea"
    V(co_net)$color[V(co_net)$color == "Phenotype"] <- "#007de9"
    V(co_net)$frame.color <- V(co_net)$color
    
    ### Set edges color -- Positive:red; Negative: blue
    edge_col <- c()
    for (i in 1:nrow(edge_df) )
    {
      if (edge_df$r[i] > 0 ){
        edge_col <- append(edge_col,alpha("#d7191c", alpha = 1-edge_df$p.adjust[i]))
      }else{
        edge_col <- append(edge_col,alpha("#4BACC6", alpha = 1-edge_df$p.adjust[i]))
      }
    }
    co_net <- co_net %>% set_edge_attr("color", value = edge_col)
    
    ### Set edges size according to the value of r
    edge_size <- abs(edge_df$r * edge_df$r * edge_df$r * 10)
    co_net <- co_net %>% set_edge_attr("weight", value = edge_size)
    
    # ### Set edges transparency according to the value of p.adjust
    # edge_size <- edge_df$r * edge_df$r * 5
    # co_net <- co_net %>% set_edge_attr("weight", value = edge_size)
    
    ### Set nodes size according to the degree centrality,betweenness centrality or closeness centrality
    if(centrality == "degree"){
      nodes_size = log2(centr_degree(co_net)$res + 2) * 3
    }else if(centrality == "betweenness"){
      nodes_size = centr_betw(co_net, directed = FALSE)$res^(1/3) * 3
    }else if(centrality == "closeness"){
      nodes_size = centr_clo(co_net, mode = "all")$res * 10
    }
    
    ### The nodes in the network are clustered using fast greedy
    cfg <- cluster_fast_greedy(as.undirected(co_net))
    cluster_modules <- sort(table(membership(cfg)),decr=T)
    
    ### If there are more than 5 modules, take the first 5 modules
    if(length(cluster_modules) > 5 ){
      cluster_modules <- cluster_modules[1:5]
    }
    sm_plot <- cluster_modules
    names(sm_plot) <- as.factor(1:length(cluster_modules))
    
    ### Vectorize the above modules
    cluster_modules_points <- membership(cfg)[membership(cfg) %in% names(cluster_modules)]
    cluster_points <- NULL
    for(i in cluster_modules_points){
      tx <- which(names(cluster_modules )== i)
      cluster_points <- c(cluster_points, tx)
    }
    names(cluster_points) <- names(cluster_modules_points)
    
    ### Categorize microbe, metabolite and phenotype nodes into modules
    microbe_nodes <- rownames(nodeattrib[nodeattrib$group=="Microbe",])
    metabolite_nodes <- rownames(nodeattrib[nodeattrib$group=="Metabolite",])
    phenotype_nodes <- rownames(nodeattrib[nodeattrib$group=="Phenotype",])
    cs_nodes_all <- c(microbe_nodes,metabolite_nodes,phenotype_nodes)
    
    mods_list_cs <- list()
    edge_list_cs <- list()
    for (i in names(cluster_modules)){
      x1 <- names(membership(cfg)[membership(cfg) == i])
      x2 <- x1[x1 %in% cs_nodes_all]
      mods_list_cs[[i]] <- as.numeric(V(co_net)[x2])
      edge_list_cs[[i]] <- V(co_net)[x2]
    }

    
    ### pdf.1
    pdf(paste0("./results/Inter-Cor/pairwise/network3M_label.pdf"),width = 7,height = 5)
    ### modules colors
    # mark_cols <- c()
    # for (i in 1:length(cluster_modules)) {
    #   mark_cols_name <- paste0("antiquewhite",i)
    #   mark_cols <- append(mark_cols,mark_cols_name)
    # }
    # 
    mark_cols <- colorRampPalette(c("#d4f1ff","#d9ffd4"))(length(cluster_modules))
    
    ### label size
    label_size = log10(nodes_size)/2 + 0.5
    
    ### modules mark groups is opposite to mods_list_cs
    mark_groups <- list()
    edge_groups <- list()
    for (i in 1:length(cluster_modules)) {
      mark_groups[[i]] <- mods_list_cs[[length(cluster_modules) + 1 - i]]
      edge_groups[[i]] <- edge_list_cs[[length(cluster_modules) + 1 - i]]
    }
    
    ### Set the number of random seeds, and start from the same number of 
    ### random seeds in subsequent drawing to ensure that the shapes of the 
    ### drawing before and after are corresponding
    set.seed(123)
    coords <- layout_(co_net,with_fr(niter=9999, grid="nogrid"))
    plot(co_net,vertex.label = V(co_net)$name,vertex.label.cex = label_size,vertex.label.color = "black",
         layout=coords,vertex.size = nodes_size,
         edge.lty = 1,edge.curved=TRUE,
         margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight),
         mark.groups = mark_groups,
         mark.col = mark_cols, mark.border = mark_cols
    )
    modules_leg <- c()
    for (i in 1:length(cluster_modules)) {
      modules_leg <- append(modules_leg,paste0("Cluster",i))
    }
    legend("topright",legend = modules_leg,col = mark_cols,
           bty = "n",fill = mark_cols,border = mark_cols)
    leg_cols <- c("#00b692","#6300ea","#007de9")
    legend("bottomright",legend = c("Microbe","Metabolite", "Phenotype"),col = leg_cols,
           bty = "n",pch = c(16,17,15))
    dev.off()
    
    ### pdf.2
    pdf(paste0("./results/Inter-Cor/pairwise/network3M.pdf"),width = 7,height = 5)
    set.seed(123)
    coords <- layout_(co_net,with_fr(niter=9999, grid="nogrid"))
    plot(co_net,layout=coords,vertex.label = NA,
         vertex.size = nodes_size,
         edge.lty = 1,edge.curved=TRUE,
         margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight),
         mark.groups = mark_groups,
         mark.col = mark_cols, mark.border = mark_cols
    )
    legend("topright",legend = modules_leg,col = mark_cols,
           bty = "n",fill = mark_cols,border = mark_cols)
    leg_cols <- c("#00b692","#6300ea","#007de9")
    legend("bottomright",legend = c("Microbe","Metabolite", "Phenotype"),col = leg_cols,
           bty = "n",pch = c(16,17,15))
    dev.off()
    
    
    ### Topological index table of the graphics level
    ### 1 - Average nearest neighbor degree
    Average_nearest_neighbor_degree <- mean(knn(co_net)$knn)
    ### 2 - Average path length
    Average_path_length <- mean_distance(co_net, directed = F)
    ### 3 - Betweenness centrality
    Betweenness_centrality <- centr_betw(co_net,directed = FALSE)$centralization
    ### 4 - Closeness centrality
    Closeness_centrality <- centr_clo(co_net,mode = "all")$centralization
    ### 5 - Degree assortativity
    Degree_assortativity <- assortativity_degree(co_net, directed = F)
    ### 6 - Degree centralization
    Degree_centralization <- centr_degree(co_net)$centralization
    ### 7 - Density
    Density <- edge_density(co_net)
    ### 8 - Transitivity
    Transitivity <- transitivity(co_net, type="global")
    ### 9 - Number of vertice
    Number_of_vertice <- vcount(co_net)
    ### 10 - Number of edge
    Number_of_edge <- ecount(co_net)
    ### 11 - Modularity
    Modularity <- modularity(co_net, membership(cfg))
    ### 12 - Diameter
    Diameter <- diameter(co_net, unconnected=FALSE)
    ### 13 - Cluster counts
    Cluster_counts <- length(cluster_modules)
    
    topo_table <- cbind(Average_nearest_neighbor_degree,Average_path_length,
                        Betweenness_centrality,Closeness_centrality,
                        Degree_assortativity,Degree_centralization,
                        Density,Transitivity,Number_of_vertice,
                        Number_of_edge,Modularity,Diameter,Cluster_counts)
    topo_table <- as.data.frame(topo_table)
    write.table(topo_table, file = "./results/Inter-Cor/pairwise/Topological_index_table.txt", sep = "\t", row.names = FALSE, quote = F)
    
    ### Zi-Pi plot
    ### Zi:Within-module connectivity; Pi:Among-module connectivity
    
    ### Generate subgraph
    subgraph_list <- list()
    for (i in 1:length(cluster_modules)) {
      subgraph_list[[i]] <- induced_subgraph(co_net, edge_groups[[i]])
    }
    
    zipiSores_table <- matrix(data = NA, nrow = 1, ncol = 3)
    colnames(zipiSores_table) <- c("node","Zi","Pi")
    ### n:module number
    for (n in 1:length(cluster_modules)) {
      ### i: nodes number in subgraph 
      for(i in 1:vcount(subgraph_list[[n]])){
        Node_name <- names(igraph::degree(subgraph_list[[n]])[i])
        ### Kin is the number of edges between node i and other nodes in module n
        Kin <- igraph::degree(subgraph_list[[n]])[i]
        ### Kan is the average of the degrees of all nodes in module n
        Kan <- mean(igraph::degree(subgraph_list[[n]]))
        ### Sdn is the standard deviation of the degree of all nodes in module n
        Sdn <- sd(igraph::degree(subgraph_list[[n]]))
        ### Ki is the number of edges between node i and other nodes in the whole network
        Ki <- as.numeric(igraph::degree(co_net)[Node_name]) 
        ### Filter the edges that contain node i
        temp_row_names <- row.names(edge_df[rbind(which(Node_name == edge_df$from),which(Node_name == edge_df$to)),])
        nodei_edge_df <- edge_df[temp_row_names,]
        cum_sum <- 0
        for (m in 1:length(cluster_modules)) {
          ### Kic is the number of edges between node i and other nodes in each module
          Kic <- 0
          for (j in 1:length(V(subgraph_list[[m]])$name)) {
            ### V(subgraph_list[[m]])$name[j]) can not be Node_name
            ### because Kic is the number of edges between node i and *other* nodes in each module
            if(V(subgraph_list[[m]])$name[j] == Node_name){ 
              next
            }
            if((V(subgraph_list[[m]])$name[j] %in% nodei_edge_df$from)||(V(subgraph_list[[m]])$name[j] %in% nodei_edge_df$to)){
              Kic <- Kic + 1
            }
            # if(j >= length(V(subgraph_list[[m]])$name)){
            #   break
            # }
          }
          cum_sum <- cum_sum + (Kic/Ki)^2
        }
        ### Cumulative Sums
        zi <- (Kin - Kan) / Sdn
        pi <- 1 - cum_sum
        temp_c <- c(Node_name,zi,pi)
        zipiSores_table <- rbind(zipiSores_table,temp_c)
        
      }
    }
    zipiSores_table <- as.data.frame(zipiSores_table[-1,])
    rownames(zipiSores_table) <- zipiSores_table[,1]
    
    ## plot
    ### Filter NaN value
    n = grep("NaN",zipiSores_table$Zi)
    m = grep("NaN",zipiSores_table$Pi)
    All <- c(n,m)
    if(length(All) != 0){
      zi_pi <- zipiSores_table[-All,]
    }else{
      zi_pi <- zipiSores_table
    }
    zi_pi <- merge(zi_pi,nodeattrib,all.x=TRUE,by='node')
    zi_pi$Zi <- as.numeric(zi_pi$Zi)
    zi_pi$Pi <- as.numeric(zi_pi$Pi)
    if(!is.na(max(zi_pi$Zi)) & max(zi_pi$Zi) < 2.5){
      max_y <- 4
    }else if(!is.na(max(zi_pi$Zi)) & max(zi_pi$Zi) > 2.5){
      max_y <- max(zi_pi$Zi) + 2
    }else{
      max_y <- 4
    }
    
    colnames(zi_pi)[4] <- "Node_type"
    
    zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Peripherals'
    zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Connectors'
    zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Module hubs'
    zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Network hubs'
    
    label_points_row <- c()
    label_points_row <- append(label_points_row,which(zi_pi$type == 'Connectors'))
    label_points_row <- append(label_points_row,which(zi_pi$type == 'Module hubs'))
    label_points_row <- append(label_points_row,which(zi_pi$type == 'Network hubs'))
    show_label_points <- zi_pi[label_points_row,]
    
    p <- ggplot(data = zi_pi, aes(x = Pi, y = Zi, color = Node_type)) +
      geom_point(size = 2) +
      scale_color_manual(values = c('#00b692','#6300ea','#007de9'),
                         limits = c('Microbe','Metabolite','Phenotype'))+
      theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
            axis.text = element_text(size = 10,face = "bold"),
            panel.background = element_blank(), legend.key = element_blank(),
            panel.border = element_blank()) +
      labs(x = 'Among-module connectivities(Pi)', y = 'Within-module connectivities(Zi)', size = 6) +
      geom_vline(xintercept = 0.62,linetype = "dashed") +
      geom_hline(yintercept = 2.5 ,linetype = "dashed") + 
      scale_x_continuous(limits = c(0,1),expand = c(0,0)) +
      scale_y_continuous(limits = c(min(zi_pi$Zi),max_y)) +
      # geom_text(aes(label = node)) +
      annotate("text",x = c(0.25, 0.85, 0.25,0.85), y = c(min(zi_pi$Zi)+0.2, min(zi_pi$Zi)+0.2, max_y-0.2,max_y-0.2), 
               label = c('Peripherals', 'Connectors', 'Module hubs','Network hubs'), size = 4, colour = 'black') +
      geom_text_repel(data = show_label_points,aes(x = Pi, y = Zi,label = node),size = 3,show.legend = FALSE)
    ggsave("./results/Inter-Cor/pairwise/zipi_plot.pdf",width = 7,height = 5)
    
    return(edge_df)
  }
  
  ############ No phenotypes or exist categorical phenotype variables ############
  else{
    network_func <- function(indiCor,rThres,fdrThres,centr,pheno){
      mic_meta_cor <- indiCor[[paste0("r_",pheno)]]
      mic_meta_p <- indiCor[[paste0("p_",pheno)]]
      mic_meta_p.adjust <- indiCor[[paste0("p.adjust_",pheno)]]
      ### Generate edge table
      CorrDF <- function(cormat, pmat) {
        ut = matrix (TRUE, nrow = nrow (cormat), ncol = ncol (pmat))
        data.frame(
          from = rownames(cormat)[row(cormat)[ut]],
          to = colnames(cormat)[col(cormat)[ut]],
          r =(cormat)[ut],
          p.adjust = pmat[ut]
        )
      }
      
      mic_meta_cor_df <- CorrDF(mic_meta_cor,mic_meta_p.adjust)
      
      ### The relationship of r > rThres and p-value < fdrThres is taken as the edge 
      ### of the selected co-occurrence network
      mic_meta_cor_df <- mic_meta_cor_df[which(abs(mic_meta_cor_df$r) > rThres),]
      mic_meta_cor_df <- mic_meta_cor_df[which(mic_meta_cor_df$p.adjust < fdrThres),]
      edge_df <- mic_meta_cor_df
      
      
      ### Generate nodes attribute table 
      ### The two edge columns merge into nodes
      nodeattrib <- data.frame(node = union(edge_df$from,edge_df$to))
      
      if(nrow(nodeattrib) == 0){
        stop("There are no nodes that satisfy these thresholds!")
      }
      
      ### Categorize nodes
      nodeattrib$group <- 0
      for (i in as.character(nodeattrib$node)){
        if (i %in% rownames(mic_meta_cor) == TRUE){
          nodeattrib[nodeattrib$node == i,"group"] <- "Microbe"
        }else if(i %in% colnames(mic_meta_cor) == TRUE){ 
          nodeattrib[nodeattrib$node == i,"group"] <- "Metabolite"
        }
      }
      rownames(nodeattrib) <- nodeattrib$node
      
      ### Draw a co-occurrence network diagram with igraph
      co_net <- graph_from_data_frame(edge_df,direct = F, vertices = nodeattrib) 
      
      ### The above coloring information is mapped into the nodes attribute table
      net_nodes <- rownames(nodeattrib[nodeattrib$group %in% co_net,])
      
      ### Add a Hexagon shape 
      Hexagon <- function(coords, v=NULL, params) {
        vertex.color <- params("vertex", "color")
        if (length(vertex.color) != 1 && !is.null(v)) {
          vertex.color <- "#6300ea"
        }
        vertex.size <- 1/200 * params("vertex", "size")
        if (length(vertex.size) != 1 && !is.null(v)) {
          vertex.size <- vertex.size[v]
        }
        symbols(x = coords[,1], y = coords[,2], bg = vertex.color, fg = vertex.color,
                stars = cbind(vertex.size, vertex.size, vertex.size,vertex.size,vertex.size,vertex.size),
                add = TRUE, inches = FALSE)
      }
      add_shape("Hexagon",plot = Hexagon)
      
      triangle <- function(coords, v=NULL, params) {
        vertex.color <- params("vertex", "color")
        if (length(vertex.color) != 1 && !is.null(v)) {
          vertex.color <- "#6300ea"
        }
        vertex.size <- 1/200 * params("vertex", "size")
        if (length(vertex.size) != 1 && !is.null(v)) {
          vertex.size <- vertex.size[v]
        }
        symbols(x = coords[,1], y = coords[,2], bg = vertex.color, fg = vertex.color,
                stars = cbind(vertex.size, vertex.size, vertex.size),
                add = TRUE, inches = FALSE)
      }
      add_shape("triangle",plot = triangle)
      
      ### Set the node shape
      V(co_net)$shape <- V(co_net)$group
      V(co_net)$shape[V(co_net)$shape == "Microbe"] <- "circle"
      V(co_net)$shape[V(co_net)$shape == "Metabolite"] <- "triangle"
      
      ### Color different categories of nodes
      nodes_cat <- c("Microbe","Metabolite")
      unique(V(co_net)$group)
      V(co_net)$color <- V(co_net)$group
      V(co_net)$color[!V(co_net)$color %in% nodes_cat] <- "gray30"
      V(co_net)$color[V(co_net)$color == "Microbe"] <- "#00b692"
      V(co_net)$color[V(co_net)$color == "Metabolite"] <- "#6300ea"
      V(co_net)$frame.color <- V(co_net)$color
      
      ### Set edges color -- Positive:red; Negative: blue
      edge_col <- c()
      for (i in 1:nrow(edge_df) )
      {
        if (edge_df$r[i] > 0 ){
          edge_col <- append(edge_col,alpha("#d7191c", alpha = 1-edge_df$p.adjust[i]))
        }else{
          edge_col <- append(edge_col,alpha("#4BACC6", alpha = 1-edge_df$p.adjust[i]))
        }
      }
      co_net <- co_net %>% set_edge_attr("color", value = edge_col)
      
      ### Set edges size according to the value of r
      edge_size <- abs(edge_df$r * edge_df$r * edge_df$r * 10)
      co_net <- co_net %>% set_edge_attr("weight", value = edge_size)
      
      # ### Set edges transparency according to the value of p.adjust
      # edge_size <- edge_df$r * edge_df$r * 5
      # co_net <- co_net %>% set_edge_attr("weight", value = edge_size)
      
      ### Set nodes size according to the degree centrality,betweenness centrality or closeness centrality
      if(centr == "degree"){
        nodes_size = log2(centr_degree(co_net)$res + 2) * 3
      }else if(centr == "betweenness"){
        nodes_size = centr_betw(co_net, directed = FALSE)$res^(1/3) * 3
      }else if(centr == "closeness"){
        nodes_size = centr_clo(co_net, mode = "all")$res * 10
      }
      
      ### The nodes in the network are clustered using fast greedy
      cfg <- cluster_fast_greedy(as.undirected(co_net))
      cluster_modules <- sort(table(membership(cfg)),decr=T)
      
      ### If there are more than 5 modules, take the first 5 modules
      if(length(cluster_modules) > 5 ){
        cluster_modules <- cluster_modules[1:5]
      }
      sm_plot <- cluster_modules
      names(sm_plot) <- as.factor(1:length(cluster_modules))
      
      ### Vectorize the above modules
      cluster_modules_points <- membership(cfg)[membership(cfg) %in% names(cluster_modules)]
      cluster_points <- NULL
      for(i in cluster_modules_points){
        tx <- which(names(cluster_modules )== i)
        cluster_points <- c(cluster_points, tx)
      }
      names(cluster_points) <- names(cluster_modules_points)
      
      ### Categorize microbe and metabolite nodes into modules
      microbe_nodes <- rownames(nodeattrib[nodeattrib$group=="Microbe",])
      metabolite_nodes <- rownames(nodeattrib[nodeattrib$group=="Metabolite",])
      cs_nodes_all <- c(microbe_nodes,metabolite_nodes)
      
      mods_list_cs <- list()
      edge_list_cs <- list()
      for (i in names(cluster_modules)){
        x1 <- names(membership(cfg)[membership(cfg) == i])
        x2 <- x1[x1 %in% cs_nodes_all]
        mods_list_cs[[i]] <- as.numeric(V(co_net)[x2])
        edge_list_cs[[i]] <- V(co_net)[x2]
      }
      
      ### pdf.1
      pdf(paste0("./results/Inter-Cor/pairwise/network3M_label_",pheno,".pdf"),width = 7,height = 5)
      ### modules colors
      # mark_cols <- c()
      # for (i in 1:length(cluster_modules)) {
      #   mark_cols_name <- paste0("antiquewhite",i)
      #   mark_cols <- append(mark_cols,mark_cols_name)
      # }
      # 
      mark_cols <- colorRampPalette(c("#d4f1ff","#d9ffd4"))(length(cluster_modules))
      
      ### label size
      label_size = log10(nodes_size)/2 + 0.5
      
      ### modules mark groups is opposite to mods_list_cs and modules edge groups is opposite to edge_list_cs
      mark_groups <- list()
      edge_groups <- list()
      for (i in 1:length(cluster_modules)) {
        mark_groups[[i]] <- mods_list_cs[[length(cluster_modules) + 1 - i]]
        edge_groups[[i]] <- edge_list_cs[[length(cluster_modules) + 1 - i]]
      }
      
      ### Set the number of random seeds, and start from the same number of 
      ### random seeds in subsequent drawing to ensure that the shapes of the 
      ### drawing before and after are corresponding
      set.seed(8051)
      coords <- layout_(co_net,with_fr(niter=9999, grid="nogrid"))
      plot(co_net,vertex.label = V(co_net)$name,vertex.label.cex = label_size,vertex.label.color = "black",
           layout=coords,vertex.size = nodes_size,
           edge.lty = 1,edge.curved=TRUE,
           margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight),
           mark.groups = mark_groups,
           mark.col = mark_cols, mark.border = mark_cols
      )
      modules_leg <- c()
      for (i in 1:length(cluster_modules)) {
        modules_leg <- append(modules_leg,paste0("Cluster",i))
      }
      legend("topright",legend = modules_leg,col = mark_cols,
             bty = "n",fill = mark_cols,border = mark_cols)
      leg_cols <- c("#00b692","#6300ea")
      legend("bottomright",legend = c("Microbe","Metabolite"),col = leg_cols,
             bty = "n",pch = c(16,17))
      dev.off()
      
      ### pdf.2
      pdf(paste0("./results/Inter-Cor/pairwise/network3M_",pheno,".pdf"),width = 7,height = 5)
      set.seed(8051)
      coords <- layout_(co_net,with_fr(niter=9999, grid="nogrid"))
      plot(co_net,layout=coords,vertex.label = NA,
           vertex.size = nodes_size,
           edge.lty = 1,edge.curved=TRUE,
           margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight),
           mark.groups = mark_groups,
           mark.col = mark_cols, mark.border = mark_cols
      )
      legend("topright",legend = modules_leg,col = mark_cols,
             bty = "n",fill = mark_cols,border = mark_cols)
      leg_cols <- c("#00b692","#6300ea")
      legend("bottomright",legend = c("Microbe","Metabolite"),col = leg_cols,
             bty = "n",pch = c(16,17))
      dev.off()
      
      ### Topological index table of the graphics level
      ### 1 - Average nearest neighbor degree
      Average_nearest_neighbor_degree <- mean(knn(co_net)$knn)
      ### 2 - Average path length
      Average_path_length <- mean_distance(co_net, directed = F)
      ### 3 - Betweenness centrality
      Betweenness_centr <- centr_betw(co_net,directed = FALSE)$centralization
      ### 4 - Closeness centrality
      Closeness_centr <- centr_clo(co_net,mode = "all")$centralization
      ### 5 - Degree assortativity
      Degree_assortativity <- assortativity_degree(co_net, directed = F)
      ### 6 - Degree centralization
      Degree_centralization <- centr_degree(co_net)$centralization
      ### 7 - Density
      Density <- edge_density(co_net)
      ### 8 - Transitivity
      Transitivity <- transitivity(co_net, type="global")
      ### 9 - Number of vertice
      Number_of_vertice <- vcount(co_net)
      ### 10 - Number of edge
      Number_of_edge <- ecount(co_net)
      ### 11 - Modularity
      Modularity <- modularity(co_net, membership(cfg))
      ### 12 - Diameter
      Diameter <- diameter(co_net, unconnected=FALSE)
      ### 13 - Cluster counts
      Cluster_counts <- length(cluster_modules)
      
      topo_table <- cbind(Average_nearest_neighbor_degree,Average_path_length,
                          Betweenness_centr,Closeness_centr,
                          Degree_assortativity,Degree_centralization,
                          Density,Transitivity,Number_of_vertice,
                          Number_of_edge,Modularity,Diameter,Cluster_counts)
      topo_table <- as.data.frame(topo_table)
      write.table(topo_table, file = paste0("./results/Inter-Cor/pairwise/Topological_index_table_",pheno,".txt"), sep = "\t", row.names = FALSE, quote = F)
      
      
      ### Zi-Pi plot
      ### Zi:Within-module connectivity; Pi:Among-module connectivity
      
      ### Generate subgraph
      subgraph_list <- list()
      for (i in 1:length(cluster_modules)) {
        subgraph_list[[i]] <- induced_subgraph(co_net, edge_groups[[i]])
      }
      
      zipiSores_table <- matrix(data = NA, nrow = 1, ncol = 3)
      colnames(zipiSores_table) <- c("node","Zi","Pi")
      ### n:module number
      for (n in 1:length(cluster_modules)) {
        ### i: nodes number in subgraph 
        for(i in 1:vcount(subgraph_list[[n]])){
          Node_name <- names(igraph::degree(subgraph_list[[n]])[i])
          ### Kin is the number of edges between node i and other nodes in module n
          Kin <- igraph::degree(subgraph_list[[n]])[i]
          ### Kan is the average of the degrees of all nodes in module n
          Kan <- mean(igraph::degree(subgraph_list[[n]]))
          ### Sdn is the standard deviation of the degree of all nodes in module n
          Sdn <- sd(igraph::degree(subgraph_list[[n]]))
          ### Ki is the number of edges between node i and other nodes in the whole network
          Ki <- as.numeric(igraph::degree(co_net)[Node_name]) 
          ### Filter the edges that contain node i
          temp_row_names <- row.names(edge_df[rbind(which(Node_name == edge_df$from),which(Node_name == edge_df$to)),])
          nodei_edge_df <- edge_df[temp_row_names,]
          cum_sum <- 0
          for (m in 1:length(cluster_modules)) {
            ### Kic is the number of edges between node i and other nodes in each module
            Kic <- 0
            for (j in 1:length(V(subgraph_list[[m]])$name)) {
              ### V(subgraph_list[[m]])$name[j]) can not be Node_name
              ### because Kic is the number of edges between node i and *other* nodes in each module
              if(V(subgraph_list[[m]])$name[j] == Node_name){ 
                next
              }
              if((V(subgraph_list[[m]])$name[j] %in% nodei_edge_df$from)||(V(subgraph_list[[m]])$name[j] %in% nodei_edge_df$to)){
                Kic <- Kic + 1
              }
              # if(j >= length(V(subgraph_list[[m]])$name)){
              #   break
              # }
            }
            cum_sum <- cum_sum + (Kic/Ki)^2
          }
          ### Cumulative Sums
          zi <- (Kin - Kan) / Sdn
          pi <- 1 - cum_sum
          temp_c <- c(Node_name,zi,pi)
          zipiSores_table <- rbind(zipiSores_table,temp_c)
          
        }
      }
      zipiSores_table <- as.data.frame(zipiSores_table[-1,])
      rownames(zipiSores_table) <- zipiSores_table[,1]
      
      ## plot
      ### Filter NaN value
      n = grep("NaN",zipiSores_table$Zi)
      m = grep("NaN",zipiSores_table$Pi)
      All <- c(n,m)
      if(length(All) != 0){
        zi_pi <- zipiSores_table[-All,]
      }else{
        zi_pi <- zipiSores_table
      }
      zi_pi <- merge(zi_pi,nodeattrib,all.x=TRUE,by='node')
      zi_pi$Zi <- as.numeric(zi_pi$Zi)
      zi_pi$Pi <- as.numeric(zi_pi$Pi)
      if(!is.na(max(zi_pi$Zi)) & max(zi_pi$Zi) < 2.5){
        max_y <- 4
      }else{
        max_y <- max(zi_pi$Zi) + 2
      }
      
      colnames(zi_pi)[4] <- "Node_type"
      
      zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Peripherals'
      zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Connectors'
      zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Module hubs'
      zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Network hubs'
      
      label_points_row <- c()
      label_points_row <- append(label_points_row,which(zi_pi$type == 'Connectors'))
      label_points_row <- append(label_points_row,which(zi_pi$type == 'Module hubs'))
      label_points_row <- append(label_points_row,which(zi_pi$type == 'Network hubs'))
      show_label_points <- zi_pi[label_points_row,]
      
      p <- ggplot(data = zi_pi, aes(x = Pi, y = Zi,color = Node_type)) +
        geom_point(size = 2) +
        scale_color_manual(values = c('#00b692','#6300ea'),
                           limits = c('Microbe','Metabolite'))+
        theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
              axis.text = element_text(size = 10,face = "bold"),
              panel.background = element_blank(), legend.key = element_blank(),
              panel.border = element_blank()) +
        labs(x = 'Among-module connectivities(Pi)', y = 'Within-module connectivities(Zi)', size = 6) +
        geom_vline(xintercept = 0.62,linetype = "dashed") +
        geom_hline(yintercept = 2.5 ,linetype = "dashed") + 
        scale_x_continuous(limits = c(0,1),expand = c(0,0)) +
        scale_y_continuous(limits = c(min(zi_pi$Zi),max_y)) +
        # geom_text(aes(label = node)) +
        annotate("text",x = c(0.25, 0.85, 0.25,0.85), y = c(min(zi_pi$Zi)+0.2, min(zi_pi$Zi)+0.2, max_y-0.2,max_y-0.2), 
                 label = c('Peripherals', 'Connectors', 'Module hubs','Network hubs'), size = 4, colour = 'black')+
        geom_text_repel(data = show_label_points,aes(x = Pi, y = Zi,label = node),size = 3,show.legend = FALSE)
      ggsave(paste0("./results/Inter-Cor/pairwise/zipi_plot_",pheno,".pdf"),width = 7,height = 5)
      
      if(pheno == "All"){
        return(edge_df)
      }
    }
    # No phenotypes
    if(length(individualCorData)== 3){  # results of non-gramm
      names(individualCorData) <- c("r_All","p_All","p.adjust_All")
      network_func(indiCor = individualCorData,rThres = rThreshold,fdrThres = fdrThreshold,centr = centrality,pheno = "All")
    }else if(length(individualCorData)== 4){ # results of gramm 
      names(individualCorData) <- c("r_All","p_All","type_All","p.adjust_All")
      network_func(indiCor = individualCorData,rThres = rThreshold,fdrThres = fdrThreshold,centr = centrality,pheno = "All")
    }else{
      # Categorical phenotype variables
      for (i in unique(phenoData[,1])) {
        tempindividualCorData <- list()
        tempindividualCorData[[paste0("r_",i)]] <- individualCorData[[paste0("r_",i)]]
        tempindividualCorData[[paste0("p_",i)]] <- individualCorData[[paste0("p_",i)]]
        tempindividualCorData[[paste0("p.adjust_",i)]] <- individualCorData[[paste0("p.adjust_",i)]]
        network_func(indiCor = tempindividualCorData,rThres = rThreshold,fdrThres = fdrThreshold,centr = centrality,pheno = i)
      }
      # All samples
      network_edge_df <- network_func(indiCor = individualCorData,rThres = rThreshold,fdrThres = fdrThreshold,centr = centrality,pheno = "All")
    }  
  return(network_edge_df)
  }
}
