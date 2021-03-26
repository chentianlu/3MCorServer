# function - sankeyPlot #
########################################################################
# File: sankeyPlot.R
# Aim : Draw a sankey diagram with the selected correlation modules, the 
#       important metabolites and functions, and their affiliations
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
##    metaCluData --- A list included the classified metabolome modules results from  ##
##                   metabolitesCluster.R                                             ##
##    micCluData --- A list included the classified microbiome modules results from   ##
##                   microbesCluster.R                                                ##
##    cluCorData --- A list from clustersCor.R included coefficient, p value          ##
##                        and p.adjust value of the correlation between metabolome    ##
##                        modules and microbiome modules from clustersCor.R           ##
##    phenoData --- A dataframe of phenotype data                                     ##
##                Default: NA  (rows: samples,columns: phenotype)                     ##
##    phenoDataType --- Phenotype data type                                           ##
##                Default: continuous                                                 ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
##    Graph results in PDF format                                                     ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
sankeyPlot <- function(metaCluData,micCluData,cluCorData,phenoData = NA,phenoDataType = "continuous")
{
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  sankey_func <- function(metaClu,micClu,cluCor,phenotype){
    meta_cluster_data <- metaClu[[1]]
    mic_cluster_data <- micClu[[1]]
    cluster_cor_data <- cluCor[[1]]
    ### PCA,PCoA
    if(all(meta_cluster_data$cluster_name == "all")){
      ### Select metabolite clusters according to cluster_cor_data
      meta_clu_names <- colnames(cluster_cor_data)
      meta_filter_data <- as.data.frame(meta_cluster_data[,colnames(meta_cluster_data) %in% meta_clu_names]) 
      colnames(meta_filter_data) <- meta_clu_names
      rownames(meta_filter_data) <- rownames(meta_cluster_data)
      meta_links <- data.frame(source = 1,target = 1)
      for (i in colnames(meta_filter_data)) {
        temp <- as.data.frame(meta_filter_data[,i])
        colnames(temp) <- i
        rownames(temp) <- rownames(meta_filter_data)
        temp <- abs(temp)
        temp <- as.data.frame(apply(temp,2,sort,decreasing = T)) 
        ### Top 3 metabolites
        if(nrow(temp) > 3){
          temp_top3 <- as.data.frame(temp[1:3,])
          colnames(temp_top3) <- colnames(temp)
          rownames(temp_top3) <- rownames(temp)[1:3]
        }else{
          temp_top3 <- temp
        }
        temp_top3$modulte_name <- colnames(temp)[1]
        temp_meta_links <- data.frame(source = rownames(temp_top3),target = colnames(temp)[1])
        meta_links <- rbind(meta_links,temp_meta_links)
      }
      meta_links <- meta_links[-1,]
      meta_links$group <- "meta"
      
      ### meta_mic_links
      meta_mic_links <- data.frame(source = 1,target = 1)
      for (i in 1:ncol(cluster_cor_data)) {
        temp_meta_mic_links <- data.frame(source = colnames(cluster_cor_data)[i],target = rownames(cluster_cor_data))
        meta_mic_links <- rbind(meta_mic_links,temp_meta_mic_links)
      }
      meta_mic_links <- meta_mic_links[-1,]
      meta_mic_links$group <- "meta_mic"
      for (i in 1:ncol(cluster_cor_data)) {
        for (j in 1:nrow(cluster_cor_data)) {
          if(cluster_cor_data[j,i] < 0){
            meta_mic_links$group[(i-1)*nrow(cluster_cor_data) + j] <- "meta_mic_neg"
          }else{
            meta_mic_links$group[(i-1)*nrow(cluster_cor_data) + j] <- "meta_mic_pos"
          }
        }
      }
      
      ### Select microbes clusters according to cluster_cor_data
      mic_clu_names <- rownames(cluster_cor_data)
      mic_filter_data <- as.data.frame(mic_cluster_data[,colnames(mic_cluster_data) %in% mic_clu_names]) 
      colnames(mic_filter_data) <- mic_clu_names
      rownames(mic_filter_data) <- rownames(mic_cluster_data)
      mic_links <- data.frame(source = 1,target = 1)
      for (i in colnames(mic_filter_data)) {
        temp <- as.data.frame(mic_filter_data[,i])
        colnames(temp) <- i
        rownames(temp) <- rownames(mic_filter_data)
        temp <- abs(temp)
        temp <- as.data.frame(apply(temp,2,sort,decreasing = T)) 
        ### Top 3 microbes
        if(nrow(temp) > 3){
          temp_top3 <- as.data.frame(temp[1:3,])
          colnames(temp_top3) <- colnames(temp)
          rownames(temp_top3) <- rownames(temp)[1:3]
        }else{
          temp_top3 <- temp
        }
        temp_top3$modulte_name <- colnames(temp)[1]
        temp_mic_links <- data.frame(source = colnames(temp)[1],target = rownames(temp_top3))
        mic_links <- rbind(mic_links,temp_mic_links)
      }
      mic_links <- mic_links[-1,]
      mic_links$group <- "mic"
      
      links <- rbind(meta_links,meta_mic_links,mic_links)
      links$weight <- 10
      
      nodes <- data.frame(name = c(links$source, links$target) %>% unique()) 
      nodes$group <- "mic"
      nodes$group[1:nrow(meta_links)] <- "meta"
      nodes$group[nrow(meta_links) + 1:ncol(cluster_cor_data)] <- "meta_clusters"
      nodes$group[nrow(meta_links) + ncol(cluster_cor_data) + 1:nrow(cluster_cor_data)] <- "mic_clusters"
      
      
      links$IDsource <- match(links$source, nodes$name)-1 
      links$IDtarget <- match(links$target, nodes$name)-1
      colorJS <- 'd3.scaleOrdinal(["#6300ea","#4bacc6","#d7191c","#00b692","#6300ea","#00b692"])'
      sn <- sankeyNetwork(Links = links, 
                          Nodes = nodes,
                          Source = "IDsource", 
                          Target = "IDtarget",
                          Value = "weight",  
                          # NodeID = "name", 
                          LinkGroup = "group",
                          NodeGroup = "group",
                          colourScale = colorJS, 
                          fontSize = 12,
                          sinksRight = TRUE)
      # setwd("./results/Inter-Cor/hierarchical/")
      saveNetwork(sn,file = paste0("sankeyNetwork_",phenotype,".html"))
      webshot(paste0("sankeyNetwork_",phenotype,".html"), paste0("sankeyNetwork_",phenotype,".pdf"))
      pdf_subset(paste0("sankeyNetwork_",phenotype,".pdf"),pages = 1:1, output = paste0("sn_",phenotype,".pdf"))
    }
    ### WGCNA
    else{
      ### Select metabolite clusters according to cluster_cor_data
      meta_clu_names <- colnames(cluster_cor_data)
      meta_filter_data <- meta_cluster_data[meta_cluster_data$cluster_name %in% meta_clu_names,] 
      meta_links <- data.frame(source = 1,target = 1)
      for (i in meta_clu_names) {
        meta_module_i <- meta_filter_data[meta_filter_data$cluster_name %in% i,]
        ### Top 3 metabolites
        if(nrow(meta_module_i) > 3){
          temp_meta_links <- meta_module_i[1:3,]
          temp_meta_links <- data.frame(source = rownames(temp_meta_links),target = i)
          meta_links <- rbind(meta_links,temp_meta_links)
        }else{
          temp_meta_links <- data.frame(source = rownames(meta_module_i),target = i)
          meta_links <- rbind(meta_links,temp_meta_links)
        }
      }
      meta_links <- meta_links[-1,]
      
      meta_links$group <- "meta"
      
      ### meta_mic_links
      meta_mic_links <- data.frame(source = 1,target = 1)
      for (i in 1:ncol(cluster_cor_data)) {
        temp_meta_mic_links <- data.frame(source = colnames(cluster_cor_data)[i],target = rownames(cluster_cor_data))
        meta_mic_links <- rbind(meta_mic_links,temp_meta_mic_links)
      }
      meta_mic_links <- meta_mic_links[-1,]
      meta_mic_links$group <- "meta_mic"
      for (i in 1:ncol(cluster_cor_data)) {
        for (j in 1:nrow(cluster_cor_data)) {
          if(cluster_cor_data[j,i] < 0){
            meta_mic_links$group[(i-1)*nrow(cluster_cor_data) + j] <- "meta_mic_neg"
          }else{
            meta_mic_links$group[(i-1)*nrow(cluster_cor_data) + j] <- "meta_mic_pos"
          }
        }
      }
      
      ### Select microbes clusters according to cluster_cor_data
      mic_clu_names <- rownames(cluster_cor_data)
      mic_filter_data <- mic_cluster_data[mic_cluster_data$cluster_name %in% mic_clu_names,] 
      mic_links <- data.frame(source = 1,target = 1)
      for (i in mic_clu_names) {
        mic_module_i <- mic_filter_data[mic_filter_data$cluster_name %in% i,]
        ### Top 3 microbes
        if(nrow(mic_module_i) > 3){
          temp_mic_links <- mic_module_i[1:3,]
          temp_mic_links <- data.frame(source = i,target = rownames(temp_mic_links))
          mic_links <- rbind(mic_links,temp_mic_links)
        }else{
          temp_mic_links <- data.frame(source = i,target = rownames(mic_module_i))
          mic_links <- rbind(mic_links,temp_mic_links)
        }
      }
      mic_links <- mic_links[-1,]
      mic_links$group <- "mic"
      
      links <- rbind(meta_links,meta_mic_links,mic_links)
      links$weight <- 10
      
      nodes <- data.frame(name = c(links$source, links$target) %>% unique()) 
      nodes$group <- "mic"
      nodes$group[1:nrow(meta_links)] <- "meta"
      nodes$group[nrow(meta_links) + 1:ncol(cluster_cor_data)] <- "meta_clusters"
      nodes$group[nrow(meta_links) + ncol(cluster_cor_data) + 1:nrow(cluster_cor_data)] <- "mic_clusters"
      
      
      links$IDsource <- match(links$source, nodes$name)-1 
      links$IDtarget <- match(links$target, nodes$name)-1
      colorJS <- 'd3.scaleOrdinal(["#6300ea","#4bacc6","#d7191c","#00b692","#6300ea","#00b692"])'
      sn <- sankeyNetwork(Links = links, 
                          Nodes = nodes,
                          Source = "IDsource", 
                          Target = "IDtarget",
                          Value = "weight",  
                          NodeID = "name",
                          LinkGroup = "group",
                          NodeGroup = "group",
                          colourScale = colorJS, 
                          fontSize = 12,
                          sinksRight = TRUE)
      
      saveNetwork(sn,file = paste0("sankeyNetwork_",phenotype,".html"))
      webshot(paste0("sankeyNetwork_",phenotype,".html"), paste0("sankeyNetwork_",phenotype,".pdf"))
      pdf_subset(paste0("sankeyNetwork_",phenotype,".pdf"),pages = 1:1, output = paste0("sn_",phenotype,".pdf"))
      }
  }
  
  ## categorical phenotypes
  if(!is.na(phenoData) && phenoDataType == "categorical"){
    # sankeyPlot for each group
    for (i in 1:length(unique(phenoData[,1]))) {
      tempmetaClu <- list(metaCluData[[2*i+1]],metaCluData[[2*i+2]])
      tempmicClu <- list(micCluData[[2*i+1]],micCluData[[2*i+2]])
      if("type_All" %in% names(cluCorData)){  # results from gramm
        tempcluCor <- list(cluCorData[[4*i+1]],cluCorData[[4*i+2]],cluCorData[[4*i+3]])
      }else{
        tempcluCor <- list(cluCorData[[3*i+1]],cluCorData[[3*i+2]],cluCorData[[3*i+3]])
      }
      
      sankey_func(metaClu = tempmetaClu,micClu = tempmicClu,cluCor = tempcluCor, phenotype = unique(phenoData[,1])[i])
    }
    # sankeyPlot for all samples
    sankey_func(metaClu = metaCluData,micClu = micCluData,cluCor = cluCorData, phenotype = "All")
  }
  ## continuous or no phenotypes
  else{
    sankey_func(metaClu = metaCluData,micClu = micCluData,cluCor = cluCorData, phenotype = "All")
  }
}

