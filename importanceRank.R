# function - importanceRank #
########################################################################
# File: importanceRank.R
# Aim : Draw a histogram of the importance of variables in each module 
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
##    cluData --- A list included the classified metabolome or microbiome modules     ##
##                results from metabolitesCluster.R or microbesCluster.R              ##
##    clusterName --- The module name(e.g. "Meta1")                                   ##
##    clustergroup --- The module name belong to which kinds of module                ##
##                Default: All                                                        ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
##    Graph results in PDF format                                                     ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
importanceRank <- function(cluData,clusterName,clustergroup = "All")
{
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  # All group
  if(clustergroup == "All"){
    cluster_data <- cluData[[1]]
  }
  # sub group
  else{
    cluster_data <- cluData[[paste0("individual_metabolites_",clustergroup)]]
    
    # If this group does not exist, stop.
    if(is.null(cluster_data)){
      cluster_data <- cluData[[paste0("individual_microbes_",clustergroup)]]
    }
    if(is.null(cluster_data)){
      stop("This group does not exist!")
    }
  }
  
  if(grepl("Meta",cluster_data$cluster_name[1])){
    fill_cor <- "#6300ea"
  }else{
    fill_cor <- "#00b692"
  }
  
  ### PCA¡¢PCoA
  if(all(cluster_data$cluster_name == "all")){
    plot_data <- as.data.frame(cluster_data[,clusterName]) 
    rownames(plot_data) <- rownames(cluster_data)
    colnames(plot_data) <- clusterName
    plot_data <- as.data.frame(apply(plot_data,2,sort,decreasing = F)) 

    name <- rownames(plot_data)
    # name <- factor(name,levels=name)
    df <- data.frame('name' = name,'importance' = plot_data[,1])
    df$importance <- abs(df$importance)
    df <- df[order(df$importance,decreasing = F),]
    df$importance <- as.numeric(df$importance)
    ### Top 10 metabolites or microbes
    if(nrow(df) > 10){
      df <- df[(nrow(df)-10):nrow(df),]
    }
    df$name <- factor(df$name,levels = df$name)
    
    ### plot
    plot <- ggplot(data = df, aes(x = name, y = importance)) +
      geom_bar(stat = "identity", fill = fill_cor) + 
      labs(title = paste0(clusterName, " module")  ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(limits = c(0,(max(df$importance) + 0.1)),expand = c(0,0))+ 
      coord_flip() 
    ggsave(paste0("./results/Inter-Cor/hierarchical/",clusterName,"_module_importanceRank_",clustergroup,".pdf"),width = 7,height = 5)
  }
  ### WGCNA
  else{
    # cluster_data <- cluData[[3]]
    # clusterName <- "Mic6"
    plot_data <- cluster_data[which(cluster_data$cluster_name == clusterName),]
    name <- rownames(plot_data)
    # name <- factor(name,levels=name)
    df <- data.frame('name' = name,'importance' = plot_data$kIN)
    df$importance <- abs(df$importance)
    df <- df[order(df$importance,decreasing = F),]
    df$importance <- as.numeric(df$importance)
    ### Top 10 metabolites or microbes
    if(nrow(df) > 10){
      df <- df[(nrow(df)-9):nrow(df),]
    }
    df$name <- factor(df$name,levels = df$name)
    ### plot
    plot <- ggplot(data = df, aes(x = name, y = importance)) +
      geom_bar(stat = "identity", fill = fill_cor) + 
      labs(title = paste0(clusterName, " module")  ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_y_continuous(limits = c(0,(max(df$importance) + 0.05)),expand = c(0,0)) +
      coord_flip() 
    ggsave(paste0("./results/Inter-Cor/hierarchical/",clusterName,"_module_importanceRank_",clustergroup,".pdf"),width = 7,height = 5)
  }
}

