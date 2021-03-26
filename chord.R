# function - chord #
########################################################################
# File: chord.R
# Aim : Draw the chord of correlation between metabolome and microbiome modules
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
##    cluData --- A list from clustersCor.R included coefficient, p value and         ##
##                p.adjust value of the correlation between metabolome modules and    ##
##                microbiome modules from clustersCor.R                               ##
##    phenoData --- A dataframe of phenotype data                                     ##
##                Default: NA  (rows: samples,columns: phenotype)                     ##
##    phenoDataType --- Phenotype data type                                           ##
##                Default: continuous                                                 ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
## 	  Graph results in PDF format                                                     ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
chord <- function(cluData,phenoData = NA,phenoDataType = NA)
{
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  chord_func <- function(cludf,phenotype){
    cluster_data <- cludf[[1]]
    CairoPDF(paste0('./results/Inter-Cor/hierarchical/Chord_',phenotype,'.pdf'),width = 8,height = 8)
    migration <- as.matrix(cluster_data) 
    # sectorcol <- rainbow(nrow(migration) * ncol(migration))
    linkcol <- c()
    for (i in 1:ncol(cluster_data)) {
      for(j in 1:nrow(cluster_data)){
        if(cluster_data[j,i] > 0) {
          linkcol <- append(linkcol,"#d7191c")
        }
        else{
          linkcol <- append(linkcol,"#4bacc6")
        }
      }
    }
    cluster_data2 <- cludf[[2]]
    linktran <- c()
    for (i in 1:ncol(cluster_data2)) {
      for(j in 1:nrow(cluster_data2)){
        temp_tran <- cluster_data2[j,i]
        linktran <- append(linktran,temp_tran)
      }
    }
    
    circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
    if(ncol(cluster_data) <= 2){
      chordDiagram(migration, col = linkcol,annotationTrack = c("name", "grid"))
    }else{
      chordDiagram(migration, col = linkcol,transparency = linktran,annotationTrack = c("name", "grid"))
    }
    legend("topright",legend = c("Positive association","Negative association"),col = c("#d7191c","#4bacc6"),
           bty = "n",fill = c("#d7191c","#4bacc6"),border = c("#d7191c","#4bacc6"))
    # legend(-1, 1.9, c("sin", "cos", "tan"), col = c(3, 4, 6),
    #        text.col = "green4", lty = c(2, -1, 1), pch = c(NA, 3, 4),
    #        merge = TRUE, bg = "gray90")
    circos.clear()
    dev.off()
  }
  
  ## categorical phenotypes
  if(!is.na(phenoData) && phenoDataType == "categorical"){
    # chord for each group
    for (i in 1:length(unique(phenoData[,1]))) {
      if("type_All" %in% names(cluData)){  # results from gramm
        new_clustersCor_res <- list(clustersCor_res[[4*i+1]],clustersCor_res[[4*i+2]],clustersCor_res[[4*i+3]])
      }else{
        new_clustersCor_res <- list(clustersCor_res[[3*i+1]],clustersCor_res[[3*i+2]],clustersCor_res[[3*i+3]])
      }
      chord_func(cludf = new_clustersCor_res,phenotype = unique(phenoData[,1])[i])
    }
    # chord for all samples
    new_clustersCor_res <- list(clustersCor_res[[1]],clustersCor_res[[2]],clustersCor_res[[3]])
    chord_func(cludf = new_clustersCor_res,phenotype = "All")
  }
  
  ## continuous or no phenotypes
  else{
    chord_func(cluData,phenotype = "All")
  }
}

