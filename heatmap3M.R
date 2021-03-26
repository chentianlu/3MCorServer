# function - heatmap3M #
########################################################################
# File: heatmap3M.R
# Aim : Draw the heatmap of correlation between metabolome and microbiome modules
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
########################################################################################
## Input:                                                                             ##
##		clustersCorData --- A list from clustersCor.R included coefficient, p value     ##
##                        and p.adjust value of the correlation between metabolome    ##
##                        modules and microbiome modules from clustersCor.R           ##
##		metaCluData --- A list included the classified metabolome modules results from  ##
##                   metabolitesCluster.R                                             ##
##		micCluData --- A list included the classified microbiome modules results from   ##
##                   microbesCluster.R                                                ##
##    phenotype --- The number of nearest neighbours to use in knn imputation         ##
##                Default: NA                                                         ##
##    phenoData --- A dataframe of phenotype data                                     ##
##                Default: NA  (rows: samples,columns: phenotype)                     ##
##    phenoDataType --- Phenotype data type                                           ##
##                Default: continuous                                                 ##
##    fdrThreshold --- FDR threshold for filtering modules by correlation analysis    ##
##                  between modules                                                   ##
##                  Default: 0.99                                                     ##
##    cluster2pheno --- Which kind of module is correlated to and phenotype data      ##
##                  Default: mic                                                      ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
## 	  Graph results in PDF format                                                     ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
heatmap3M <- function(clustersCorData,metaCluData,micCluData,phenotype = NA,phenoData = NA,phenoDataType = "continuous",fdrThreshold = 0.99,cluster2pheno = "mic")
{
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  if(cluster2pheno == "mic"){
    # clustersCorData = clustersCor_res
    # micCluData = microbesCluster_res
    # phenotype = "BMI"
    # fdrThreshold = 0.1
    
    ### Select only significant associations to show in heatmap (to make it visually
    ### comprehensible) I.e. exclude microbes modules and metabolite modules with not at
    ### least one significant association. In this example, all rows/columns are
    ### kept.
    mic_meta_cor <- as.data.frame(clustersCorData[["r_All"]])
    mic_meta_p.adjust <- as.data.frame(clustersCorData[["p.adjust_All"]]) 
    tmp = mic_meta_p.adjust < fdrThreshold
    issig = rowSums (tmp, na.rm = T)
    nmic = na.omit (names (issig [issig > 0]))
    length (nmic)
    issig = colSums (tmp, na.rm=T)
    nmeta = na.omit (names (issig [issig > 0]))
    length (nmeta)
    
    ### create matrix for heatmap
    plotmat = as.data.frame(mic_meta_cor[nmic,nmeta]) 
    rownames (plotmat) <- nmic
    colnames (plotmat) <- nmeta
    
    ### Plot heatmap
    
    if(nrow(plotmat) < 2 && ncol(plotmat) >= 2){
      row_dendrogram = FALSE
      col_dendrogram = TRUE
    }else if(nrow(plotmat) >= 2 && ncol(plotmat) < 2){
      row_dendrogram = FALSE
      col_dendrogram = FALSE
    }else if(nrow(plotmat) < 2 && ncol(plotmat) < 2){
      row_dendrogram = FALSE
      col_dendrogram = FALSE
      stop("Unable to generate heat map due to too few significant correlation modules.")
    }else{
      row_dendrogram = FALSE
      col_dendrogram = TRUE
    }
    
    ### create matrix with significance stars
    plotmat_p = as.data.frame(mic_meta_p.adjust [nmic, nmeta])### p-values to make stars for heatmap
    rownames (plotmat_p) <- nmic
    colnames (plotmat_p) <- nmeta
    stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
    rownames (stars) <- nmic
    colnames (stars) <- nmeta
    for (z in 1:ncol (stars)) {
      for (j in 1:nrow (stars)) {
        if (plotmat_p [j, z] < 0.05) {
          stars [j, z] = "+"
        }
        if (plotmat_p [j, z] < 0.01) {
          stars [j, z] = "*"
        }
        if (plotmat_p [j, z] < 0.001) {
          stars [j, z] = "**"
        }
      }
    }
    
    
    
    ############ Exist phenotypes ############
    if(!is.na(phenotype)){
      ### make sidebar with phenotype associations
      phenobar = matrix (NA, ncol = length (phenotype), nrow = length (nmic))
      colnames (phenobar) = phenotype
      rownames (phenobar) = nmic
      pheno_cluster_microbes <- micCluData[[paste0(phenotype,"_cluster_microbes")]]
      for (j in phenotype) {
        tmp_p.adjust = pheno_cluster_microbes[,3]
        names(tmp_p.adjust) = rownames(pheno_cluster_microbes)
        for (k in nmic) {
          if (tmp_p.adjust [k] >= 0.1) {        ### No significant association
            phenobar [k, j] = "No_significant_association"
          } else {
            if ( pheno_cluster_microbes[k,paste0(phenotype,"_estimate")] > 0) {        ### Positive association
              phenobar [k, j] = "Positive_association"
            } else if (pheno_cluster_microbes[k,paste0(phenotype,"_estimate")] < 0){        ### Negative association
              phenobar[k, j] = "Negative_association"
            } else {
              print(paste0(phenobar[k, j]," errors."))
            }
          }
          
        }
      }
      
      ### Rename column names of phenobar
      colnames (phenobar) = phenotype
      
      ### annotation_row and ann_colors
      annotation_row = data.frame(phenobar)
      ann_colors = list(
        assign(colnames(annotation_row)[1],c(No_significant_association = "grey40",Positive_association = "#d7191c",Negative_association = "#4bacc6"))
      )
      names(ann_colors) = phenotype
      
      # if(){
      #   map_color
      # }
      
      ### Note, normally one would cluster the rows and columns in the heatmap as
      ### shown below but for the paper (Pedersen et al, 2016) we needed to group the
      ### KEGG modules by biological similarity to make a higher-level annotation for
      ### the figure and thus made a manual arrangement.
      
      ### The breaks parameter redefines the color bar range and divides the color range according to the break range
      bk <- c(seq(-1,1,by=0.01))
      pheatmap(plotmat,
               # scale="row",
               clustering_method = "average",
               color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
               breaks = bk,
               border = FALSE,
               show_rownames = T,
               show_colnames = T,
               display_numbers = stars,
               fontsize_number = 20,
               number_color = "black",
               cluster_row = row_dendrogram,
               cluster_cols = col_dendrogram,
               annotation_row = annotation_row,
               annotation_colors = ann_colors,
               angle_col = "90",filename = "./results/Inter-Cor/hierarchical/heatmap3M_All.pdf")
    }
    ############ No phenotypes ############
    else{
      
      ### The breaks parameter redefines the color bar range and divides the color range according to the break range
      bk <- c(seq(-1,1,by=0.01))
      
      pheatmap(plotmat,
               # scale="row",
               clustering_method = "average",
               color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
               breaks = bk,
               border = FALSE,
               show_rownames = T,
               show_colnames = T,
               display_numbers = stars,
               fontsize_number = 20,
               number_color = "black",
               cluster_row = row_dendrogram,
               cluster_cols = col_dendrogram,
               # annotation_row = annotation_row,
               # annotation_colors = ann_colors,
               angle_col = "90",filename = "./results/Inter-Cor/hierarchical/heatmap3M_All.pdf")
    }
  }
  else if(cluster2pheno == "meta"){
    # clustersCorData = clustersCor_res
    # micCluData = microbesCluster_res
    # phenotype = "BMI"
    # fdrThreshold = 0.1
    
    ### Select only significant associations to show in heatmap (to make it visually
    ### comprehensible) I.e. exclude microbes modules and metabolite modules with not at
    ### least one significant association. In this example, all rows/columns are
    ### kept.
    mic_meta_cor <- as.data.frame(clustersCorData[["r_All"]])
    mic_meta_p.adjust <- as.data.frame(clustersCorData[["p.adjust_All"]]) 
    tmp = mic_meta_p.adjust < fdrThreshold
    issig = rowSums (tmp, na.rm = T)
    nmic = na.omit (names (issig [issig > 0]))
    length (nmic)
    issig = colSums (tmp, na.rm=T)
    nmeta = na.omit (names (issig [issig > 0]))
    length (nmeta)
    
    ### create matrix for heatmap
    plotmat = as.data.frame(t(mic_meta_cor[nmic,nmeta]))
    rownames (plotmat) <- nmeta
    colnames (plotmat) <- nmic
    
    ### Plot heatmap
    
    if(nrow(plotmat) < 2 && ncol(plotmat) >= 2){
      row_dendrogram = FALSE
      col_dendrogram = TRUE
    }else if(nrow(plotmat) >= 2 && ncol(plotmat) < 2){
      row_dendrogram = FALSE
      col_dendrogram = FALSE
    }else if(nrow(plotmat) < 2 && ncol(plotmat) < 2){
      row_dendrogram = FALSE
      col_dendrogram = FALSE
      stop("Unable to generate heat map due to too few significant correlation modules.")
    }else{
      row_dendrogram = FALSE
      col_dendrogram = TRUE
    }
    
    ### create matrix with significance stars
    plotmat_p = as.data.frame(t(mic_meta_p.adjust [nmic, nmeta]))### p-values to make stars for heatmap
    rownames (plotmat_p) <- nmeta
    colnames (plotmat_p) <- nmic
    stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
    rownames (stars) <- nmeta
    colnames (stars) <- nmic
    for (z in 1:ncol (stars)) {
      for (j in 1:nrow (stars)) {
        if (plotmat_p [j, z] < 0.05) {
          stars [j, z] = "+"
        }
        if (plotmat_p [j, z] < 0.01) {
          stars [j, z] = "*"
        }
        if (plotmat_p [j, z] < 0.001) {
          stars [j, z] = "**"
        }
      }
    }
    
    
    
    ############ Exist phenotypes ############
    if(!is.na(phenotype)){
      ### make sidebar with phenotype associations
      phenobar = matrix (NA, ncol = length (phenotype), nrow = length (nmeta))
      colnames (phenobar) = phenotype
      rownames (phenobar) = nmeta
      pheno_cluster_metabolites <- metaCluData[[paste0(phenotype,"_cluster_metabolites")]]
      for (j in phenotype) {
        tmp_p.adjust = pheno_cluster_metabolites[,3]
        names(tmp_p.adjust) = rownames(pheno_cluster_metabolites)
        for (k in nmeta) {
          if (tmp_p.adjust [k] >= 0.1) {        ### No significant association
            phenobar [k, j] = "No_significant_association"
          } else {
            if ( pheno_cluster_metabolites[k,paste0(phenotype,"_estimate")] > 0) {        ### Positive association
              phenobar [k, j] = "Positive_association"
            } else if (pheno_cluster_metabolites[k,paste0(phenotype,"_estimate")] < 0){        ### Negative association
              phenobar[k, j] = "Negative_association"
            } else {
              print(paste0(phenobar[k, j]," errors."))
            }
          }
          
        }
      }
      
      ### Rename column names of phenobar
      colnames (phenobar) = phenotype
      
      ### annotation_row and ann_colors
      annotation_row = data.frame(phenobar)
      ann_colors = list(
        assign(colnames(annotation_row)[1],c(No_significant_association = "grey40",Positive_association = "#d7191c",Negative_association = "#4bacc6"))
      )
      names(ann_colors) = phenotype
      
      # if(){
      #   map_color
      # }
      
      ### Note, normally one would cluster the rows and columns in the heatmap as
      ### shown below but for the paper (Pedersen et al, 2016) we needed to group the
      ### KEGG modules by biological similarity to make a higher-level annotation for
      ### the figure and thus made a manual arrangement.
      
      ### The breaks parameter redefines the color bar range and divides the color range according to the break range
      bk <- c(seq(-1,1,by=0.01))
      pheatmap(plotmat,
               # scale="row",
               clustering_method = "average",
               color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
               breaks = bk,
               border = FALSE,
               show_rownames = T,
               show_colnames = T,
               display_numbers = stars,
               fontsize_number = 20,
               number_color = "black",
               cluster_row = row_dendrogram,
               cluster_cols = col_dendrogram,
               annotation_row = annotation_row,
               annotation_colors = ann_colors,
               angle_col = "90",filename = "./results/Inter-Cor/hierarchical/heatmap3M_All.pdf")
    }
    ############ No phenotypes ############
    else{
      
      ### The breaks parameter redefines the color bar range and divides the color range according to the break range
      bk <- c(seq(-1,1,by=0.01))
      
      pheatmap(plotmat,
               # scale="row",
               clustering_method = "average",
               color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
               breaks = bk,
               border = FALSE,
               show_rownames = T,
               show_colnames = T,
               display_numbers = stars,
               fontsize_number = 20,
               number_color = "black",
               cluster_row = row_dendrogram,
               cluster_cols = col_dendrogram,
               # annotation_row = annotation_row,
               # annotation_colors = ann_colors,
               angle_col = "90",filename = "./results/Inter-Cor/hierarchical/heatmap3M_All.pdf")
    }
  }
  
  if(phenoDataType == "categorical"){
    if(cluster2pheno == "mic"){
      # clustersCorData = clustersCor_res
      # micCluData = microbesCluster_res
      # phenotype = "BMI"
      # fdrThreshold = 0.1
      
      ### Select only significant associations to show in heatmap (to make it visually
      ### comprehensible) I.e. exclude microbes modules and metabolite modules with not at
      ### least one significant association. In this example, all rows/columns are
      ### kept.
      for (i in unique(phenoData[,1])) {
        mic_meta_cor <- as.data.frame(clustersCorData[[paste0("r_",i)]])
        mic_meta_p.adjust <- as.data.frame(clustersCorData[[paste0("p.adjust_",i)]]) 
        tmp = mic_meta_p.adjust < fdrThreshold
        issig = rowSums (tmp, na.rm = T)
        nmic = na.omit (names (issig [issig > 0]))
        length (nmic)
        issig = colSums (tmp, na.rm=T)
        nmeta = na.omit (names (issig [issig > 0]))
        length (nmeta)
        
        ### create matrix for heatmap
        plotmat = as.data.frame(mic_meta_cor[nmic,nmeta]) 
        rownames (plotmat) <- nmic
        colnames (plotmat) <- nmeta
        
        ### Plot heatmap
        
        if(nrow(plotmat) < 2 && ncol(plotmat) >= 2){
          row_dendrogram = FALSE
          col_dendrogram = TRUE
        }else if(nrow(plotmat) >= 2 && ncol(plotmat) < 2){
          row_dendrogram = FALSE
          col_dendrogram = FALSE
        }else if(nrow(plotmat) < 2 && ncol(plotmat) < 2){
          row_dendrogram = FALSE
          col_dendrogram = FALSE
          stop("Unable to generate heat map due to too few significant correlation modules.")
        }else{
          row_dendrogram = FALSE
          col_dendrogram = TRUE
        }
        
        ### create matrix with significance stars
        plotmat_p = as.data.frame(mic_meta_p.adjust [nmic, nmeta])### p-values to make stars for heatmap
        rownames (plotmat_p) <- nmic
        colnames (plotmat_p) <- nmeta
        stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
        rownames (stars) <- nmic
        colnames (stars) <- nmeta
        for (z in 1:ncol (stars)) {
          for (j in 1:nrow (stars)) {
            if (plotmat_p [j, z] < 0.05) {
              stars [j, z] = "+"
            }
            if (plotmat_p [j, z] < 0.01) {
              stars [j, z] = "*"
            }
            if (plotmat_p [j, z] < 0.001) {
              stars [j, z] = "**"
            }
          }
        }
        ### The breaks parameter redefines the color bar range and divides the color range according to the break range
        bk <- c(seq(-1,1,by=0.01))
        
        pheatmap(plotmat,
                 # scale="row",
                 clustering_method = "average",
                 color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
                 breaks = bk,
                 border = FALSE,
                 show_rownames = T,
                 show_colnames = T,
                 display_numbers = stars,
                 fontsize_number = 20,
                 number_color = "black",
                 cluster_row = row_dendrogram,
                 cluster_cols = col_dendrogram,
                 # annotation_row = annotation_row,
                 # annotation_colors = ann_colors,
                 angle_col = "90",filename = paste0("./results/Inter-Cor/hierarchical/heatmap3M_",i,".pdf"))
      }
    }
    else if(cluster2pheno == "meta"){
      # clustersCorData = clustersCor_res
      # micCluData = microbesCluster_res
      # phenotype = "BMI"
      # fdrThreshold = 0.1
      
      ### Select only significant associations to show in heatmap (to make it visually
      ### comprehensible) I.e. exclude microbes modules and metabolite modules with not at
      ### least one significant association. In this example, all rows/columns are
      ### kept.
      for (i in unique(phenoData[,1])) {
        mic_meta_cor <- as.data.frame(clustersCorData[[paste0("r_",i)]])
        mic_meta_p.adjust <- as.data.frame(clustersCorData[[paste0("p.adjust_",i)]]) 
        tmp = mic_meta_p.adjust < fdrThreshold
        issig = rowSums (tmp, na.rm = T)
        nmic = na.omit (names (issig [issig > 0]))
        length (nmic)
        issig = colSums (tmp, na.rm=T)
        nmeta = na.omit (names (issig [issig > 0]))
        length (nmeta)
        
        ### create matrix for heatmap
        plotmat = as.data.frame(t(mic_meta_cor[nmic,nmeta]))
        rownames (plotmat) <- nmeta
        colnames (plotmat) <- nmic
        
        ### Plot heatmap
        
        if(nrow(plotmat) < 2 && ncol(plotmat) >= 2){
          row_dendrogram = FALSE
          col_dendrogram = TRUE
        }else if(nrow(plotmat) >= 2 && ncol(plotmat) < 2){
          row_dendrogram = FALSE
          col_dendrogram = FALSE
        }else if(nrow(plotmat) < 2 && ncol(plotmat) < 2){
          row_dendrogram = FALSE
          col_dendrogram = FALSE
          stop("Unable to generate heat map due to too few significant correlation modules.")
        }else{
          row_dendrogram = FALSE
          col_dendrogram = TRUE
        }
        
        ### create matrix with significance stars
        plotmat_p = as.data.frame(t(mic_meta_p.adjust [nmic, nmeta]))### p-values to make stars for heatmap
        rownames (plotmat_p) <- nmeta
        colnames (plotmat_p) <- nmic
        stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
        rownames (stars) <- nmeta
        colnames (stars) <- nmic
        for (z in 1:ncol (stars)) {
          for (j in 1:nrow (stars)) {
            if (plotmat_p [j, z] < 0.05) {
              stars [j, z] = "+"
            }
            if (plotmat_p [j, z] < 0.01) {
              stars [j, z] = "*"
            }
            if (plotmat_p [j, z] < 0.001) {
              stars [j, z] = "**"
            }
          }
        }
        
        ### The breaks parameter redefines the color bar range and divides the color range according to the break range
        bk <- c(seq(-1,1,by=0.01))
        
        pheatmap(plotmat,
                 # scale="row",
                 clustering_method = "average",
                 color = colorRampPalette(c("#4bacc6", "#FFFFFF" ,"#d7191c"))(length(bk)),
                 breaks = bk,
                 border = FALSE,
                 show_rownames = T,
                 show_colnames = T,
                 display_numbers = stars,
                 fontsize_number = 20,
                 number_color = "black",
                 cluster_row = row_dendrogram,
                 cluster_cols = col_dendrogram,
                 # annotation_row = annotation_row,
                 # annotation_colors = ann_colors,
                 angle_col = "90",filename = paste0("./results/Inter-Cor/hierarchical/heatmap3M_",i,".pdf"))
      }
    }
  }
}