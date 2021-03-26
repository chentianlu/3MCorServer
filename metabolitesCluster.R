# function - metabolitesCluster #
########################################################################
# File: metabolitesCluster.R
# Aim : Clustering for metabolome data
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
######################################################################################
## Input:                                                                           ##
##		metaData --- A dataframe of microbes                                          ##
##                 (rows: samples,columns: microbes)                                ## 
##    phenoData --- A dataframe of phenotype data                                   ##
##                 (rows: samples,columns: phenotype)                               ##
##    phenoDataType --- Phenotype data type                                         ##
##    dimReduMethod --- Dimension reduction method                                  ##
##                      Default: WGCNA                                              ##
##    corMethod --- Correlation method of associating individual microbes or        ##
##                  microbiome cluster with phenotype                               ##
##                  Default: spearman                                               ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
## 	  list_res --- A list included the classified metabolome modules results        ##                                 
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
metabolitesCluster <- function(metaData,phenoData,phenoDataType,dimReduMethod = "WGCNA",corMethod = "spearman")
{
  corFuntion = "bicor"
  clusterMethod = "average"
  corOptionsList = list(use = 'pairwise.complete.obs')
  corOptionsStr = "use = 'pairwise.complete.obs'"
  networkType = "signed"
  mergingThresh = 0.20
  maxdissimilarity = 0.20
  minModuleSize = 3
  softPower = 13
  
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  metabolites <- metaData
  phenotypes <- phenoData
  
  #---------------------------------------------#
  #--------------------WGCNA--------------------#
  #---------------------------------------------#
  if(dimReduMethod == "WGCNA"){
    #######################################################
    ###### Step 1 - Identify clusters of metabolites ######
    #######################################################
    
    ###Specify data and parameters
    # dat_tmp = log2(metabolites)  ### work in logarithmic space
    dat_tmp = metabolites
    
    ### Calculate weighted adjacency matrix
    network_adjacency = adjacency(dat_tmp, power = softPower, type = networkType, corFnc = corFuntion, corOptions = corOptionsStr)
    colnames(network_adjacency) = rownames(network_adjacency) = colnames(dat_tmp)
    
    ### Define dissimilarity based on topological overlap
    dissTOM = TOMdist(network_adjacency, TOMType = networkType)
    colnames(dissTOM) = rownames(dissTOM) = colnames(dat_tmp)
    
    ### Hierarchical clustering
    metaTree = flashClust(as.dist(dissTOM), method = clusterMethod)
    
    ### Define modules by cutting branches
    moduleLabels1 = cutreeDynamic(dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
    moduleLabels1 = labels2colors(moduleLabels1)
    
    ### Automatically merge highly correlated modules
    merge = mergeCloseModules(dat_tmp, moduleLabels1, corFnc = corFuntion, corOptions = corOptionsList, cutHeight = mergingThresh)
    
    ### Determine resulting merged module colors
    moduleLabels2 = merge$colors
    
    ### Establish eigengenes of the newly merged modules, used for cluster overall abundances
    MEs = merge$newMEs
    
    ### Choose final module assignments
    moduleColorsMeta = moduleLabels2
    names(moduleColorsMeta) = colnames(dat_tmp)
    MEsMeta = orderMEs(MEs)
    rownames(MEsMeta) = rownames(dat_tmp)
    
    ##Create the cluster_mapping_file file
    cluster_mapping_file_m <- data.frame(New_Name=paste0("Meta",1:length(MEsMeta)),Description=c("None"))
    row.names(cluster_mapping_file_m) <- paste0("Meta_",colnames(MEsMeta))
    cluster_mapping_file_m$label = sapply(rownames(cluster_mapping_file_m), function(x) paste(cluster_mapping_file_m[x, "New_Name"], cluster_mapping_file_m[x, "Description"], sep = ": "))
    
    ### Determine relevant descriptive statistics of established clusters
    ### kIN: within-module connectivity, determined by summing connectivity with all
    ###      other metabolites in the given cluster.
    ### kME: bicor-correlation between the metabolite profile and module eigenvector; 
    ### both measures of intramodular hub-metabolite status.
    kIN <-      vector(length = ncol(dat_tmp)); names(kIN) = colnames(dat_tmp)
    kME <-      vector(length = ncol(dat_tmp)); names(kME) = colnames(dat_tmp)
    modules <-  vector(length = ncol(dat_tmp)); names(modules) = colnames(dat_tmp)
    
    for(module in names(table(moduleColorsMeta))) {   
      all.metabolites = names(dat_tmp)
      inModule =(moduleColorsMeta == module)
      module.metabolites = names(moduleColorsMeta[inModule])
      modules[module.metabolites] = module 
      kIN[module.metabolites] = sapply(module.metabolites, function(x) sum(network_adjacency[x, module.metabolites]) - 1)
      datKME = signedKME(dat_tmp, MEsMeta, corFnc = corFuntion, corOptions = corOptionsStr)
      rownames(datKME) = colnames(dat_tmp)
      kME[module.metabolites] = datKME[module.metabolites, paste("kME", module, sep = "")]   
    }
    output = data.frame("module" = modules, "kME" = kME, "kIN" = kIN,"cluster_name" = sapply(modules, function(m) cluster_mapping_file_m[paste0("Meta_ME", m), "New_Name"]))
    
    
    ############ Exist continuous phenotype variables ############
    if(is.data.frame(phenotypes) & phenoDataType == "continuous"){
      ###########################################################################
      ###### Step 2 - Link individual metabolite to phenotypes of interest ######
      ###########################################################################
      
      tmpMat = array(NA, c(ncol(metabolites), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = colnames(metabolites)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating individual metabolites with phenotypes
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(metabolites, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      ### Sort by cluster_name and then decreasing values of kIN
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        output <- cbind(output, tmp_output)
      }
      output = output[with(output,  order(cluster_name, -kIN)), ]
      
      ### Write to file
      write.table(output, file = "./results/Inter-Cor/hierarchical/individual_metabolites.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      
      ### return 
      res1 <- output
      rm(output, tmpMat, dat_tmp)
      
      ########################################
      ###### step3 - Merge the results #######
      ########################################
      
      ### Create a joint data frame of metabolite cluster eigengene
      ### equivalents/effective abundances('MEsMeta')
      MEsMeta = MEsMeta[rownames(MEsMeta),]
      
      ### Rename module names(columnames) from 'colors' to numbers
      colnames(MEsMeta) = cluster_mapping_file_m[paste0("Meta_", colnames(MEsMeta)), "New_Name"]
      
      ### Save MEsMeta to file, order by module number
      write.table(MEsMeta[,order(colnames(MEsMeta))], file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      
      ### return 
      res2 <- MEsMeta[,order(colnames(MEsMeta))]
      
      ########################################################################
      ###### Step 4 - Link metabolite clusters to phenotype of interest ######
      ########################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_pheno <- list() ### Data structure for storing results of correlation tests under different setups
      tmpMat = array(NA, c(ncol(MEsMeta), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = names(MEsMeta)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating metabolite clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMeta, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_pheno[["meta"]] <- tmpMat
      rm(tmpMat)
      
      ######################################################
      ######## Step 5 - Save phenotype associations ######## 
      ######################################################
      
      ### Save the age association of metabolite clusters calculated in Step 4.
      tmpMat = cor_pheno[["meta"]]
      list_res <- list(individual_metabolites = res1, MEs_metabolites_clusters = res2)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_metabolites")]] <- tmp_output
      }
      ### return 
      return(list_res)
    }
    
    
    ############ Exist categorical phenotype variables ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############
      
      output = output[with(output,  order(cluster_name, -kIN)), ]
      write.table(output, file = "./results/Inter-Cor/hierarchical/individual_metabolites.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      res1 <- output
      
      ### Create a joint data frame of metabolite cluster eigengene
      ### equivalents/effective abundances('MEsMeta')
      MEsMeta = MEsMeta[rownames(MEsMeta),]
      
      ### Rename module names(columnames) from 'colors' to numbers
      colnames(MEsMeta) = cluster_mapping_file_m[paste0("Meta_", colnames(MEsMeta)), "New_Name"]
      
      ### Save MEsMeta to file, order by module number
      write.table(MEsMeta[,order(colnames(MEsMeta))], file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMeta[,order(colnames(MEsMeta))]
      
      list_res <- list(individual_metabolites = res1, MEs_metabolites_clusters = res2)
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(metabolites)  ### work in logarithmic space
        
        temp_metabolites_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        
        dat_tmp = metabolites[temp_metabolites_rownames,]
        
        ### Calculate weighted adjacency matrix
        network_adjacency = adjacency(dat_tmp, power = softPower, type = networkType, corFnc = corFuntion, corOptions = corOptionsStr)
        colnames(network_adjacency) = rownames(network_adjacency) = colnames(dat_tmp)
        
        ### Define dissimilarity based on topological overlap
        dissTOM = TOMdist(network_adjacency, TOMType = networkType)
        colnames(dissTOM) = rownames(dissTOM) = colnames(dat_tmp)
        
        ### Hierarchical clustering
        metaTree = flashClust(as.dist(dissTOM), method = clusterMethod)
        
        ### Define modules by cutting branches
        moduleLabels1 = cutreeDynamic(dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
        moduleLabels1 = labels2colors(moduleLabels1)
        
        ### Automatically merge highly correlated modules
        merge = mergeCloseModules(dat_tmp, moduleLabels1, corFnc = corFuntion, corOptions = corOptionsList, cutHeight = mergingThresh)
        
        ### Determine resulting merged module colors
        moduleLabels2 = merge$colors
        
        ### Establish eigengenes of the newly merged modules, used for cluster overall abundances
        MEs = merge$newMEs
        
        ### Choose final module assignments
        moduleColorsMeta = moduleLabels2
        names(moduleColorsMeta) = colnames(dat_tmp)
        MEsMeta = orderMEs(MEs)
        rownames(MEsMeta) = rownames(dat_tmp)
        
        ##Create the cluster_mapping_file file
        cluster_mapping_file_m <- data.frame(New_Name=paste0("Meta",1:length(MEsMeta)),Description=c("None"))
        row.names(cluster_mapping_file_m) <- paste0("Meta_",colnames(MEsMeta))
        cluster_mapping_file_m$label = sapply(rownames(cluster_mapping_file_m), function(x) paste(cluster_mapping_file_m[x, "New_Name"], cluster_mapping_file_m[x, "Description"], sep = ": "))
        
        ### Determine relevant descriptive statistics of established clusters
        ### kIN: within-module connectivity, determined by summing connectivity with all
        ###      other metabolites in the given cluster.
        ### kME: bicor-correlation between the metabolite profile and module eigenvector; 
        ### both measures of intramodular hub-metabolite status.
        kIN <-      vector(length = ncol(dat_tmp)); names(kIN) = colnames(dat_tmp)
        kME <-      vector(length = ncol(dat_tmp)); names(kME) = colnames(dat_tmp)
        modules <-  vector(length = ncol(dat_tmp)); names(modules) = colnames(dat_tmp)
        
        for(module in names(table(moduleColorsMeta))) {   
          all.metabolites = names(dat_tmp)
          inModule =(moduleColorsMeta == module)
          module.metabolites = names(moduleColorsMeta[inModule])
          modules[module.metabolites] = module 
          kIN[module.metabolites] = sapply(module.metabolites, function(x) sum(network_adjacency[x, module.metabolites]) - 1)
          datKME = signedKME(dat_tmp, MEsMeta, corFnc = corFuntion, corOptions = corOptionsStr)
          rownames(datKME) = colnames(dat_tmp)
          kME[module.metabolites] = datKME[module.metabolites, paste("kME", module, sep = "")]   
        }
        temp_output = data.frame("module" = modules, "kME" = kME, "kIN" = kIN,"cluster_name" = sapply(modules, function(m) cluster_mapping_file_m[paste0("Meta_ME", m), "New_Name"]))
        temp_output = temp_output[with(temp_output,  order(cluster_name, -kIN)), ]
        write.table(temp_output, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- temp_output
        
        ### Create a joint data frame of metabolite cluster eigengene
        ### equivalents/effective abundances('MEsMeta')
        MEsMeta = MEsMeta[rownames(MEsMeta),]
        
        ### Rename module names(columnames) from 'colors' to numbers
        colnames(MEsMeta) = cluster_mapping_file_m[paste0("Meta_", colnames(MEsMeta)), "New_Name"]
        
        ### Save MEsMeta to file, order by module number
        write.table(MEsMeta[,order(colnames(MEsMeta))], file = paste0("./results/Inter-Cor/hierarchical/MEs_metabolites_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMeta[,order(colnames(MEsMeta))]
        list_res[[paste0("individual_metabolites_",i)]] <- temp_res3
        list_res[[paste0("MEs_metabolites_clusters_",i)]] <- temp_res4
      }
      ### return
      return(list_res)
    }
    
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      output = output[with(output,  order(cluster_name, -kIN)), ]
      write.table(output, file = "./results/Inter-Cor/hierarchical/individual_metabolites.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      res1 <- output
      
      ### Create a joint data frame of metabolite cluster eigengene
      ### equivalents/effective abundances('MEsMeta')
      MEsMeta = MEsMeta[rownames(MEsMeta),]
      
      ### Rename module names(columnames) from 'colors' to numbers
      colnames(MEsMeta) = cluster_mapping_file_m[paste0("Meta_", colnames(MEsMeta)), "New_Name"]
      
      ### Save MEsMeta to file, order by module number
      write.table(MEsMeta[,order(colnames(MEsMeta))], file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMeta[,order(colnames(MEsMeta))]
      list_res <- list(individual_metabolites = res1, MEs_metabolites_clusters = res2)
      ### return
      return(list_res)
    }
    
    else{
      stop("Something error in metabolites cluster!")
    }
  }
  
  
  #---------------------------------------------#
  #---------------------PCA---------------------#
  #---------------------------------------------#
  else if(dimReduMethod == "PCA"){
    
    #######################################################
    ###### Step 1 - Identify clusters of metabolites ######
    #######################################################
    
    meta.pr <- tryCatch(
      {
        ### normalization
        tmp_metabolites <- as.data.frame(scale(metabolites)) 
        ### PCA
        princomp(t(tmp_metabolites),cor=TRUE,scores = T) 
      },
      error = function(cond) {
        ### normalization
        tmp_metabolites <- as.data.frame(metabolites) 
        ### PCA
        meta.pr <- princomp(t(tmp_metabolites),cor=TRUE,scores = T) 
        return(meta.pr)
      }
    )    
    
    ### calculate the proportion of variance
    sd <- meta.pr$sdev
    var_proportion  <- sd^2/sum(sd^2)
    
    ### calculate the cumulative proportion
    cum_proportion <- 0
    top_num <- 0
    for (i in 1:length(var_proportion)) {
      cum_proportion <- cum_proportion + var_proportion[i]
      if(cum_proportion >= 0.9){
        if(i == 1){
          top_num <- i
        }
        else{
          top_num <- i - 1
        }
        break 
      }
    }
  
    ############ Exist continuous phenotype variables ############
    if(is.data.frame(phenotypes) & phenoDataType == "continuous"){
      
      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(meta.pr[["loadings"]][,1:top_num])
      colnames(MEsMeta) <- gsub("Comp.","Meta",temp_names)
      res1 <- MEsMeta
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      
      ### Save scores of metabolites in top n clusters 
      indi_metabolite <- as.data.frame(meta.pr$scores[,1:top_num]) 
      colnames(indi_metabolite) <- gsub("Comp.","Meta",colnames(indi_metabolite))
      
      ###########################################################################
      ###### Step 3 - Link individual metabolite to phenotypes of interest ######
      ###########################################################################

      tmpMat = array(NA, c(ncol(metabolites), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = colnames(metabolites)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating individual metabolites with phenotypes
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(metabolites, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      
      for (i in 1:length(colnames(phenotypes))){
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        indi_metabolite <- cbind(indi_metabolite, tmp_output)
      }
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      ### Write to file
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_metabolite
      
      ########################################################################
      ###### Step 4 - Link metabolite clusters to phenotype of interest ######
      ########################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_pheno <- list() ### Data structure for storing results of correlation tests under different setups
      tmpMat = array(NA, c(ncol(MEsMeta), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = names(MEsMeta)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating metabolite clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMeta, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_pheno[["meta"]] <- tmpMat
      rm(tmpMat)
      
      ######################################################
      ######## Step 5 - Save phenotype associations ######## 
      ######################################################
      
      ### Save the age association of metabolite clusters calculated in Step 4.
      tmpMat = cor_pheno[["meta"]]
      list_res <- list(individual_metabolites = res2,MEs_metabolites_clusters = res1)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_metabolites")]] <- tmp_output
      }
      ### return 
      return(list_res)
      
      # ### plot Bubble chart of scores 
      # col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
      #                            "cyan", "#007FFF", "blue", "#00007F"))
      # corrplot(meta.pr$scores[,1:top_num],is.corr = FALSE,method = "circle",tl.col = "black",tl.cex = 0.5,col = col1(100))
      
    }
    
    
    ############ Exist categorical phenotype variables ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############
      
      ### Save scores of metabolites in top n clusters 
      indi_metabolite <- as.data.frame(meta.pr$scores[,1:top_num]) 
      colnames(indi_metabolite) <- gsub("Comp.","Meta",colnames(indi_metabolite))
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- indi_metabolite
      
      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(meta.pr[["loadings"]][,1:top_num])
      colnames(MEsMeta) <- gsub("Comp.","Meta",temp_names)
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMeta
      list_res <- list(individual_metabolites = res1,MEs_metabolites_clusters = res2)
      
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(metabolites)  ### work in logarithmic space
        
        temp_metabolites_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        temp_metabolites = metabolites[temp_metabolites_rownames,]
        meta.pr <- tryCatch(
          {
            ### normalization
            tmp_metabolites <- as.data.frame(scale(temp_metabolites)) 
            ### PCA
            princomp(t(tmp_metabolites),cor=TRUE,scores = T) 
          },
          error = function(cond) {
            ### normalization
            tmp_metabolites <- as.data.frame(temp_metabolites) 
            ### PCA
            meta.pr <- princomp(t(tmp_metabolites),cor=TRUE,scores = T) 
            return(meta.pr)
          }
        )    
        
        ### calculate the proportion of variance
        sd <- meta.pr$sdev
        var_proportion  <- sd^2/sum(sd^2)
        
        ### calculate the cumulative proportion
        cum_proportion <- 0
        top_num <- 0
        for (j in 1:length(var_proportion)) {
          cum_proportion <- cum_proportion + var_proportion[j]
          if(cum_proportion >= 0.9){
            if(i == 1){
              top_num <- j
            }
            else{
              top_num <- j - 1
            }
            break 
          }
        }
        
        ### Save scores of metabolites in top n clusters 
        indi_metabolite <- as.data.frame(meta.pr$scores[,1:top_num]) 
        colnames(indi_metabolite) <- gsub("Comp.","Meta",colnames(indi_metabolite))
        
        ### cluster names of PCA result include all components
        cluster_name <- rep("all",nrow(indi_metabolite))
        cluster_name <- data.frame(cluster_name)
        indi_metabolite <- cbind(cluster_name,indi_metabolite)
        write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- indi_metabolite
        
        ### select top n clusters 
        MEsMeta <- as.data.frame(meta.pr[["loadings"]][,1:top_num]) 
        temp_names <- colnames(meta.pr[["loadings"]][,1:top_num])
        colnames(MEsMeta) <- gsub("Comp.","Meta",temp_names)
        write.table(MEsMeta, file = paste0("./results/Inter-Cor/hierarchical/MEs_metabolites_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMeta

        ### Save MEsMeta to file, order by module number
        list_res[[paste0("individual_metabolites_",i)]] <- temp_res3
        list_res[[paste0("MEs_metabolites_clusters_",i)]] <- temp_res4
      }
      return(list_res)
    }
    
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      ### Save scores of metabolites in top n clusters 
      indi_metabolite <- as.data.frame(meta.pr$scores[,1:top_num]) 
      colnames(indi_metabolite) <- gsub("Comp.","Meta",colnames(indi_metabolite))
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- indi_metabolite
      
      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(meta.pr[["loadings"]][,1:top_num])
      colnames(MEsMeta) <- gsub("Comp.","Meta",temp_names)
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMeta
      list_res <- list(individual_metabolites = res1,MEs_metabolites_clusters = res2)
      return(list_res)
    }
    
    else{
      stop("Something error in metabolites cluster!")
    }
  }
  
  
  #---------------------------------------------#
  #---------------------PCoA--------------------#
  #---------------------------------------------#
  else if(dimReduMethod == "PCoA"){
    #######################################################
    ###### Step 1 - Identify clusters of metabolites ######
    #######################################################
    
    ### PCoA
    meta.D <- vegdist(metabolites, "bray")
    meta.pcoa <- pcoa(meta.D,correction = "none")
    var_proportion <- meta.pcoa[["values"]][["Relative_eig"]]
    ### calculate the cumulative proportion
    cum_proportion <- 0
    top_num <- 0
    for (i in 1:length(var_proportion)) {
      cum_proportion <- cum_proportion + var_proportion[i]
      if(cum_proportion >= 0.9){
        if(i == 1){
          top_num <- i
        }
        else{
          top_num <- i - 1
        }
        break 
      }
    }
    
    ############ Exist continuous phenotype variables ############
    if(is.data.frame(phenotypes) & phenoDataType == "continuous"){
      
      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(meta.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMeta) <- gsub("Axis.","Meta",temp_names)
      res1 <- MEsMeta
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      
      ### Save scores of metabolites in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_metabolite <- as.data.frame(wascores(MEsMeta, metabolites))
      
      ###########################################################################
      ###### Step 3 - Link individual metabolite to phenotypes of interest ######
      ###########################################################################
      
      tmpMat = array(NA, c(ncol(metabolites), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = colnames(metabolites)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating individual metabolites with phenotypes
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(metabolites, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      
      for (i in 1:length(colnames(phenotypes))){
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        indi_metabolite <- cbind(indi_metabolite, tmp_output)
      }
      
      ### cluster names of PCoA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      ### Write to file
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_metabolite
      
      ########################################################################
      ###### Step 4 - Link metabolite clusters to phenotype of interest ######
      ########################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_pheno <- list() ### Data structure for storing results of correlation tests under different setups
      tmpMat = array(NA, c(ncol(MEsMeta), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = names(MEsMeta)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating metabolite clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMeta, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_pheno[["meta"]] <- tmpMat
      rm(tmpMat)
      
      ######################################################
      ######## Step 5 - Save phenotype associations ######## 
      ######################################################
      
      ### Save the age association of metabolite clusters calculated in Step 4.
      tmpMat = cor_pheno[["meta"]]
      list_res <- list(individual_metabolites = res2,MEs_metabolites_clusters = res1)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_metabolites")]] <- tmp_output
      }
      ### return 
      return(list_res)
      
      # ### plot Bubble chart of scores 
      # col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
      #                            "cyan", "#007FFF", "blue", "#00007F"))
      # corrplot(meta.pr$scores[,1:top_num],is.corr = FALSE,method = "circle",tl.col = "black",tl.cex = 0.5,col = col1(100))
      
    }
    
    
    ############ Exist categorical phenotype variables ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############

      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(meta.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMeta) <- gsub("Axis.","Meta",temp_names)
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- MEsMeta
      ### Save scores of metabolites in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_metabolite <- as.data.frame(wascores(MEsMeta, metabolites))
      
      ### cluster names of PCoA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      ### Write to file
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_metabolite
      list_res <- list(individual_metabolites = res2,MEs_metabolites_clusters = res1)
      
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(metabolites)  ### work in logarithmic space
        
        temp_metabolites_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        temp_metabolites = metabolites[temp_metabolites_rownames,]
        ### PCoA
        meta.D <- vegdist(temp_metabolites, "bray")
        meta.pcoa <- pcoa(meta.D,correction = "none")
        var_proportion <- meta.pcoa[["values"]][["Relative_eig"]]
        ### calculate the cumulative proportion
        cum_proportion <- 0
        top_num <- 0
        for (j in 1:length(var_proportion)) {
          cum_proportion <- cum_proportion + var_proportion[j]
          if(cum_proportion >= 0.9){
            if(j == 1){
              top_num <- j
            }
            else{
              top_num <- j - 1
            }
            break 
          }
        }
        
        ### select top n clusters 
        MEsMeta <- as.data.frame(meta.pcoa[["vectors"]][,1:top_num])
        temp_names <- colnames(meta.pcoa[["vectors"]])[1:top_num]
        colnames(MEsMeta) <- gsub("Axis.","Meta",temp_names)
        write.table(MEsMeta, file = paste0("./results/Inter-Cor/hierarchical/MEs_metabolites_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMeta
        ### Save scores of metabolites in top n clusters 
        ### wascores:Weighted Averages Scores for Species
        indi_metabolite <- as.data.frame(wascores(MEsMeta, temp_metabolites))
        
        ### cluster names of PCoA result include all components
        cluster_name <- rep("all",nrow(indi_metabolite))
        cluster_name <- data.frame(cluster_name)
        indi_metabolite <- cbind(cluster_name,indi_metabolite)
        
        ### Write to file
        write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- indi_metabolite

        ### Save MEsMeta to file, order by module number
        list_res[[paste0("individual_metabolites_",i)]] <- temp_res3
        list_res[[paste0("MEs_metabolites_clusters_",i)]] <- temp_res4
      }
      return(list_res)
  
    }
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      ### select top n clusters 
      MEsMeta <- as.data.frame(meta.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(meta.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMeta) <- gsub("Axis.","Meta",temp_names)
      write.table(MEsMeta, file = "./results/Inter-Cor/hierarchical/MEs_metabolites_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- MEsMeta
      ### Save scores of metabolites in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_metabolite <- as.data.frame(wascores(MEsMeta, metabolites))
      
      ### cluster names of PCoA result include all components
      cluster_name <- rep("all",nrow(indi_metabolite))
      cluster_name <- data.frame(cluster_name)
      indi_metabolite <- cbind(cluster_name,indi_metabolite)
      
      ### Write to file
      write.table(indi_metabolite, file = paste0("./results/Inter-Cor/hierarchical/individual_metabolites.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_metabolite
      list_res <- list(individual_metabolites = res2,MEs_metabolites_clusters = res1)
      return(list_res)
    }
    
    else{
      stop("Something error in metabolites cluster!")
    }
  }
}


