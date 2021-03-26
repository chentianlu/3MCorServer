# function - microbesCluster #
########################################################################
# File: microbesCluster.R
# Aim : Clustering for microbiome data
#---------------------------------------------------------------------------------------------------------------------
# Author : Tianlu Chen, Tao Sun, Dandan Liang, Mengci Li
# Email  : chentianlu@sjtu.edu.cn
# Date   : 2020-08
# Version: 1.0
#---------------------------------------------------------------------------------------------------------------------
#
#
######################################################################################
## Input:                                                                           ##
##    micData --- A dataframe of microbes                                           ##
##                       (rows: samples,columns: microbes)                          ## 
##    phenoData --- A dataframe of phenotype data                                   ##
##                (rows: samples,columns: phenotype)                                ##
##    phenoDataType --- Phenotype data type                                         ##
##    dimReduMethod --- Dimension reduction method                                  ##
##                      Default: WGCNA                                              ##
##    corMethod --- Correlation method of associating individual microbes or        ##
##                  microbiome cluster with phenotype                               ##
##                  Default: spearman                                               ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    list_res --- A list included the classified microbiome modules results        ##                             
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
microbesCluster <- function(micData,phenoData,phenoDataType,dimReduMethod = "WGCNA",corMethod = "spearman")
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
  
  microbes <- micData
  phenotypes <- phenoData
  
  #---------------------------------------------#
  #--------------------WGCNA--------------------#
  #---------------------------------------------#
  if(dimReduMethod == "WGCNA"){
    #################################################
    ###### Step 1 - Identify clusters microbes ######
    #################################################
    
    ###Specify data and parameters 
    # dat_tmp = log2 (microbes) ### work in logarithmic space
    dat_tmp = microbes
    
    ### Calculate weighted adjacency matrix
    network_adjacency = adjacency (dat_tmp, power = softPower, type = networkType, corFnc = corFuntion, corOptions = corOptionsStr)
    colnames (network_adjacency) = rownames (network_adjacency) = colnames (dat_tmp)
    
    ### Define dissimilarity based on topological overlap
    dissTOM = TOMdist (network_adjacency, TOMType = networkType)
    colnames (dissTOM) = rownames (dissTOM) = colnames (dat_tmp)
    
    ### Hierarchical clustering
    micTree = flashClust (as.dist (dissTOM), method = clusterMethod)
    
    ### Define modules by cutting branches
    moduleLabels1 = cutreeDynamic (dendro = micTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
    moduleLabels1 = labels2colors (moduleLabels1)
    
    ### Automatically merge highly correlated modules
    merge = mergeCloseModules (dat_tmp, moduleLabels1, corFnc = corFuntion, corOptions = corOptionsList, cutHeight = mergingThresh)
    
    ### Determine resulting merged module colors
    moduleLabels2 = merge$colors
    
    ### Establish eigengenes of the newly merged modules, used for cluster overall abundances
    MEs = merge$newMEs
    
    ### Choose final module assignments
    moduleColorsMic = moduleLabels2
    names (moduleColorsMic) = colnames (dat_tmp)
    MEsMic = orderMEs (MEs)
    rownames (MEsMic) = rownames (dat_tmp)
    
    ##Create the cluster_mapping_file file
    cluster_mapping_file_b <- data.frame(New_Name=paste0("Mic",1:length(MEsMic)),Description="None")
    row.names(cluster_mapping_file_b) <- paste0("Mic_",colnames(MEsMic))
    cluster_mapping_file_b$label =  sapply (rownames (cluster_mapping_file_b),  function (x) paste (cluster_mapping_file_b [x, "New_Name"], cluster_mapping_file_b [x, "Description"], sep = ": "))
    
    ### Determine relevant descriptive statistics of established clusters
    ### kIN: within-module connectivity, determined by summing connectivity with all
    ###      other microbes in the given cluster.
    ### kME: bicor-correlation between the microbes profile and module eigenvector; 
    ### both measures of intramodular hub-microbes status.
    kIN <-      vector (length = ncol (dat_tmp)); names (kIN) = colnames (dat_tmp)
    kME <-      vector (length = ncol (dat_tmp)); names (kME) = colnames (dat_tmp)
    modules <-  vector (length = ncol (dat_tmp)); names (modules) = colnames (dat_tmp)
    
    for (module in names (table (moduleColorsMic))) {   
      all.microbes = names (dat_tmp)
      inModule = (moduleColorsMic == module)
      module.microbes = names (moduleColorsMic [inModule])
      modules [module.microbes] = module 
      kIN [module.microbes] = sapply (module.microbes, function (x) sum (network_adjacency[x, module.microbes]) - 1)
      datKME = signedKME (dat_tmp, MEsMic, corFnc = corFuntion, corOptions = corOptionsStr)
      rownames (datKME) = colnames (dat_tmp)
      kME [module.microbes] = datKME [module.microbes, paste ("kME", module, sep = "")]   
    }
    output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN,"cluster_name" = sapply (modules, function (m) cluster_mapping_file_b [paste0 ("Mic_ME", m), "New_Name"]))
    
    
    ############ Exist continuous phenotypes variable ############
    if(is.data.frame(phenotypes) & phenoDataType == "continuous"){
      ################################################################################
      ########## Step 2 - Link individual microbes to phenotype of interest ##########
      ################################################################################
      
      tmpMat = array (NA, c (ncol (microbes), length(colnames(phenotypes)), 2))
      dimnames (tmpMat) [[1]] = colnames (microbes)
      dimnames (tmpMat) [[2]] = colnames(phenotypes)
      dimnames (tmpMat) [[3]] = c ("estimate", "p.value")
      
      ### Associating individual microbes with age
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(microbes, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        output <- cbind(output, tmp_output)
      }
      output = output[with(output,  order(cluster_name, -kIN)), ]
      
      ### Write to file
      write.table (output, file = "./results/Inter-Cor/hierarchical/individual_microbes.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      
      ### return 
      res1 <- output
      rm (output, tmpMat, dat_tmp)
      
      ########################################
      ###### step3 - Merge the results #######
      ########################################
      
      ### Create a joint data frame of microbes cluster eigengene
      ### equivalents/effective abundances ('MEsMic')
      MEsMic = MEsMic [rownames (MEsMic),]
      
      ### Rename module names (columnames) from 'colors' to numbers
      colnames (MEsMic) = cluster_mapping_file_b [paste0 ("Mic_", colnames (MEsMic)), "New_Name"]
      
      ### Save MEsMic to file, order by module number
      write.table (MEsMic[,order(colnames(MEsMic))], file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      
      ### return 
      res2 <- MEsMic[,order(colnames(MEsMic))]
      
      ################################################################################
      ########### Step 4 - Link microbes clusters to phenotype of interest ###########
      ################################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_age <- list () ### Data structure for storing results of correlation tests under different setups
      tmpMat = array (NA, c (ncol (MEsMic), length(colnames(phenotypes)), 2))
      dimnames (tmpMat) [[1]] = names (MEsMic)
      dimnames (tmpMat) [[2]] = colnames(phenotypes)
      dimnames (tmpMat) [[3]] = c ("estimate", "p.value")
      
      ### Associating microbes clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMic, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_age [["mic"]] <- tmpMat
      rm (tmpMat)
      
      #############################################################
      ########### Step 5 - Save phenotype associations ########### 
      #############################################################
      
      ### Save the age association of microbes clusters calculated in Step 4.
      tmpMat = cor_age [["mic"]]
      list_res <- list(individual_microbes = res1, MEs_microbes_clusters = res2)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_microbes")]] <- tmp_output
      }
      ### return 
      return(list_res)
    }
    
    
    ############ Exist categorical phenotype variables ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############
      
      output = output[with(output,  order(cluster_name, -kIN)), ]
      write.table(output, file = "./results/Inter-Cor/hierarchical/individual_microbes.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      res1 <- output
      
      ### Create a joint data frame of microbe cluster eigengene
      ### equivalents/effective abundances('MEsMic')
      MEsMic = MEsMic[rownames(MEsMic),]
      
      ### Rename module names(columnames) from 'colors' to numbers
      colnames(MEsMic) = cluster_mapping_file_b[paste0("Mic_", colnames(MEsMic)), "New_Name"]
      
      ### Save MEsMic to file, order by module number
      write.table(MEsMic[,order(colnames(MEsMic))], file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMic[,order(colnames(MEsMic))]
      
      list_res <- list(individual_microbes = res1, MEs_microbes_clusters = res2)
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(microbes)  ### work in logarithmic space
        
        temp_microbes_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        
        dat_tmp = microbes[temp_microbes_rownames,]
        
        ### Calculate weighted adjacency matrix
        network_adjacency = adjacency(dat_tmp, power = softPower, type = networkType, corFnc = corFuntion, corOptions = corOptionsStr)
        colnames(network_adjacency) = rownames(network_adjacency) = colnames(dat_tmp)
        
        ### Define dissimilarity based on topological overlap
        dissTOM = TOMdist(network_adjacency, TOMType = networkType)
        colnames(dissTOM) = rownames(dissTOM) = colnames(dat_tmp)
        
        ### Hierarchical clustering
        micTree = flashClust(as.dist(dissTOM), method = clusterMethod)
        
        ### Define modules by cutting branches
        moduleLabels1 = cutreeDynamic(dendro = micTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
        moduleLabels1 = labels2colors(moduleLabels1)
        
        ### Automatically merge highly correlated modules
        merge = mergeCloseModules(dat_tmp, moduleLabels1, corFnc = corFuntion, corOptions = corOptionsList, cutHeight = mergingThresh)
        
        ### Determine resulting merged module colors
        moduleLabels2 = merge$colors
        
        ### Establish eigengenes of the newly merged modules, used for cluster overall abundances
        MEs = merge$newMEs
        
        ### Choose final module assignments
        moduleColorsMic = moduleLabels2
        names(moduleColorsMic) = colnames(dat_tmp)
        MEsMic = orderMEs(MEs)
        rownames(MEsMic) = rownames(dat_tmp)
        
        ##Create the cluster_mapping_file file
        cluster_mapping_file_b <- data.frame(New_Name=paste0("Mic",1:length(MEsMic)),Description=c("None"))
        row.names(cluster_mapping_file_b) <- paste0("Mic_",colnames(MEsMic))
        cluster_mapping_file_b$label = sapply(rownames(cluster_mapping_file_b), function(x) paste(cluster_mapping_file_b[x, "New_Name"], cluster_mapping_file_b[x, "Description"], sep = ": "))
        
        ### Determine relevant descriptive statistics of established clusters
        ### kIN: within-module connectivity, determined by summing connectivity with all
        ###      other microbes in the given cluster.
        ### kME: bicor-correlation between the microbe profile and module eigenvector; 
        ### both measures of intramodular hub-microbe status.
        kIN <-      vector(length = ncol(dat_tmp)); names(kIN) = colnames(dat_tmp)
        kME <-      vector(length = ncol(dat_tmp)); names(kME) = colnames(dat_tmp)
        modules <-  vector(length = ncol(dat_tmp)); names(modules) = colnames(dat_tmp)
        
        for(module in names(table(moduleColorsMic))) {   
          all.microbes = names(dat_tmp)
          inModule =(moduleColorsMic == module)
          module.microbes = names(moduleColorsMic[inModule])
          modules[module.microbes] = module 
          kIN[module.microbes] = sapply(module.microbes, function(x) sum(network_adjacency[x, module.microbes]) - 1)
          datKME = signedKME(dat_tmp, MEsMic, corFnc = corFuntion, corOptions = corOptionsStr)
          rownames(datKME) = colnames(dat_tmp)
          kME[module.microbes] = datKME[module.microbes, paste("kME", module, sep = "")]   
        }
        temp_output = data.frame("module" = modules, "kME" = kME, "kIN" = kIN,"cluster_name" = sapply(modules, function(m) cluster_mapping_file_b[paste0("Mic_ME", m), "New_Name"]))
        temp_output = temp_output[with(temp_output,  order(cluster_name, -kIN)), ]
        write.table(temp_output, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- temp_output
        
        ### Create a joint data frame of microbe cluster eigengene
        ### equivalents/effective abundances('MEsMic')
        MEsMic = MEsMic[rownames(MEsMic),]
        
        ### Rename module names(columnames) from 'colors' to numbers
        colnames(MEsMic) = cluster_mapping_file_b[paste0("Mic_", colnames(MEsMic)), "New_Name"]
        
        ### Save MEsMic to file, order by module number
        write.table(MEsMic[,order(colnames(MEsMic))], file = paste0("./results/Inter-Cor/hierarchical/MEs_microbes_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMic[,order(colnames(MEsMic))]
        list_res[[paste0("individual_microbes_",i)]] <- temp_res3
        list_res[[paste0("MEs_microbes_clusters_",i)]] <- temp_res4
      }
      ### return
      return(list_res)
    }
    
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      output = output[with(output,  order(cluster_name, -kIN)), ]
      ### Write to file
      write.table (output, file = "./results/Inter-Cor/hierarchical/individual_microbes.txt", sep = "\t", col.names = NA, quote = F, row.names = T)
      res1 <- output
      
      ### Create a joint data frame of microbes cluster eigengene
      ### equivalents/effective abundances ('MEsMic')
      MEsMic = MEsMic [rownames (MEsMic),]
      
      ### Rename module names (columnames) from 'colors' to numbers
      colnames (MEsMic) = cluster_mapping_file_b [paste0 ("Mic_", colnames (MEsMic)), "New_Name"]
      
      ### Save MEsMic to file, order by module number
      write.table (MEsMic[,order(colnames(MEsMic))], file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMic[,order(colnames(MEsMic))]
      list_res <- list(individual_microbes = res1, MEs_microbes_clusters = res2)
      ### return
      return(list_res)
    }
    
    else{
      stop("Something error in microbes cluster!")
    }
  }
  
  
  #---------------------------------------------#
  #---------------------PCA---------------------#
  #---------------------------------------------#
  else if(dimReduMethod == "PCA"){
    
    #######################################################
    ###### Step 1 - Identify clusters of microbes ######
    #######################################################
    
    mic.pr <- tryCatch(
      {
        ### normalization
        tmp_microbes <- as.data.frame(scale(microbes)) 
        ### PCA
        princomp(t(tmp_microbes),cor=TRUE,scores = T) 
      },
      error = function(cond) {
        ### normalization
        tmp_microbes <- as.data.frame(microbes) 
        ### PCA
        mic.pr <- princomp(t(tmp_microbes),cor=F,scores = T) 
        return(mic.pr)
      }
    )    
    
    ### calculate the proportion of variance
    sd <- mic.pr$sdev
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
      MEsMic <- as.data.frame(mic.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(mic.pr[["loadings"]][,1:top_num])
      colnames(MEsMic) <- gsub("Comp.","Mic",temp_names)
      res1 <- MEsMic
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      
      ### Save scores of microbes in top n clusters 
      indi_microbe <- as.data.frame(mic.pr$scores[,1:top_num]) 
      colnames(indi_microbe) <- gsub("Comp.","Mic",colnames(indi_microbe))
      
      ###########################################################################
      ###### Step 3 - Link individual microbe to phenotypes of interest ######
      ###########################################################################
      
      tmpMat = array(NA, c(ncol(microbes), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = colnames(microbes)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating individual microbes with phenotypes
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(microbes, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      
      for (i in 1:length(colnames(phenotypes))){
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        indi_microbe <- cbind(indi_microbe, tmp_output)
      }
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      ### Write to file
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_microbe
      
      ########################################################################
      ###### Step 4 - Link microbe clusters to phenotype of interest ######
      ########################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_pheno <- list() ### Data structure for storing results of correlation tests under different setups
      tmpMat = array(NA, c(ncol(MEsMic), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = names(MEsMic)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating microbe clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMic, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_pheno[["mic"]] <- tmpMat
      rm(tmpMat)
      
      ######################################################
      ######## Step 5 - Save phenotype associations ######## 
      ######################################################
      
      ### Save the age association of microbe clusters calculated in Step 4.
      tmpMat = cor_pheno[["mic"]]
      list_res <- list(individual_microbes = res2,MEs_microbes_clusters = res1)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_microbes")]] <- tmp_output
      }
      ### return 
      return(list_res)
      
      # ### plot Bubble chart of scores 
      # col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
      #                            "cyan", "#007FFF", "blue", "#00007F"))
      # corrplot(mic.pr$scores[,1:top_num],is.corr = FALSE,method = "circle",tl.col = "black",tl.cex = 0.5,col = col1(100))
      
    }
    
    
    ############ Exist categorical phenotypes variable ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############
      
      ### Save scores of microbes in top n clusters 
      indi_microbe <- as.data.frame(mic.pr$scores[,1:top_num]) 
      colnames(indi_microbe) <- gsub("Comp.","Mic",colnames(indi_microbe))
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- indi_microbe
      
      ### select top n clusters 
      MEsMic <- as.data.frame(mic.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(mic.pr[["loadings"]][,1:top_num])
      colnames(MEsMic) <- gsub("Comp.","Mic",temp_names)
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMic
      list_res <- list(individual_microbes = res1,MEs_microbes_clusters = res2)
      
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(microbes)  ### work in logarithmic space
        
        temp_microbes_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        temp_microbes = microbes[temp_microbes_rownames,]
        mic.pr <- tryCatch(
          {
            ### normalization
            tmp_microbes <- as.data.frame(scale(temp_microbes)) 
            ### PCA
            princomp(t(tmp_microbes),cor=TRUE,scores = T) 
          },
          error = function(cond) {
            ### normalization
            tmp_microbes <- as.data.frame(temp_microbes) 
            ### PCA
            mic.pr <- princomp(t(tmp_microbes),cor=TRUE,scores = T) 
            return(mic.pr)
          }
        )    
        
        ### calculate the proportion of variance
        sd <- mic.pr$sdev
        var_proportion  <- sd^2/sum(sd^2)
        
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
        
        ### Save scores of microbes in top n clusters 
        indi_microbe <- as.data.frame(mic.pr$scores[,1:top_num]) 
        colnames(indi_microbe) <- gsub("Comp.","Mic",colnames(indi_microbe))
        
        ### cluster names of PCA result include all components
        cluster_name <- rep("all",nrow(indi_microbe))
        cluster_name <- data.frame(cluster_name)
        indi_microbe <- cbind(cluster_name,indi_microbe)
        write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- indi_microbe
        
        ### select top n clusters 
        MEsMic <- as.data.frame(mic.pr[["loadings"]][,1:top_num]) 
        temp_names <- colnames(mic.pr[["loadings"]][,1:top_num])
        colnames(MEsMic) <- gsub("Comp.","Mic",temp_names)
        write.table(MEsMic, file = paste0("./results/Inter-Cor/hierarchical/MEs_microbes_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMic
        
        ### Save MEsMic to file, order by module number
        list_res[[paste0("individual_microbes_",i)]] <- temp_res3
        list_res[[paste0("MEs_microbes_clusters_",i)]] <- temp_res4
      }
      return(list_res)
    }
    
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      ### Save scores of microbes in top n clusters 
      indi_microbe <- as.data.frame(mic.pr$scores[,1:top_num]) 
      colnames(indi_microbe) <- gsub("Comp.","Mic",colnames(indi_microbe))
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- indi_microbe
      
      ### select top n clusters 
      MEsMic <- as.data.frame(mic.pr[["loadings"]][,1:top_num]) 
      temp_names <- colnames(mic.pr[["loadings"]][,1:top_num])
      colnames(MEsMic) <- gsub("Comp.","Mic",temp_names)
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- MEsMic
      list_res <- list(individual_microbes = res1,MEs_microbes_clusters = res2)
      return(list_res)
    }
    
    else{
      stop("Something error in microbes cluster!")
    }
  }
  
  
  #---------------------------------------------#
  #---------------------PCoA--------------------#
  #---------------------------------------------#
  else if(dimReduMethod == "PCoA"){
    #######################################################
    ###### Step 1 - Identify clusters of microbes ######
    #######################################################
    
    ### PCoA
    mic.D <- vegdist(microbes, "bray")
    mic.pcoa <- pcoa(mic.D,correction = "none")
    var_proportion <- mic.pcoa[["values"]][["Relative_eig"]]
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
      MEsMic <- as.data.frame(mic.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(mic.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMic) <- gsub("Axis.","Mic",temp_names)
      res1 <- MEsMic
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      ### Save scores of microbes in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_microbe <- as.data.frame(wascores(MEsMic, microbes))
      
      ###########################################################################
      ###### Step 3 - Link individual microbe to phenotypes of interest ######
      ###########################################################################
      
      tmpMat = array(NA, c(ncol(microbes), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = colnames(microbes)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating individual microbes with phenotypes
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(microbes, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      
      for (i in 1:length(colnames(phenotypes))){
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        indi_microbe <- cbind(indi_microbe, tmp_output)
      }
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      ### Write to file
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_microbe
      
      ########################################################################
      ###### Step 4 - Link microbe clusters to phenotype of interest ######
      ########################################################################
      
      ### This is a core analysis step generating associations between the
      ### integrated/clustered ¨Comics data and a clinically interesting phenotype.
      
      cor_pheno <- list() ### Data structure for storing results of correlation tests under different setups
      tmpMat = array(NA, c(ncol(MEsMic), length(colnames(phenotypes)), 2))
      dimnames(tmpMat)[[1]] = names(MEsMic)
      dimnames(tmpMat)[[2]] = colnames(phenotypes)
      dimnames(tmpMat)[[3]] = c("estimate", "p.value")
      
      ### Associating microbe clusters with age     
      for(i in 1:length(colnames(phenotypes))) {
        tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")] =
          t(apply(MEsMic, MARGIN = 2, FUN = function(x)
            unlist(cor.test(phenotypes[,i], x,
                            method = corMethod, use = "pairwise.complete.obs")[c("estimate", "p.value")])))
      }
      cor_pheno[["mic"]] <- tmpMat
      rm(tmpMat)
      
      ######################################################
      ######## Step 5 - Save phenotype associations ######## 
      ######################################################
      
      ### Save the age association of microbe clusters calculated in Step 4.
      tmpMat = cor_pheno[["mic"]]
      list_res <- list(individual_microbes = res2,MEs_microbes_clusters = res1)
      for (i in 1:length(colnames(phenotypes))) {
        tmp_output = cbind(tmpMat[, colnames(phenotypes)[i], c("estimate", "p.value")], "p.adjust" = p.adjust(tmpMat[, colnames(phenotypes)[i], "p.value"], method = "BH"))
        colnames(tmp_output) = paste(rep(colnames(phenotypes)[i], 3), colnames(tmp_output), sep = "_")
        write.table(tmp_output, file = paste0("./results/Inter-Cor/hierarchical/",colnames(phenotypes)[i],"_cluster_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        list_res[[paste0(colnames(phenotypes)[i],"_cluster_microbes")]] <- tmp_output
      }
      ### return 
      return(list_res)
      
      # ### plot Bubble chart of scores 
      # col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
      #                            "cyan", "#007FFF", "blue", "#00007F"))
      # corrplot(mic.pr$scores[,1:top_num],is.corr = FALSE,method = "circle",tl.col = "black",tl.cex = 0.5,col = col1(100))
      
    }
    
    
    ############ Exist categorical phenotype variables ############
    else if(is.data.frame(phenotypes) & phenoDataType == "categorical"){
      ############ All groups ############
      
      ### select top n clusters 
      MEsMic <- as.data.frame(mic.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(mic.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMic) <- gsub("Axis.","Mic",temp_names)
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- MEsMic
      ### Save scores of microbes in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_microbe <- as.data.frame(wascores(MEsMic, microbes))
      
      ### cluster names of PCoA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      ### Write to file
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_microbe
      list_res <- list(individual_microbes = res2,MEs_microbes_clusters = res1)
      
      ############ Each groups ############
      temp_list_res <- list()
      for (i in unique(phenoData[,1])) {
        ###Specify data and parameters
        # dat_tmp = log2(microbes)  ### work in logarithmic space
        
        temp_microbes_rownames <- rownames(phenoData)[which(phenoData[,1] == i)]  
        temp_microbes = microbes[temp_microbes_rownames,]
        ### PCoA
        mic.D <- vegdist(temp_microbes, "bray")
        mic.pcoa <- pcoa(mic.D,correction = "none")
        var_proportion <- mic.pcoa[["values"]][["Relative_eig"]]
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
        MEsMic <- as.data.frame(mic.pcoa[["vectors"]][,1:top_num])
        temp_names <- colnames(mic.pcoa[["vectors"]])[1:top_num]
        colnames(MEsMic) <- gsub("Axis.","Mic",temp_names)
        write.table(MEsMic, file = paste0("./results/Inter-Cor/hierarchical/MEs_microbes_clusters_",i,".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
        temp_res4 <- MEsMic
        ### Save scores of microbes in top n clusters 
        ### wascores:Weighted Averages Scores for Species
        indi_microbe <- as.data.frame(wascores(MEsMic, temp_microbes))
        
        ### cluster names of PCoA result include all components
        cluster_name <- rep("all",nrow(indi_microbe))
        cluster_name <- data.frame(cluster_name)
        indi_microbe <- cbind(cluster_name,indi_microbe)
        
        ### Write to file
        write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes_",i,".txt"), sep = "\t", col.names = NA, quote = F, row.names = T)
        temp_res3 <- indi_microbe
        
        ### Save MEsMic to file, order by module number
        list_res[[paste0("individual_microbes_",i)]] <- temp_res3
        list_res[[paste0("MEs_microbes_clusters_",i)]] <- temp_res4
      }
      return(list_res)
    }
    
    
    ############ No phenotypes ############
    else if(!is.data.frame(phenotypes)){
      ### select top n clusters 
      MEsMic <- as.data.frame(mic.pcoa[["vectors"]][,1:top_num])
      temp_names <- colnames(mic.pcoa[["vectors"]])[1:top_num]
      colnames(MEsMic) <- gsub("Axis.","Mic",temp_names)
      write.table(MEsMic, file = "./results/Inter-Cor/hierarchical/MEs_microbes_clusters.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
      res1 <- MEsMic
      ### Save scores of microbes in top n clusters 
      ### wascores:Weighted Averages Scores for Species
      indi_microbe <- as.data.frame(wascores(MEsMic, microbes))
      
      ### cluster names of PCA result include all components
      cluster_name <- rep("all",nrow(indi_microbe))
      cluster_name <- data.frame(cluster_name)
      indi_microbe <- cbind(cluster_name,indi_microbe)
      
      ### Write to file
      write.table(indi_microbe, file = paste0("./results/Inter-Cor/hierarchical/individual_microbes.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
      res2 <- indi_microbe
      list_res <- list(individual_microbes = res2,MEs_microbes_clusters = res1)
      return(list_res)
    }
    
    else{
      stop("Something error in microbes cluster!")
    }
  }
}


