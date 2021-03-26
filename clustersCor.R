# function - clustersCor #
########################################################################
# File: clustersCor.R
# Aim : The correlation between metabolome modules and microbiome modules
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
##		metaCluData --- A list included the classified metabolome modules results from  ##
##                   metabolitesCluster.R                                             ##
##		micCluData --- A list included the classified microbiome modules results from   ##
##                   microbesCluster.R                                                ##
##    phenoData --- A dataframe of phenotype data                                     ##
##                  (rows: samples,columns: phenotype)                                ##
##    phenotype --- Phenotype name for filtering modules by correlation analysis      ##
##                  between modules and specific phenotype data                       ##
##    phenoDataType --- Phenotype data type                                           ##
##                Default: continuous                                                 ##
##    confounderData --- A dataframe of confounder                                    ##
##                       (rows:samples,columns:confounder)                            ##
##                Default: NA                                                         ##
##    fdrThreshold --- FDR threshold for filtering modules by correlation analysis    ##
##                  between modules and specific phenotype data                       ##
##                  Default: 0.1                                                      ##
##    corMethod --- Method of the correlation between metabolome modules and          ##
##                  microbiome modules                                                ##
##                  Default: spearman                                                 ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
## 	  cor_res ---  A list included coefficient, p value and p.adjust value            ##
##                 of the correlation between metabolome modules and microbiome       ##
##                 modules                                                            ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
clustersCor <- function(metaCluData,micCluData,phenoData,phenotype,phenoDataType,
                       confounderData = NA,fdrThreshold = 0.1,corMethod = "spearman")
{
  if(!file.exists("./results/Inter-Cor/hierarchical")){
    dir.create("./results/Inter-Cor/hierarchical",recursive = T)
  }
  
  pheno2metaCluData <- as.data.frame(metaCluData[[paste0(phenotype,"_cluster_metabolites")]])
  pheno2micCluData <- as.data.frame(micCluData[[paste0(phenotype,"_cluster_microbes")]])

  ############ Exist continuous phenotypes variable############
  if(!is.na(phenoDataType) & phenoDataType == "continuous"){
    ### Identify the phenotype
    temp_phenotype <- substring(colnames(pheno2metaCluData)[1], 1, nchar(colnames(pheno2metaCluData)[1])-9) 
    if(is.na(temp_phenotype)){
      print(paste0("ERROR:",phenotype," does not exist!"))
    }
    else{
      ################################################################
      #### Step 1 - Select features with significant differences #####
      ################################################################
      
      ### Here, combine and integrate those functional, taxonomic and metabolomics
      ### features which reliably correspond to the host phenotype of interest, then
      ### later determine their inter-correlations.
      
      ### Extract metabolites clusters with significant characteristics 
      sig_name_metabolites <- row.names(pheno2metaCluData)[pheno2metaCluData[,3] < fdrThreshold]
      MEsmeta <- metaCluData
      # MEsmeta <- MEsmeta[["MEs_metabolites_clusters"]]
      MEsmeta <- as.data.frame(MEsmeta[["MEs_metabolites_clusters"]][,sig_name_metabolites])
      colnames(MEsmeta) <- sig_name_metabolites
      rownames(MEsmeta) <- rownames(metaCluData[["MEs_metabolites_clusters"]])
      ### Extract microbes clusters with significant characteristics 
      sig_name_microbes <- row.names(pheno2micCluData)[pheno2micCluData[,3] < fdrThreshold]
      MEsmic <- micCluData
      # MEsmic <- MEsmic[["MEs_microbes_clusters"]]
      MEsmic <- as.data.frame(MEsmic[["MEs_microbes_clusters"]][,sig_name_microbes])
      colnames(MEsmic) <- sig_name_microbes
      rownames(MEsmic) <- rownames(micCluData[["MEs_microbes_clusters"]])
      
      #####################################################################################
      ############ Step 2 - Correlate metabolite clusters to microbes clusters ############
      #####################################################################################
      
      ### The set of metabolite clusters associated with the phenotype of interest
      ### should next be tested for association with the set of microbes
      ### clusters likewise so associated (identified in Step 1). 
      ###
      ### In this step, the correlations between each metabolites cluster and
      ### microbes cluster are determined.
      
      
      ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
      if(corMethod == "gramm"){
        naivegramm <- naivegramm(metaData = MEsmeta,micData = MEsmic,covdata = confounderData,pheno = "All")
        write.csv(naivegramm[["r_All"]],paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,".csv"))
        write.csv(naivegramm[["p_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,".csv"))
        if(nrow(naivegramm[["p_All"]]) == 1 & ncol(naivegramm[["p_All"]]) != 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,],method = "BH")
          mic_meta_p.adjust <- as.data.frame(t(mic_meta_p.adjust))
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else if(nrow(naivegramm[["p_All"]]) != 1 & ncol(naivegramm[["p_All"]]) == 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][,1],method = "BH")
          mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else if(nrow(naivegramm[["p_All"]]) == 1 & ncol(naivegramm[["p_All"]]) == 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,1],method = "BH")
          mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else{
          mic_meta_p.adjust = apply(naivegramm[["p_All"]], 2, FUN = function(x){
            p.adjust(x,method = "BH")})
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }
        write.csv(naivegramm[["r_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,".csv"))
        write.csv(naivegramm[["p_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,".csv"))
        write.csv(naivegramm[["p.adjust_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,".csv"))
        write.csv(naivegramm[["type_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_type_by_",phenotype,".csv"))
        gramm_res <- naivegramm
        return(gramm_res)
      }
      
      ### Method 2 - Spearman,Pearson,Kendall
      if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
          mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }

        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,".csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        return(cor_res)
      }
      
      ### Method 3 - Partial correlation
      if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
        corMethod = strsplit(corMethod," ")[[1]][2]
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          # exist confounder
          if(is.na(confounderData)){
            stop("Lack of confounders!")
          }else{
            mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
              pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$estimate)
            mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
              pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$p.value)
          }
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }
        
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,".csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        return(cor_res)
      }
      
      ### Method 4 - Generalized linear models
      if (corMethod == "glm") {
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          # exist confounder
          if(is.na(confounderData)){
            mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,4])
          }else{
            temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
            mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
          }
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }
        
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,".csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        return(cor_res)
      }
    }
    
  }
  
  
  ############ Exist categorical phenotypes variable############
  else if(!is.na(phenoDataType) & phenoDataType == "categorical"){
    phenoData$SampleID <- rownames(phenoData)
    
    ### Identify the phenotype
    temp_phenotype <- colnames(phenoData)
    if(!(phenotype %in% temp_phenotype)){
      print(paste0("ERROR:",phenotype," does not exist!"))
    }
    else{
      #############################################################################################################
      #### Step 1 - Select features with significant differences by Mann-Whitney U test or Kruskal-Wallis test ####
      #############################################################################################################
      
      ### Identify phenotype groups
      temp <- cbind(metaCluData[["MEs_metabolites_clusters"]],rownames(metaCluData[["MEs_metabolites_clusters"]]))
      colnames(temp)[ncol(temp)] <- "SampleID"
      test_df <- merge(temp,phenoData,all.x = T)
      
      ### Test with each metabolites cluster between groups
      colnames(test_df)[ncol(test_df)] <- "group"
      sig_name_metabolites <- c()
      
      ### Mann-Whitney U test for two groups
      if(length(unique(test_df$group)) == 2){
        rownames1 <- rownames(test_df)[which(test_df$group == unique(test_df$group)[1])]
        group1 <- test_df[rownames1,]
        
        rownames2 <- rownames(test_df)[which(test_df$group == unique(test_df$group)[2])]
        group2 <- test_df[rownames2,]
        
        for (i in 3:(ncol(test_df))-1) {
          module1 <- group1[,i]
          module2 <- group2[,i]
          wt_res <- wilcox.test(module2,module1)
          if(wt_res$p.value < fdrThreshold){
            sig_name_metabolites <- append(sig_name_metabolites,colnames(test_df)[i])
          }
        }
      }else if(length(unique(test_df$group)) > 2){     ### Kruskal-Wallis test for three or more groups
        test_df$group <- as.factor(test_df$group)
        for (i in 3:(ncol(test_df))-1) {
          kt_res <- kruskal.test(test_df[,i] ~ group, data = test_df) 
          if(kt_res$p.value < fdrThreshold){
            sig_name_metabolites <- append(sig_name_metabolites,colnames(test_df)[i])
          }
        }
      }else{
        stop("Error in group counts!")
      }
      
      ### Extract metabolites clusters with significant characteristics 
      MEsmeta <- metaCluData
      # MEsmeta <- MEsmeta[["MEs_metabolites_clusters"]]
      MEsmeta <- as.data.frame(MEsmeta[["MEs_metabolites_clusters"]][,sig_name_metabolites])
      colnames(MEsmeta) <- sig_name_metabolites
      rownames(MEsmeta) <- rownames(metaCluData[["MEs_metabolites_clusters"]])
      
      ### Identify phenotype groups
      temp <- cbind(micCluData[["MEs_microbes_clusters"]],rownames(micCluData[["MEs_microbes_clusters"]]))
      colnames(temp)[ncol(temp)] <- "SampleID"
      test_df <- merge(temp,phenoData,all.x = T)
      
      ### Test with each microbes cluster between groups
      colnames(test_df)[ncol(test_df)] <- "group"
      sig_name_microbes <- c()
      
      ### Mann-Whitney U test for two groups
      if(length(unique(test_df$group)) == 2){
        rownames1 <- rownames(test_df)[which(test_df$group == unique(test_df$group)[1])]
        group1 <- test_df[rownames1,]
        
        rownames2 <- rownames(test_df)[which(test_df$group == unique(test_df$group)[2])]
        group2 <- test_df[rownames2,]
        
        for (i in 3:(ncol(test_df))-1) {
          module1 <- group1[,i]
          module2 <- group2[,i]
          wt_res <- wilcox.test(module2,module1)
          if(wt_res$p.value < fdrThreshold){
            sig_name_microbes <- append(sig_name_microbes,colnames(test_df)[i])
          }
        }
      }else if(length(unique(test_df$group)) > 2){  ### Kruskal-Wallis test for three or more groups
        test_df$group <- as.factor(test_df$group)
        for (i in 3:(ncol(test_df))-1) {
          kt_res <- kruskal.test(test_df[,i] ~ group, data = test_df) 
          if(kt_res$p.value < fdrThreshold){
            sig_name_microbes <- append(sig_name_microbes,colnames(test_df)[i])
          }
        }
      }else{
        stop("Error in group counts!")
      }
      
      ### Extract microbes clusters with significant characteristics 
      MEsmic <- micCluData
      # MEsmic <- MEsmic[["MEs_microbes_clusters"]]
      MEsmic <- as.data.frame(MEsmic[["MEs_microbes_clusters"]][,sig_name_microbes])
      colnames(MEsmic) <- sig_name_microbes
      rownames(MEsmic) <- rownames(micCluData[["MEs_microbes_clusters"]])
      
      
      #####################################################################################
      ############ Step 2 - Correlate metabolite clusters to microbes clusters ############
      #####################################################################################
      
      ### The set of metabolite clusters associated with the phenotype of interest
      ### should next be tested for association with the set of microbes
      ### clusters likewise so associated (identified in Step 1). 
      ###
      ### In this step, the correlations between each metabolites cluster and
      ### microbes cluster are determined.
      
      
      ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
      if(corMethod == "gramm"){
        ### All groups
        naivegramm <- naivegramm(MEsmeta,MEsmic,covdata = confounderData,pheno = "All")
        if(nrow(naivegramm[["p_All"]]) == 1 & ncol(naivegramm[["p_All"]]) != 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,],method = "BH")
          mic_meta_p.adjust <- as.data.frame(t(mic_meta_p.adjust))
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else if(nrow(naivegramm[["p_All"]]) != 1 & ncol(naivegramm[["p_All"]]) == 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][,1],method = "BH")
          mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else if(nrow(naivegramm[["p_All"]]) == 1 & ncol(naivegramm[["p_All"]]) == 1){
          mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,1],method = "BH")
          mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
          colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }else{
          mic_meta_p.adjust = apply(naivegramm[["p_All"]], 2, FUN = function(x){
            p.adjust(x,method = "BH")})
          naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
        }
        write.csv(naivegramm[["r_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_all.csv"))
        write.csv(naivegramm[["p_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_all.csv"))
        write.csv(naivegramm[["p.adjust_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_all.csv"))
        write.csv(naivegramm[["type_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_type_by_",phenotype,"_all.csv"))
        
        ### Each groups
        for (i in unique(phenoData[,1])) {
          temp_MEsmeta <- metaCluData[[paste0("MEs_metabolites_clusters_",i)]]
          temp_MEsmic <- micCluData[[paste0("MEs_microbes_clusters_",i)]]
          if(is.na(confounderData)){
            temp_confounderData <- NA
          }else{
            temp_confounderData <- as.data.frame(confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]) 
            rownames(temp_confounderData) <- rownames(temp_MEsmeta)
            colnames(temp_confounderData) <- colnames(confounderData)
          }
          temp_naivegramm <- naivegramm(temp_MEsmeta,temp_MEsmic,covdata = temp_confounderData,pheno = i)
          
          if(nrow(temp_naivegramm[[paste0("p_",i)]]) == 1 & ncol(temp_naivegramm[[paste0("p_",i)]]) != 1){
            mic_meta_p.adjust = p.adjust(temp_naivegramm[[paste0("p_",i)]][1,],method = "BH")
            mic_meta_p.adjust <- as.data.frame(t(mic_meta_p.adjust))
            rownames(mic_meta_p.adjust) <- rownames(temp_naivegramm[[paste0("p_",i)]])
            colnames(mic_meta_p.adjust) <- colnames(temp_naivegramm[[paste0("p_",i)]])
            temp_naivegramm[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
          }else if(nrow(temp_naivegramm[[paste0("p_",i)]]) != 1 & ncol(temp_naivegramm[[paste0("p_",i)]]) == 1){
            mic_meta_p.adjust = p.adjust(temp_naivegramm[[paste0("p_",i)]][,1],method = "BH")
            mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_p.adjust) <- rownames(temp_naivegramm[[paste0("p_",i)]])
            colnames(mic_meta_p.adjust) <- colnames(temp_naivegramm[[paste0("p_",i)]])
            temp_naivegramm[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
          }else if(nrow(temp_naivegramm[[paste0("p_",i)]]) == 1 & ncol(temp_naivegramm[[paste0("p_",i)]]) == 1){
            mic_meta_p.adjust = p.adjust(temp_naivegramm[[paste0("p_",i)]][1,1],method = "BH")
            mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_p.adjust) <- rownames(temp_naivegramm[[paste0("p_",i)]])
            colnames(mic_meta_p.adjust) <- colnames(temp_naivegramm[[paste0("p_",i)]])
            temp_naivegramm[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
          }else{
            mic_meta_p.adjust = apply(temp_naivegramm[[paste0("p_",i)]], 2, FUN = function(x){
              p.adjust(x,method = "BH")})
            temp_naivegramm[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
          }
          write.csv(temp_naivegramm[[paste0("r_",i)]],paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_",i,".csv"))
          write.csv(temp_naivegramm[[paste0("p_",i)]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_",i,".csv"))
          write.csv(temp_naivegramm[[paste0("p.adjust_",i)]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_",i,".csv"))
          write.csv(temp_naivegramm[[paste0("type_",i)]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_type_by_",phenotype,"_",i,".csv"))
          naivegramm[[paste0("r_",i)]] <- temp_naivegramm[[paste0("r_",i)]]
          naivegramm[[paste0("p_",i)]] <- temp_naivegramm[[paste0("p_",i)]]
          naivegramm[[paste0("p.adjust_",i)]] <- temp_naivegramm[[paste0("p.adjust_",i)]]
          naivegramm[[paste0("type_",i)]] <- temp_naivegramm[[paste0("type_",i)]]
        }
        gramm_res <- naivegramm
        return(gramm_res)
      }
      
      
      ### Method 2 - Method 2 - Spearman,Pearson,Kendall
      if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
        ### All groups
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
          mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }
        
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_all.csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        ### Each groups
        for (i in unique(phenoData[,1])) {
          temp_MEsmeta <- metaCluData[[paste0("MEs_metabolites_clusters_",i)]]
          temp_MEsmic <- micCluData[[paste0("MEs_microbes_clusters_",i)]]
          na_flag <- apply(is.na(temp_MEsmeta), 2, sum)
          temp_MEsmeta <- temp_MEsmeta[,which(na_flag == 0)]
          na_flag <- apply(is.na(temp_MEsmic), 2, sum)
          temp_MEsmic <- temp_MEsmic[,which(na_flag == 0)]

          mic_meta_cor = matrix (NA, nrow = ncol (temp_MEsmic), ncol = ncol (temp_MEsmeta))
          dim(mic_meta_cor)
          rownames (mic_meta_cor) = colnames (temp_MEsmic)
          colnames (mic_meta_cor) = colnames (temp_MEsmeta)
          mic_meta_p = mic_meta_cor
          mic_meta_p.adjust = mic_meta_cor
          for (m in colnames (temp_MEsmeta)) {
            mic_meta_cor [ , m] = apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
              cor.test (x, temp_MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
            mic_meta_p[,m] = apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
              cor.test (x, temp_MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
          }
          if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
            mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
            {p.adjust(x,method = "BH")})
            mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
            # mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
            mic_meta_p = as.data.frame(mic_meta_p)
            mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmic)
            rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmeta)
          }else{
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
          }
          
          write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_",i,".csv"))
          cor_res[[paste0("r_",i)]] <- as.data.frame(mic_meta_cor)
          cor_res[[paste0("p_",i)]] <- as.data.frame(mic_meta_p)
          cor_res[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
        }
        return(cor_res)
      }
      
      
      ### Method 3 - Partial correlation
      if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
        corMethod = strsplit(corMethod," ")[[1]][2]
        ### All groups
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          # exist confounder
          if(is.na(confounderData)){
            stop("Lack of confounders!")
          }else{
            mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
              pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$estimate)
            mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
              pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$p.value)
          }
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }
        
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_all.csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        ### Each groups
        for (i in unique(phenoData[,1])) {
          temp_MEsmeta <- metaCluData[[paste0("MEs_metabolites_clusters_",i)]]
          temp_MEsmic <- micCluData[[paste0("MEs_microbes_clusters_",i)]]
          temp_confounderData <- confounderData[rownames(phenoData[which(phenoData[,1] == i),]),]
          # temp_confounderData <- confounderData
          mic_meta_cor = matrix (NA, nrow = ncol (temp_MEsmic), ncol = ncol (temp_MEsmeta))
          dim(mic_meta_cor)
          rownames (mic_meta_cor) = colnames (temp_MEsmic)
          colnames (mic_meta_cor) = colnames (temp_MEsmeta)
          mic_meta_p = mic_meta_cor
          mic_meta_p.adjust = mic_meta_cor
          for (m in colnames (temp_MEsmeta)) {
            # exist confounder
            if(is.na(temp_confounderData)){
              stop("Lack of confounders!")
            }else{
              mic_meta_cor [ , m] = apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                pcor.test (x, temp_MEsmeta[,m], temp_confounderData, method = corMethod)$estimate)
              mic_meta_p[,m] = apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                pcor.test (x, temp_MEsmeta[,m], temp_confounderData, method = corMethod)$p.value)
            }
          }
          if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
            mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
            {p.adjust(x,method = "BH")})
            mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
            # mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
            mic_meta_p = as.data.frame(mic_meta_p)
            mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmic)
            rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmeta)
          }else{
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
          }
          
          write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_",i,".csv"))
          cor_res[[paste0("r_",i)]] <- as.data.frame(mic_meta_cor)
          cor_res[[paste0("p_",i)]] <- as.data.frame(mic_meta_p)
          cor_res[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
        }
        return(cor_res)
      }
      
      ### Method 4 - Generalized linear models
      if (corMethod == "glm") {
        mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (MEsmic)
        colnames (mic_meta_cor) = colnames (MEsmeta)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (MEsmeta)) {
          # exist confounder
          if(is.na(confounderData)){
            mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,4])
          }else{
            temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
            mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
          }
        }
        if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
          mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
          {p.adjust(x,method = "BH")})
          mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
          # mic_meta_p.adjust = t(mic_meta_p.adjust)
        }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
          mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
          mic_meta_p = as.data.frame(mic_meta_p)
          mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
          rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
          rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
        }else{
          mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          {p.adjust(x,method = "BH")})
        }
        
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_all.csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_all.csv"))
        cor_res <- list()
        cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
        cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
        cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
        ### Each groups
        for (i in unique(phenoData[,1])) {
          temp_MEsmeta <- metaCluData[[paste0("MEs_metabolites_clusters_",i)]]
          temp_MEsmic <- micCluData[[paste0("MEs_microbes_clusters_",i)]]
          # temp_confounderData <- confounderData[rownames(phenoData[which(phenoData[,1] == i),]),]
          if(is.na(confounderData)){
            temp_confounderData <- NA
          }else{
            # temp_confounderData <- confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]
            # rownames(temp_confounderData) <- rownames(temp_MEsmeta)
            # colnames(temp_confounderData) <- colnames(confounderData)
            temp_confounderData <- confounderData
          }
          
          mic_meta_cor = matrix (NA, nrow = ncol (temp_MEsmic), ncol = ncol (temp_MEsmeta))
          dim(mic_meta_cor)
          rownames (mic_meta_cor) = colnames (temp_MEsmic)
          colnames (mic_meta_cor) = colnames (temp_MEsmeta)
          mic_meta_p = mic_meta_cor
          mic_meta_p.adjust = mic_meta_cor
          for (m in colnames (temp_MEsmeta)) {
            # exist confounder
            if(is.na(confounderData)){
              mic_meta_cor [ , m] <- apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                summary(glm(temp_MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,1])
              mic_meta_p [ , m] <- apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                summary(glm(temp_MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,4])
            }else{
              temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
              mic_meta_cor [ , m] <- apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                summary(glm(as.formula(paste0("temp_MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
              mic_meta_p [ , m] <- apply (temp_MEsmic, MARGIN = 2, FUN = function (x) 
                summary(glm(as.formula(paste0("temp_MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
            }
          }
          if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
            mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
            {p.adjust(x,method = "BH")})
            mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
            # mic_meta_p.adjust = t(mic_meta_p.adjust)
          }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
            mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
            mic_meta_p = as.data.frame(mic_meta_p)
            mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmic)
            rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(temp_MEsmeta)
          }else{
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
            {p.adjust(x,method = "BH")})
          }
          
          write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p_by_",phenotype,"_",i,".csv"))
          write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust_by_",phenotype,"_",i,".csv"))
          cor_res[[paste0("r_",i)]] <- as.data.frame(mic_meta_cor)
          cor_res[[paste0("p_",i)]] <- as.data.frame(mic_meta_p)
          cor_res[[paste0("p.adjust_",i)]] <- as.data.frame(mic_meta_p.adjust)
        }
        return(cor_res)
      }
    }
    
  }  
  
    
  ############ No phenotypes ############
  else if(is.na(phenoDataType)){
    MEsmeta <- metaCluData[["MEs_metabolites_clusters"]]
    MEsmic <- micCluData[["MEs_microbes_clusters"]]
    #####################################################################################
    ############ Step 1 - Correlate metabolite clusters to microbes clusters ############
    #####################################################################################
    
    ### The set of metabolite clusters associated with the phenotype of interest
    ### should next be tested for association with the set of microbes
    ### clusters likewise so associated (identified in Step 1). 
    ###
    ### In this step, the correlations between each metabolites cluster and
    ### microbes cluster are determined.
    
    
    ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
    if(corMethod == "gramm"){
      naivegramm <- naivegramm(MEsmeta,MEsmic,covdata = confounderData,pheno = "All")
      write.csv(naivegramm[["r_All"]],"./results/Inter-Cor/hierarchical/mic_meta_cor.csv")
      write.csv(naivegramm[["p_All"]], "./results/Inter-Cor/hierarchical/mic_meta_p.csv")
      if(nrow(naivegramm[["p_All"]]) == 1 & ncol(naivegramm[["p_All"]]) != 1){
        mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,],method = "BH")
        mic_meta_p.adjust <- as.data.frame(t(mic_meta_p.adjust))
        rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
        colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
        naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
      }else if(nrow(naivegramm[["p_All"]]) != 1 & ncol(naivegramm[["p_All"]]) == 1){
        mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][,1],method = "BH")
        mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
        rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
        colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
        naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
      }else if(nrow(naivegramm[["p"]]) == 1 & ncol(naivegramm[["p"]]) == 1){
        mic_meta_p.adjust = p.adjust(naivegramm[["p_All"]][1,1],method = "BH")
        mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
        rownames(mic_meta_p.adjust) <- rownames(naivegramm[["p_All"]])
        colnames(mic_meta_p.adjust) <- colnames(naivegramm[["p_All"]])
        naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
      }else{
        mic_meta_p.adjust = apply(naivegramm[["p_All"]], 2, FUN = function(x){
          p.adjust(x,method = "BH")})
        naivegramm[["p.adjust_All"]] <- as.data.frame(mic_meta_p.adjust)
      }
      write.csv(naivegramm[["r_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor.csv"))
      write.csv(naivegramm[["p_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p.csv"))
      write.csv(naivegramm[["p.adjust_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust.csv"))
      write.csv(naivegramm[["type_All"]], paste0("./results/Inter-Cor/hierarchical/mic_meta_cor_type_by.csv"))
      gramm_res <- naivegramm
      return(gramm_res)
    }
    
    ### Method 2 - Method 2 - Spearman,Pearson,Kendall
    if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
      mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (MEsmic)
      colnames (mic_meta_cor) = colnames (MEsmeta)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (MEsmeta)) {
        mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
          cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
        mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
          cor.test (x, MEsmeta[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
      }
      if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
        mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
        {p.adjust(x,method = "BH")})
        mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
        mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
        mic_meta_p = as.data.frame(mic_meta_p)
        mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
        rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
        rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
      }else{
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
      }
      write.csv(mic_meta_cor,"./results/Inter-Cor/hierarchical/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/hierarchical/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/hierarchical/mic_meta_p.adjust.csv")
      cor_res <- list()
      cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
      cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
      cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
      return(cor_res)
    }
    
    ### Method 3 - Partial correlation
    if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
      corMethod = strsplit(corMethod," ")[[1]][2]
      mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (MEsmic)
      colnames (mic_meta_cor) = colnames (MEsmeta)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (MEsmeta)) {
        # exist confounder
        if(is.na(confounderData)){
          stop("Lack of confounders!")
        }else{
          mic_meta_cor [ , m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$estimate)
          mic_meta_p[,m] = apply (MEsmic, MARGIN = 2, FUN = function (x) 
            pcor.test (x, MEsmeta[,m], confounderData, method = corMethod)$p.value)
        }
      }
      if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
        mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
        {p.adjust(x,method = "BH")})
        mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
        mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
        mic_meta_p = as.data.frame(mic_meta_p)
        mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
        rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
        rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
      }else{
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
      }
      write.csv(mic_meta_cor,"./results/Inter-Cor/hierarchical/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/hierarchical/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/hierarchical/mic_meta_p.adjust.csv")
      cor_res <- list()
      cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
      cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
      cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
      return(cor_res)
    }
    
    ### Method 4 - Generalized linear models
    if (corMethod == "glm") {
      mic_meta_cor = matrix (NA, nrow = ncol (MEsmic), ncol = ncol (MEsmeta))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (MEsmic)
      colnames (mic_meta_cor) = colnames (MEsmeta)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (MEsmeta)) {
        # exist confounder
        if(is.na(confounderData)){
          mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
            summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
            summary(glm(MEsmeta[,m] ~ x,family = gaussian))$coefficients[2,4])
        }else{
          temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
          mic_meta_cor [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (MEsmic, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("MEsmeta[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
        }
      }
      if(nrow(mic_meta_p) == 1 & ncol(mic_meta_p) != 1){
        mic_meta_p.adjust = apply(mic_meta_p, 1, FUN = function(x)
        {p.adjust(x,method = "BH")})
        mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) != 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
        # mic_meta_p.adjust = t(mic_meta_p.adjust)
      }else if (nrow(mic_meta_p) == 1 & ncol(mic_meta_p) == 1){
        mic_meta_p.adjust = p.adjust(mic_meta_p,method = "BH")
        mic_meta_p = as.data.frame(mic_meta_p)
        mic_meta_p.adjust = as.data.frame(mic_meta_p.adjust)
        rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(MEsmic)
        rownames(mic_meta_p.adjust) = rownames(mic_meta_p.adjust) = colnames(MEsmeta)
      }else{
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        {p.adjust(x,method = "BH")})
      }
      
      write.csv(mic_meta_cor,paste0("./results/Inter-Cor/hierarchical/mic_meta_cor.csv"))
      write.csv(mic_meta_p, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.csv"))
      write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/hierarchical/mic_meta_p.adjust.csv"))
      cor_res <- list()
      cor_res[[paste0("r_All")]] <- as.data.frame(mic_meta_cor)
      cor_res[[paste0("p_All")]] <- as.data.frame(mic_meta_p)
      cor_res[[paste0("p.adjust_All")]] <- as.data.frame(mic_meta_p.adjust)
      return(cor_res)
    }
  }
}
