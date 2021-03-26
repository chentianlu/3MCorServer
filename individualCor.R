# function - individualCor #
########################################################################
# File: individualCor.R
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
##    metaIndData --- A dataframe of metabolites                                      ##
##                       (rows: samples,columns: metabolites)                         ##
##    micIndData --- A dataframe of microbes                                          ##
##                       (rows: samples,columns: microbes)                            ##
##    phenoData --- A dataframe of phenotype data                                     ##
##                  (rows: samples,columns: phenotype)                                ##
##    phenoDataType --- Phenotype data type                                           ##
##                Default: NA                                                         ##
##    confounderData --- A dataframe of confounder                                    ##
##                       Default: NA(rows:samples,columns:confounder)                 ##
##    corMethod --- Method of the correlation between metabolome and microbiome data  ##
##                  Default: spearman                                                 ##
##    venn.rthreshold --- The correlation coefficient threshold of a Venn diagram     ##
##    venn.pthreshold --- The p value threshold of a Venn diagram                     ##
##                                                                                    ##
########################################################################################
## Output:                                                                            ## 
##    cor_res ---  A list included coefficient, p value and p.adjust value            ##
##                 of the correlation between metabolome data and microbiome          ##
##                 data                                                               ##
##                                                                                    ##
########################################################################################
#---------------------------------------------------------------------------------------------------------------------
individualCor <- function(metaIndData,micIndData,phenoData,phenoDataType = NA,confounderData = NA,corMethod = "spearman",venn.rthreshold = 0.5,venn.pthreshold = 0.05)
{
  if(!file.exists("./results/Inter-Cor/pairwise")){
    dir.create("./results/Inter-Cor/pairwise",recursive = T)
  }
  
  ############ Exist continuous phenotypes variable############
  if(!is.na(phenoDataType) & phenoDataType == "continuous"){
    ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
    if(corMethod == "gramm"){
      ### metabolites to microbes
      naivegramm_res1 <- naivegramm(metaIndData,micIndData,covdata = confounderData)
      write.csv(naivegramm_res1[["r"]],"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(naivegramm_res1[["p"]], "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      mic_meta_p.adjust = apply(naivegramm_res1[["p"]], 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      naivegramm_res1[["p.adjust" ]] = as.data.frame(mic_meta_p.adjust) 
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      write.csv(naivegramm_res1[["type"]], "./results/Inter-Cor/pairwise/mic_meta_cor_type.csv")
      
      
      ### metabolites to phenotypes
      naivegramm_res2 <- naivegramm(metaIndData,phenoData,covdata = confounderData)
      write.csv(naivegramm_res2[["r"]],"./results/Inter-Cor/pairwise/pheno_meta_cor.csv")
      write.csv(naivegramm_res2[["p"]], "./results/Inter-Cor/pairwise/pheno_meta_p.csv")
      pheno_meta_p.adjust = apply(naivegramm_res2[["p"]], 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      naivegramm_res2[["p.adjust" ]] = as.data.frame(pheno_meta_p.adjust)
      write.csv(pheno_meta_p.adjust, "./results/Inter-Cor/pairwise/pheno_meta_p.adjust.csv")
      write.csv(naivegramm_res2[["type"]], "./results/Inter-Cor/pairwise/pheno_meta_cor_type.csv")
      
      
      ### microbes to phenotypes
      naivegramm_res3 <- naivegramm(micIndData,phenoData,covdata = confounderData)
      write.csv(naivegramm_res3[["r"]],"./results/Inter-Cor/pairwise/pheno_mic_cor.csv")
      write.csv(naivegramm_res3[["p"]], "./results/Inter-Cor/pairwise/pheno_mic_p.csv")
      pheno_mic_p.adjust = apply(naivegramm_res3[["p"]], 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      naivegramm_res3[["p.adjust" ]] = as.data.frame(pheno_mic_p.adjust)
      write.csv(pheno_mic_p.adjust, "./results/Inter-Cor/pairwise/pheno_mic_p.adjust.csv")
      write.csv(naivegramm_res3[["type"]],"./results/Inter-Cor/pairwise/pheno_mic_cor_type.csv")
      
      
      res_list <- list(mic2metaCor = naivegramm_res1, pheno2metaCor = naivegramm_res2, pheno2micCor = naivegramm_res3)
      return(res_list)
      
    }
    
    ### Method 2 - Spearman,Pearson,Kendall
    if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall" ){
      ### metabolites to microbes
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        mic_meta_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
        mic_meta_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      mic_meta_result <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      
      ### metabolites to phenotypes
      pheno_meta_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (metaIndData))
      dim(pheno_meta_cor)
      rownames (pheno_meta_cor) = colnames (phenoData)
      colnames (pheno_meta_cor) = colnames (metaIndData)
      pheno_meta_p = pheno_meta_cor
      pheno_meta_p.adjust = pheno_meta_cor
      for (m in colnames (metaIndData)) {
        pheno_meta_cor [ , m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
        pheno_meta_p[,m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
      }
      pheno_meta_p.adjust = apply(pheno_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_meta_cor,"./results/Inter-Cor/pairwise/pheno_meta_cor.csv")
      write.csv(pheno_meta_p, "./results/Inter-Cor/pairwise/pheno_meta_p.csv")
      write.csv(pheno_meta_p.adjust, "./results/Inter-Cor/pairwise/pheno_meta_p.adjust.csv")
      pheno_meta_result <- list(r = as.data.frame(pheno_meta_cor), p = as.data.frame(pheno_meta_p), p.adjust = as.data.frame(pheno_meta_p.adjust))
      
      ### microbes to phenotypes
      pheno_mic_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (micIndData))
      dim(pheno_mic_cor)
      rownames (pheno_mic_cor) = colnames (phenoData)
      colnames (pheno_mic_cor) = colnames (micIndData)
      pheno_mic_p = pheno_mic_cor
      pheno_mic_p.adjust = pheno_mic_cor
      for (m in colnames (micIndData)) {
        pheno_mic_cor [ , m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
          cor.test (x, micIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
        pheno_mic_p[,m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
          cor.test (x, micIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
      }
      pheno_mic_p.adjust = apply(pheno_mic_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_mic_cor,"./results/Inter-Cor/pairwise/pheno_mic_cor.csv")
      write.csv(pheno_mic_p, "./results/Inter-Cor/pairwise/pheno_mic_p.csv")
      write.csv(pheno_mic_p.adjust, "./results/Inter-Cor/pairwise/pheno_mic_p.adjust.csv")
      pheno_mic_result <- list(r = as.data.frame(pheno_mic_cor), p = as.data.frame(pheno_mic_p), p.adjust = as.data.frame(pheno_mic_p.adjust))
      
      res_list <- list(mic2metaCor = mic_meta_result, pheno2metaCor = pheno_meta_result, pheno2micCor = pheno_mic_result)
      return(res_list)
    }
    
    ### Method 3 - Partial correlation
    if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
      corMethod = strsplit(corMethod," ")[[1]][2]
      ### metabolites to microbes
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          stop("Lack of confounders!")
        }else{
          mic_meta_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, metaIndData[,m], confounderData, method = corMethod)$estimate)
          mic_meta_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, metaIndData[,m], confounderData, method = corMethod)$p.value)
        }
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      mic_meta_result <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      
      ### metabolites to phenotypes
      pheno_meta_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (metaIndData))
      dim(pheno_meta_cor)
      rownames (pheno_meta_cor) = colnames (phenoData)
      colnames (pheno_meta_cor) = colnames (metaIndData)
      pheno_meta_p = pheno_meta_cor
      pheno_meta_p.adjust = pheno_meta_cor
      for (m in colnames (metaIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          stop("Lack of confounders!")
        }else{
          pheno_meta_cor [ , m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, metaIndData[,m], confounderData, method = corMethod)$estimate)
          pheno_meta_p[,m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, metaIndData[,m], confounderData, method = corMethod)$p.value)
        }
      }
      pheno_meta_p.adjust = apply(pheno_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_meta_cor,"./results/Inter-Cor/pairwise/pheno_meta_cor.csv")
      write.csv(pheno_meta_p, "./results/Inter-Cor/pairwise/pheno_meta_p.csv")
      write.csv(pheno_meta_p.adjust, "./results/Inter-Cor/pairwise/pheno_meta_p.adjust.csv")
      pheno_meta_result <- list(r = as.data.frame(pheno_meta_cor), p = as.data.frame(pheno_meta_p), p.adjust = as.data.frame(pheno_meta_p.adjust))
      
      ### microbes to phenotypes
      pheno_mic_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (micIndData))
      dim(pheno_mic_cor)
      rownames (pheno_mic_cor) = colnames (phenoData)
      colnames (pheno_mic_cor) = colnames (micIndData)
      pheno_mic_p = pheno_mic_cor
      pheno_mic_p.adjust = pheno_mic_cor
      for (m in colnames (micIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          stop("Lack of confounders!")
        }else{
          pheno_mic_cor [ , m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, micIndData[,m], confounderData, method = corMethod)$estimate)
          pheno_mic_p[,m] = apply (phenoData, MARGIN = 2, FUN = function (x) 
            ppcor::pcor.test (x, micIndData[,m], confounderData, method = corMethod)$p.value)
        }
      }
      pheno_mic_p.adjust = apply(pheno_mic_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_mic_cor,"./results/Inter-Cor/pairwise/pheno_mic_cor.csv")
      write.csv(pheno_mic_p, "./results/Inter-Cor/pairwise/pheno_mic_p.csv")
      write.csv(pheno_mic_p.adjust, "./results/Inter-Cor/pairwise/pheno_mic_p.adjust.csv")
      pheno_mic_result <- list(r = as.data.frame(pheno_mic_cor), p = as.data.frame(pheno_mic_p), p.adjust = as.data.frame(pheno_mic_p.adjust))
      
      res_list <- list(mic2metaCor = mic_meta_result, pheno2metaCor = pheno_meta_result, pheno2micCor = pheno_mic_result)
      return(res_list)
    }
    
    ### Method 4 - Generalized linear models
    if (corMethod == "glm") {
      ### metabolites to microbes
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          mic_meta_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
        }else{
          temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
          
          mic_meta_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
        }
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      mic_meta_result <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      
      ### metabolites to phenotypes
      pheno_meta_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (metaIndData))
      dim(pheno_meta_cor)
      rownames (pheno_meta_cor) = colnames (phenoData)
      colnames (pheno_meta_cor) = colnames (metaIndData)
      pheno_meta_p = pheno_meta_cor
      pheno_meta_p.adjust = pheno_meta_cor
      for (m in colnames (metaIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          pheno_meta_cor [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
          pheno_meta_p [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
        }else{
          temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
          pheno_meta_cor [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
          pheno_meta_p [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
        }
      }
      pheno_meta_p.adjust = apply(pheno_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_meta_cor,"./results/Inter-Cor/pairwise/pheno_meta_cor.csv")
      write.csv(pheno_meta_p, "./results/Inter-Cor/pairwise/pheno_meta_p.csv")
      write.csv(pheno_meta_p.adjust, "./results/Inter-Cor/pairwise/pheno_meta_p.adjust.csv")
      pheno_meta_result <- list(r = as.data.frame(pheno_meta_cor), p = as.data.frame(pheno_meta_p), p.adjust = as.data.frame(pheno_meta_p.adjust))
      
      ### microbes to phenotypes
      pheno_mic_cor = matrix (NA, nrow = ncol (phenoData), ncol = ncol (micIndData))
      dim(pheno_mic_cor)
      rownames (pheno_mic_cor) = colnames (phenoData)
      colnames (pheno_mic_cor) = colnames (micIndData)
      pheno_mic_p = pheno_mic_cor
      pheno_mic_p.adjust = pheno_mic_cor
      for (m in colnames (micIndData)) {
        # exist confounder
        if(is.na(confounderData)){
          pheno_mic_cor [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(micIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
          pheno_mic_p [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(micIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
        }else{
          temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
          pheno_mic_cor [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("micIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
          pheno_mic_p [ , m] <- apply (phenoData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("micIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
        }
      }
      pheno_mic_p.adjust = apply(pheno_mic_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(pheno_mic_cor,"./results/Inter-Cor/pairwise/pheno_mic_cor.csv")
      write.csv(pheno_mic_p, "./results/Inter-Cor/pairwise/pheno_mic_p.csv")
      write.csv(pheno_mic_p.adjust, "./results/Inter-Cor/pairwise/pheno_mic_p.adjust.csv")
      pheno_mic_result <- list(r = as.data.frame(pheno_mic_cor), p = as.data.frame(pheno_mic_p), p.adjust = as.data.frame(pheno_mic_p.adjust))
      
      res_list <- list(mic2metaCor = mic_meta_result, pheno2metaCor = pheno_meta_result, pheno2micCor = pheno_mic_result)
      return(res_list)
    }
  }
  
  ############ Exist categorical phenotypes variable ############ 
  else if(!is.na(phenoDataType) & phenoDataType == "categorical"){
    ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
    if(corMethod == "gramm"){
      cate_cor_func <- function(metaInddf,micInddf,confdf,phenotype){
        ### Cor between all metabolites and microbes
        naivegramm <- naivegramm(metaInddf,micInddf,covdata = confdf,pheno = phenotype)
        write.csv(naivegramm[[paste0("r_",phenotype)]],paste0("./results/Inter-Cor/pairwise/mic_meta_cor_",phenotype,".csv"))
        write.csv(naivegramm[[paste0("p_",phenotype)]], paste0("./results/Inter-Cor/pairwise/mic_meta_p_",phenotype,".csv"))
        mic_meta_p.adjust = apply(naivegramm[[paste0("p_",phenotype)]], 2, FUN = function(x)
          p.adjust(x,method = "BH"))
        naivegramm[[paste0("p.adjust_",phenotype)]] = as.data.frame(mic_meta_p.adjust)
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/pairwise/mic_meta_p.adjust_",phenotype,".csv"))
        write.csv(naivegramm[[paste0("type_",phenotype)]], paste0("./results/Inter-Cor/pairwise/mic_meta_cor_type_",phenotype,".csv"))
        res_list <- naivegramm
        
        if(phenotype == "All"){
          ### Bar plot,box plot in all groups
          mic_meta_cor = naivegramm[[paste0("r_",phenotype)]]
          dim(mic_meta_cor)
          rownames (mic_meta_cor) = colnames (micIndData)
          colnames (mic_meta_cor) = colnames (metaIndData)
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 was counted
          counts <- 0 
          for (i in 1:nrow(mic_meta_cor)) {
            for (j in 1:ncol(mic_meta_cor)) {
              if(abs(mic_meta_cor[i,j]) > 0.5 && mic_meta_p.adjust[i,j] < 0.05){
                counts <- counts + 1
                # print(paste0("i:",i,"j:",j))
              }
            }
          }
          # counts <- sum(abs(naivegramm[[paste0("r_",phenotype)]] > 0.5))
          r_df <- data.frame("All groups",counts)
          counts <- sum(mic_meta_p.adjust < 0.05)
          p.adjust_df <- data.frame("All groups",counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df <- rbind(r_df,p.adjust_df)
          colnames(bar_df)[1] <- "group_name"
          
          ### Box plot
          box_df <- naivegramm[[paste0("r_",phenotype)]]
          box_df <- as.matrix(box_df)
          dim(box_df) <- c(ncol(box_df) * nrow(box_df),1)
          box_df <- as.data.frame(box_df)
          box_df$group <- "All groups"
          colnames(box_df) <- c("r_value","Groups")
          
          ### Bar plot,box plot and Venn Diagram between each group
          
          ### Identify groups
          temp_phenoData <- phenoData
          temp_phenoData$SampleID <- rownames(temp_phenoData)
          
          temp_metaInddf <- metaInddf
          temp_metaInddf$SampleID <- rownames(temp_metaInddf)
          
          temp_micInddf <- micInddf
          temp_micInddf$SampleID <- rownames(temp_micInddf)

          if(is.na(confdf)){
            temp_covdata <- NA
          }else{
            temp_covdata <- confdf
            temp_covdata$SampleID <- rownames(temp_covdata)
            temp_covdata <- merge(temp_covdata,temp_phenoData,all.x = T)
          }
          
          temp_metaInddf <- merge(temp_metaInddf,temp_phenoData,all.x = T)
          temp_micInddf <- merge(temp_micInddf,temp_phenoData,all.x = T)
          

          ### Cor between metabolites and microbes in each group
          unique(phenoData[,1]) 
          group_res_list <- list()
          for(i in 1:length(unique(phenoData[,1]))){
            group_metaInddf <- temp_metaInddf[which(temp_metaInddf[,ncol(temp_metaInddf)]== unique(phenoData[,1])[i]),]
            group_metaInddf <- group_metaInddf[,-1]
            group_metaInddf <- group_metaInddf[,-ncol(group_metaInddf)]
            
            group_micInddf <- temp_micInddf[which(temp_micInddf[,ncol(temp_micInddf)]== unique(phenoData[,1])[i]),]
            group_micInddf <- group_micInddf[,-1]
            group_micInddf <- group_micInddf[,-ncol(group_micInddf)]
            
            if(is.na(confdf)){
              group_confdf <- NA
            }else{
              group_confdf <- temp_covdata[which(temp_covdata[,ncol(temp_covdata)]== unique(phenoData[,1])[i]),]
              group_confdf <- group_confdf[,-1]
              group_confdf <- group_confdf[,-ncol(group_confdf)]
              group_confdf <- as.data.frame(group_confdf)
              rownames(group_confdf) <- rownames(group_micInddf)
            }

            group_naivegramm <- naivegramm(group_metaInddf,group_micInddf,covdata = group_confdf,pheno = phenotype)
            p.adjust = apply(group_naivegramm[[paste0("p_",phenotype)]], 2, FUN = function(x)
              p.adjust(x,method = "BH"))
            group_naivegramm[[paste0("p.adjust_",phenotype)]] = as.data.frame(p.adjust)
            
            group_res_list[[unique(phenoData[,1])[i]]] <- group_naivegramm
          }
          
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 was counted
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            temp_counts <- 0
            for(j in 1: nrow(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])){
              for (k in 1:ncol(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])) {
                if(abs(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]][j,k]) > 0.5 && group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]][j,k] < 0.05){
                  temp_counts <- temp_counts + 1
                }
              }
            }
            counts <- append(counts,temp_counts) 
          }
          group_name <- unique(phenoData[,1])
          r_df <- data.frame(group_name,counts)
          
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            counts <- append(counts,sum(group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]] < 0.05)) 
          }
          p.adjust_df <- data.frame(group_name,counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df2 <- rbind(r_df,p.adjust_df)
          bar_df <- rbind(bar_df,bar_df2)
          ### Bar plot
          bar_plot <- ggplot(bar_df,aes(group_name,counts))+
            geom_bar(stat="identity",position="dodge",fill = "#00b692")+
            labs(x = 'Groups', y = 'Counts of correlation pairs', size = 2) +
            theme(legend.position='none')+
            facet_grid(object ~ .)+
            coord_flip()
          
          ggsave(paste0("./results/Inter-Cor/pairwise/Counts_between_groups_bar_plot.pdf"),width = 7,height = 5)
          
          ### Box plot
          box_df2 <- group_res_list[[unique(phenoData[,1])[1]]][[paste0("r_",phenotype)]]
          box_df2 <- as.matrix(box_df2)
          dim(box_df2) <- c(ncol(box_df2) * nrow(box_df2),1)
          box_df2 <- as.data.frame(box_df2)
          box_df2$group <- unique(phenoData[,1])[1]
          for (i in 2:length(unique(phenoData[,1]))){
            temp_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_df <- as.matrix(temp_df)
            dim(temp_df) <- c(ncol(temp_df) * nrow(temp_df),1)
            temp_df <- as.data.frame(temp_df)
            temp_df$group <- unique(phenoData[,1])[i]
            box_df2 <- rbind(box_df2,temp_df)
          }
          colnames(box_df2) <- c("r_value","Groups")
          
          box_df <- rbind(box_df,box_df2)
          compaired <- combn(unique(box_df$Groups),2,simplify = F)
          
          box_plot <- ggplot(box_df,aes(Groups,r_value,fill = Groups)) + 
            geom_boxplot(width = 0.6,colour = "black") + 
            theme(plot.title=element_text(size = 15),
                  axis.text.x=element_text(size=10,angle=0),
                  axis.text.y=element_text(size=10),
                  axis.title.x=element_text(size = 15),
                  axis.title.y=element_text(size = 15)) + 
            labs(x='Groups', y= 'R value') + 
            geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
          ggsave(paste0("./results/Inter-Cor/pairwise/R_value_between_groups_box_plot.pdf"),width = 7,height = 5)
          
          ### Venn Diagram
          ### Substances with |r| > venn.rthreshold && p < venn.pthreshold were screened, record the subscript 
          meta_list <- list()
          mic_list <- list()
          for (i in 1:length(unique(phenoData[,1]))){
            temp_r_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_p.adjust_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]]
            meta_df <- c()
            mic_df <- c()
            for (j in 1:nrow(temp_r_df)){
              for (k in 1:ncol(temp_r_df)){
                if(abs(temp_r_df[j,k]) >= venn.rthreshold && temp_p.adjust_df[j,k] < venn.pthreshold){
                  meta_df <- append(meta_df,k)
                  mic_df <- append(mic_df,j)
                }
              }
            }
            meta_list[[unique(phenoData[,1])[i]]] = meta_df
            mic_list[[unique(phenoData[,1])[i]]] = mic_df
          }
          if(length(meta_list) != length(unique(phenoData[,1])) || length(mic_list) != length(unique(phenoData[,1]))){
            print(paste0("There was no pair of metabolites and microbes with |r| >",venn.rthreshold," and p-value < ",venn.pthreshold," between different groups!"))
          }else{
            # Meta and MIC subscripts are combined to represent pairs names, 
            # and if the same pair names are present, a correlation pair with an intersection is indicated
            
            pairs_df_vn <- list()
            for (i in 1:length(unique(phenoData[,1]))) {
              for (j in 1:length(mic_list[[i]])) {
                pairs_df_vn[[unique(phenoData[,1])[i]]][j] <- paste0(colnames(group_res_list[[i]][[paste0("r_",phenotype)]])[meta_list[[i]][j]],"-",rownames(group_res_list[[i]][[paste0("r_",phenotype)]])[mic_list[[i]][j]])
              }
              # Output significant correlation pairs for each group
              write.csv(pairs_df_vn[[unique(phenoData[,1])[i]]],paste0("./results/Inter-Cor/pairwise/",unique(phenoData[,1])[i],"_sig_cor_pairs.csv"),row.names = F)
            }
            
            # Output intersections of significant correlation pairs for each group
            inter_pairs <- Reduce(intersect,pairs_df_vn)
            write.csv(inter_pairs,"./results/Inter-Cor/pairwise/inter_pairs.csv")
            
            fill_col <- rainbow(length(unique(phenoData[,1])))
            
            ## pairs_null_flag:Determines whether no element exists
            pairs_null_flag = 0
            for (i in 1:length(pairs_df_vn)) {
              if(length(pairs_df_vn[[i]]) == 0){
                pairs_null_flag = 1
              }
            }
            if(pairs_null_flag == 0){
              venn.plot <- venn.diagram(x = pairs_df_vn,NULL,col = "black",fill = fill_col,
                                        cat.cex = 2,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")
              pdf(file="./results/Inter-Cor/pairwise/vn_plot_of_sig_cor_pairs_between_groups.pdf")
              grid.draw(venn.plot)  # library(grDevices)
              dev.off()
            }
          }
        }
        return(res_list)
      }
      
      res_list.list <- list()
      # gramm for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        if(is.na(confounderData)){
          tempconfounderData <- NA
        }else{
          tempconfounderData <- as.data.frame(confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]) 
          rownames(tempconfounderData) <- rownames(tempmetaData)
          colnames(tempconfounderData) <- colnames(confounderData)
        }
        temp_list <- cate_cor_func(metaInddf = tempmetaData,micInddf = tempmicData,confdf = tempconfounderData, phenotype = i)
        res_list.list <- append(res_list.list,temp_list)
      }
      
      # gramm for all samples
      naivegramm <- cate_cor_func(metaInddf = metaIndData,micInddf = micIndData,confdf= confounderData, phenotype = "All")

      res_list.list <- append(res_list.list,naivegramm)
      return(res_list.list)
    }
    
    ### Method 2 - Spearman,Pearson,Kendall
    if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall" ){
      cate_cor_func <- function(metaInddf,micInddf,phenotype,method){
        mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (micInddf)
        colnames (mic_meta_cor) = colnames (metaInddf)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (metaInddf)) {
          mic_meta_cor [ , m] = apply (micInddf, MARGIN = 2, FUN = function (x) 
            cor.test (x, metaInddf[,m], method = method, use = "pairwise.complete.obs")$estimate)
          mic_meta_p[,m] = apply (micInddf, MARGIN = 2, FUN = function (x) 
            cor.test (x, metaInddf[,m], method = method, use = "pairwise.complete.obs")$p.value)
        }
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          p.adjust(x,method = "BH"))
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_cor_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/pairwise/mic_meta_p_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/pairwise/mic_meta_p.adjust_",phenotype,".csv"))
        res_list <- list()
        res_list[[paste0("r_",phenotype)]] = as.data.frame(mic_meta_cor)
        res_list[[paste0("p_",phenotype)]] = as.data.frame(mic_meta_p)
        res_list[[paste0("p.adjust_",phenotype)]] = as.data.frame(mic_meta_p.adjust)

        if(phenotype == "All"){
          ### Bar plot,box plot in all groups
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 && p<0.05 was counted
          counts <- 0 
          for (i in 1:nrow(mic_meta_cor)) {
            for (j in 1:ncol(mic_meta_cor)) {
              if(abs(mic_meta_cor[i,j]) > 0.5 && mic_meta_p.adjust[i,j] < 0.05){
                counts <- counts + 1
                # print(paste0("i:",i,"j:",j))
              }
            }
          }
          # counts <- sum(abs(mic_meta_cor > 0.5))
          r_df <- data.frame("All groups",counts)
          counts <- sum(mic_meta_p.adjust < 0.05)
          p.adjust_df <- data.frame("All groups",counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df <- rbind(r_df,p.adjust_df)
          colnames(bar_df)[1] <- "group_name"
          
          ### Box plot
          box_df <- mic_meta_cor
          box_df <- as.matrix(box_df)
          dim(box_df) <- c(ncol(box_df) * nrow(box_df),1)
          box_df <- as.data.frame(box_df)
          box_df$group <- "All groups"
          colnames(box_df) <- c("r_value","Groups")
          
          ### Bar plot,box plot and Venn Diagram between each group
          
          ### Identify groups
          temp_phenoData <- phenoData
          temp_phenoData$SampleID <- rownames(temp_phenoData)
          
          temp_metaInddf <- metaInddf
          temp_metaInddf$SampleID <- rownames(temp_metaInddf)
          
          temp_micInddf <- micInddf
          temp_micInddf$SampleID <- rownames(temp_micInddf)

          temp_metaInddf <- merge(temp_metaInddf,temp_phenoData,all.x = T)
          temp_micInddf <- merge(temp_micInddf,temp_phenoData,all.x = T)
          
          ### Cor between metabolites and microbes in each group
          unique(phenoData[,1]) 
          group_res_list <- list()
          for(i in 1:length(unique(phenoData[,1]))){
            group_metaInddf <- temp_metaInddf[which(temp_metaInddf[,ncol(temp_metaInddf)]== unique(phenoData[,1])[i]),]
            group_metaInddf <- group_metaInddf[,-1]
            group_metaInddf <- group_metaInddf[,-ncol(group_metaInddf)]
            
            group_micInddf <- temp_micInddf[which(temp_micInddf[,ncol(temp_micInddf)]== unique(phenoData[,1])[i]),]
            group_micInddf <- group_micInddf[,-1]
            group_micInddf <- group_micInddf[,-ncol(group_micInddf)]

            mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
            dim(mic_meta_cor)
            rownames (mic_meta_cor) = colnames (micInddf)
            colnames (mic_meta_cor) = colnames (metaInddf)
            mic_meta_p = mic_meta_cor
            mic_meta_p.adjust = mic_meta_cor
            for (m in colnames (group_metaInddf)) {
              mic_meta_cor [ , m] = apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                cor.test (x, group_metaInddf[,m],method = method, use = "pairwise.complete.obs")$estimate)
              mic_meta_p[,m] = apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                cor.test (x, group_metaInddf[,m],method = method, use = "pairwise.complete.obs")$p.value)
            }
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
              p.adjust(x,method = "BH"))
            
            
            
            mic_meta_cor <- as.data.frame(mic_meta_cor)
            mic_meta_p <- as.data.frame(mic_meta_p)
            mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_cor) = rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(group_micInddf)
            colnames(mic_meta_cor) = colnames(mic_meta_p) = colnames(mic_meta_p.adjust) = colnames(group_metaInddf)
            
            ### Delete NA
            mic_meta_cor <- na.omit(mic_meta_cor)
            mic_meta_p <- na.omit(mic_meta_p)
            mic_meta_p.adjust <- na.omit(mic_meta_p.adjust)
            
            # write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_cor.csv"))
            # write.csv(mic_meta_p,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.csv"))
            # write.csv(mic_meta_p.adjust,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.adjust.csv"))

            group_spearman <- list()
            group_spearman[[paste0("r_",phenotype)]]  <- mic_meta_cor
            group_spearman[[paste0("p_",phenotype)]]  <- mic_meta_p
            group_spearman[[paste0("p.adjust_",phenotype)]]  <- mic_meta_p.adjust
            
            group_res_list[[unique(phenoData[,1])[i]]] <- group_spearman
          }
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 was counted
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            temp_counts <- 0
            for(j in 1: nrow(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])){
              for (k in 1:ncol(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])) {
                if(abs(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]][j,k]) > 0.5 && group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]][j,k] < 0.05){
                  temp_counts <- temp_counts + 1
                  # print(paste0("j:",j,"k:",k))
                }
              }
            }
            counts <- append(counts,temp_counts) 
          }
          group_name <- unique(phenoData[,1])
          r_df <- data.frame(group_name,counts)
          
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            counts <- append(counts,sum(group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]] < 0.05)) 
          }
          p.adjust_df <- data.frame(group_name,counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df2 <- rbind(r_df,p.adjust_df)
          
          bar_df <- rbind(bar_df,bar_df2)
          ### Bar plot
          bar_plot <- ggplot(bar_df,aes(group_name,counts))+
            geom_bar(stat="identity",position="dodge",fill = "#00b692")+
            labs(x = 'Groups', y = 'Counts of correlation pairs', size = 2) +
            theme(legend.position='none')+
            facet_grid(object ~ .)+
            coord_flip()
          
          ggsave(paste0("./results/Inter-Cor/pairwise/Counts_between_groups_bar_plot.pdf"),width = 7,height = 5)
          
          ### Box plot
          box_df2 <- group_res_list[[unique(phenoData[,1])[1]]][[paste0("r_",phenotype)]]
          box_df2 <- as.matrix(box_df2)
          dim(box_df2) <- c(ncol(box_df2) * nrow(box_df2),1)
          box_df2 <- as.data.frame(box_df2)
          box_df2$group <- unique(phenoData[,1])[1]
          for (i in 2:length(unique(phenoData[,1]))){
            temp_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_df <- as.matrix(temp_df)
            dim(temp_df) <- c(ncol(temp_df) * nrow(temp_df),1)
            temp_df <- as.data.frame(temp_df)
            temp_df$group <- unique(phenoData[,1])[i]
            box_df2 <- rbind(box_df2,temp_df)
          }
          colnames(box_df2) <- c("r_value","Groups")
          
          box_df <- rbind(box_df,box_df2)
          compaired <- combn(unique(box_df$Groups),2,simplify = F)
          
          box_plot <- ggplot(box_df,aes(Groups,r_value,fill = Groups)) + 
            geom_boxplot(width = 0.6,colour = "black") + 
            theme(plot.title=element_text(size = 15),
                  axis.text.x=element_text(size=10,angle=0),
                  axis.text.y=element_text(size=10),
                  axis.title.x=element_text(size = 15),
                  axis.title.y=element_text(size = 15)) + 
            labs(x='Groups', y= 'R value') + 
            geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
          ggsave(paste0("./results/Inter-Cor/pairwise/R_value_between_groups_box_plot.pdf"),width = 7,height = 5)
          
          ### Venn Diagram
          ### Substances with |r| > venn.rthreshold && p < venn.pthreshold were screened, record the subscript 
          meta_list <- list()
          mic_list <- list()
          for (i in 1:length(unique(phenoData[,1]))){
            temp_r_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_p.adjust_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]]
            meta_df <- c()
            mic_df <- c()
            for (j in 1:nrow(temp_r_df)){
              for (k in 1:ncol(temp_r_df)){
                if(abs(temp_r_df[j,k]) >= venn.rthreshold && temp_p.adjust_df[j,k] < venn.pthreshold){
                  meta_df <- append(meta_df,k)
                  mic_df <- append(mic_df,j)
                }
              }
            }
            meta_list[[unique(phenoData[,1])[i]]] = meta_df
            mic_list[[unique(phenoData[,1])[i]]] = mic_df
          }
          
          if(length(meta_list) != length(unique(phenoData[,1])) || length(mic_list) != length(unique(phenoData[,1]))){
            print(paste0("There was no pair of metabolites and microbes with |r| >",venn.rthreshold," and p-value < ",venn.pthreshold," between different groups!"))
          }else{
            # Meta and MIC subscripts are combined to represent pairs names, 
            # and if the same pair names are present, a correlation pair with an intersection is indicated
            pairs_df_vn <- list()
            for (i in 1:length(unique(phenoData[,1]))) {
              for (j in 1:length(mic_list[[i]])) {
                pairs_df_vn[[unique(phenoData[,1])[i]]][j] <- paste0(colnames(group_res_list[[i]][[paste0("r_",phenotype)]])[meta_list[[i]][j]],"-",rownames(group_res_list[[i]][[paste0("r_",phenotype)]])[mic_list[[i]][j]])
              }
              # Output significant correlation pairs for each group
              write.csv(pairs_df_vn[[unique(phenoData[,1])[i]]],paste0("./results/Inter-Cor/pairwise/",unique(phenoData[,1])[i],"_sig_cor_pairs.csv"),row.names = F)
            }
            
            # Output intersections of significant correlation pairs for each group
            inter_pairs <- Reduce(intersect,pairs_df_vn)
            write.csv(inter_pairs,"./results/Inter-Cor/pairwise/inter_pairs.csv")
            
            fill_col <- rainbow(length(unique(phenoData[,1])))
            
            ## pairs_null_flag:Determines whether no element exists
            pairs_null_flag = 0
            for (i in 1:length(pairs_df_vn)) {
              if(length(pairs_df_vn[[i]]) == 0){
                pairs_null_flag = 1
              }
            }
            if(pairs_null_flag == 0){
              venn.plot <- venn.diagram(x = pairs_df_vn,NULL,col = "black",fill = fill_col,
                                        cat.cex = 2,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")
              pdf(file="./results/Inter-Cor/pairwise/vn_plot_of_sig_cor_pairs_between_groups.pdf")
              grid.draw(venn.plot)  # library(grDevices)
              dev.off()
            }
          }
        }
        return(res_list)
      }

      res_list.list <- list()
      # correlation for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        # tempconfounderData <- as.data.frame(confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]) 
        # rownames(tempconfounderData) <- rownames(tempmetaData)
        # colnames(tempconfounderData) <- colnames(confounderData)
        temp_list <- cate_cor_func(metaInddf = tempmetaData,micInddf = tempmicData,phenotype = i,method = corMethod)
        res_list.list <- append(res_list.list,temp_list)
      }
      
      # correlation for all samples
      temp_list <- cate_cor_func(metaInddf = metaIndData,micInddf = micIndData,phenotype = "All",method = corMethod)
      
      res_list.list <- append(res_list.list,temp_list)
      return(res_list.list)
    }
    
    ### Method 3 - Partial correlation
    if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
      cate_cor_func <- function(metaInddf,micInddf,confdf,phenotype,method){
        method = strsplit(method," ")[[1]][2]
        mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (micInddf)
        colnames (mic_meta_cor) = colnames (metaInddf)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (metaInddf)) {
          # exist confounder
          if(is.na(confdf)){
            stop("Lack of confounders!")
          }else{
            mic_meta_cor [ , m] = apply (micInddf, MARGIN = 2, FUN = function (x) 
              ppcor::pcor.test (x, metaInddf[,m], confdf, method = method)$estimate)
            mic_meta_p[,m] = apply (micInddf, MARGIN = 2, FUN = function (x) 
              ppcor::pcor.test (x, metaInddf[,m], confdf, method = method)$p.value)
          }
        }
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          p.adjust(x,method = "BH"))
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_cor_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/pairwise/mic_meta_p_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/pairwise/mic_meta_p.adjust_",phenotype,".csv"))
        res_list <- list()
        res_list[[paste0("r_",phenotype)]] = as.data.frame(mic_meta_cor)
        res_list[[paste0("p_",phenotype)]] = as.data.frame(mic_meta_p)
        res_list[[paste0("p.adjust_",phenotype)]] = as.data.frame(mic_meta_p.adjust)

        if(phenotype == "All"){
          ### Bar plot,box plot in all groups
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 && p<0.05 was counted
          counts <- 0 
          for (i in 1:nrow(mic_meta_cor)) {
            for (j in 1:ncol(mic_meta_cor)) {
              if(abs(mic_meta_cor[i,j]) > 0.3 && mic_meta_p.adjust[i,j] < 0.05){
                counts <- counts + 1
                # print(paste0("i:",i,"j:",j))
              }
            }
          }
          # counts <- sum(abs(mic_meta_cor > 0.5))
          r_df <- data.frame("All groups",counts)
          counts <- sum(mic_meta_p.adjust < 0.05)
          p.adjust_df <- data.frame("All groups",counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df <- rbind(r_df,p.adjust_df)
          colnames(bar_df)[1] <- "group_name"
          
          ### Box plot
          box_df <- mic_meta_cor
          box_df <- as.matrix(box_df)
          dim(box_df) <- c(ncol(box_df) * nrow(box_df),1)
          box_df <- as.data.frame(box_df)
          box_df$group <- "All groups"
          colnames(box_df) <- c("r_value","Groups")
          
          ### Bar plot,box plot and Venn Diagram between each group
          
          ### Identify groups
          temp_phenoData <- phenoData
          temp_phenoData$SampleID <- rownames(temp_phenoData)
          
          temp_metaInddf <- metaInddf
          temp_metaInddf$SampleID <- rownames(temp_metaInddf)
          
          temp_micInddf <- micInddf
          temp_micInddf$SampleID <- rownames(temp_micInddf)
          
          if(is.na(confdf)){
            temp_covdata <- NA
          }else{
            temp_covdata <- confdf
            temp_covdata$SampleID <- rownames(temp_covdata)
            temp_covdata <- merge(temp_covdata,temp_phenoData,all.x = T)
          }
          
          temp_metaInddf <- merge(temp_metaInddf,temp_phenoData,all.x = T)
          temp_micInddf <- merge(temp_micInddf,temp_phenoData,all.x = T)

          ### Cor between metabolites and microbes in each group
          unique(phenoData[,1]) 
          group_res_list <- list()
          for(i in 1:length(unique(phenoData[,1]))){
            group_metaInddf <- temp_metaInddf[which(temp_metaInddf[,ncol(temp_metaInddf)]== unique(phenoData[,1])[i]),]
            group_metaInddf <- group_metaInddf[,-1]
            group_metaInddf <- group_metaInddf[,-ncol(group_metaInddf)]
            
            group_micInddf <- temp_micInddf[which(temp_micInddf[,ncol(temp_micInddf)]== unique(phenoData[,1])[i]),]
            group_micInddf <- group_micInddf[,-1]
            group_micInddf <- group_micInddf[,-ncol(group_micInddf)]
            
            if(is.na(confdf)){
              group_confdf <- NA
            }else{
              group_confdf <- temp_covdata[which(temp_covdata[,ncol(temp_covdata)]== unique(phenoData[,1])[i]),]
              group_confdf <- group_confdf[,-1]
              group_confdf <- group_confdf[,-ncol(group_confdf)]
              group_confdf <- as.data.frame(group_confdf)
              rownames(group_confdf) <- rownames(group_micInddf)
            }
            
            mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
            dim(mic_meta_cor)
            rownames (mic_meta_cor) = colnames (micInddf)
            colnames (mic_meta_cor) = colnames (metaInddf)
            mic_meta_p = mic_meta_cor
            mic_meta_p.adjust = mic_meta_cor
            for (m in colnames (group_metaInddf)) {
              mic_meta_cor [ , m] = apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                ppcor::pcor.test (x, group_metaInddf[,m], group_confdf,method = method)$estimate)
              mic_meta_p[,m] = apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                ppcor::pcor.test (x, group_metaInddf[,m], group_confdf,method = method)$p.value)
            }
            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
              p.adjust(x,method = "BH"))
            
            mic_meta_cor <- as.data.frame(mic_meta_cor)
            mic_meta_p <- as.data.frame(mic_meta_p)
            mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_cor) = rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(group_micInddf)
            colnames(mic_meta_cor) = colnames(mic_meta_p) = colnames(mic_meta_p.adjust) = colnames(group_metaInddf)
            
            ### Delete NA
            mic_meta_cor <- na.omit(mic_meta_cor)
            mic_meta_p <- na.omit(mic_meta_p)
            mic_meta_p.adjust <- na.omit(mic_meta_p.adjust)
            
            # write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_cor.csv"))
            # write.csv(mic_meta_p,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.csv"))
            # write.csv(mic_meta_p.adjust,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.adjust.csv"))
            
            group_spearman <- list()
            group_spearman[[paste0("r_",phenotype)]]  <- mic_meta_cor
            group_spearman[[paste0("p_",phenotype)]]  <- mic_meta_p
            group_spearman[[paste0("p.adjust_",phenotype)]]  <- mic_meta_p.adjust
            group_res_list[[unique(phenoData[,1])[i]]] <- group_spearman
          }
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 was counted
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            temp_counts <- 0
            for(j in 1: nrow(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])){
              for (k in 1:ncol(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])) {
                if(abs(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]][j,k]) > 0.5 && group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]][j,k] < 0.05){
                  temp_counts <- temp_counts + 1
                  # print(paste0("j:",j,"k:",k))
                }
              }
            }
            counts <- append(counts,temp_counts) 
          }
          group_name <- unique(phenoData[,1])
          r_df <- data.frame(group_name,counts)
          
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            counts <- append(counts,sum(group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]] < 0.05)) 
          }
          p.adjust_df <- data.frame(group_name,counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df2 <- rbind(r_df,p.adjust_df)
          
          bar_df <- rbind(bar_df,bar_df2)
          ### Bar plot
          bar_plot <- ggplot(bar_df,aes(group_name,counts))+
            geom_bar(stat="identity",position="dodge",fill = "#00b692")+
            labs(x = 'Groups', y = 'Counts of correlation pairs', size = 2) +
            theme(legend.position='none')+
            facet_grid(object ~ .)+
            coord_flip()
          
          ggsave(paste0("./results/Inter-Cor/pairwise/Counts_between_groups_bar_plot.pdf"),width = 7,height = 5)
          
          ### Box plot
          box_df2 <- group_res_list[[unique(phenoData[,1])[1]]][[paste0("r_",phenotype)]]
          box_df2 <- as.matrix(box_df2)
          dim(box_df2) <- c(ncol(box_df2) * nrow(box_df2),1)
          box_df2 <- as.data.frame(box_df2)
          box_df2$group <- unique(phenoData[,1])[1]
          for (i in 2:length(unique(phenoData[,1]))){
            temp_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_df <- as.matrix(temp_df)
            dim(temp_df) <- c(ncol(temp_df) * nrow(temp_df),1)
            temp_df <- as.data.frame(temp_df)
            temp_df$group <- unique(phenoData[,1])[i]
            box_df2 <- rbind(box_df2,temp_df)
          }
          colnames(box_df2) <- c("r_value","Groups")
          
          box_df <- rbind(box_df,box_df2)
          compaired <- combn(unique(box_df$Groups),2,simplify = F)
          
          box_plot <- ggplot(box_df,aes(Groups,r_value,fill = Groups)) + 
            geom_boxplot(width = 0.6,colour = "black") + 
            theme(plot.title=element_text(size = 15),
                  axis.text.x=element_text(size=10,angle=0),
                  axis.text.y=element_text(size=10),
                  axis.title.x=element_text(size = 15),
                  axis.title.y=element_text(size = 15)) + 
            labs(x='Groups', y= 'R value') + 
            geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
          ggsave(paste0("./results/Inter-Cor/pairwise/R_value_between_groups_box_plot.pdf"),width = 7,height = 5)
          
          ### Venn Diagram
          ### Substances with |r| > venn.rthreshold && p < venn.pthreshold were screened, record the subscript 
          meta_list <- list()
          mic_list <- list()
          for (i in 1:length(unique(phenoData[,1]))){
            temp_r_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_p.adjust_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]]
            meta_df <- c()
            mic_df <- c()
            for (j in 1:nrow(temp_r_df)){
              for (k in 1:ncol(temp_r_df)){
                if(abs(temp_r_df[j,k]) >= venn.rthreshold && temp_p.adjust_df[j,k] < venn.pthreshold){
                  meta_df <- append(meta_df,k)
                  mic_df <- append(mic_df,j)
                }
              }
            }
            meta_list[[unique(phenoData[,1])[i]]] = meta_df
            mic_list[[unique(phenoData[,1])[i]]] = mic_df
          }
          if(length(meta_list) != length(unique(phenoData[,1])) || length(mic_list) != length(unique(phenoData[,1]))){
            print(paste0("There was no pair of metabolites and microbes with |r| >",venn.rthreshold," and p-value < ",venn.pthreshold," between different groups!"))
          }else{
            # Meta and MIC subscripts are combined to represent pairs names, 
            # and if the same pair names are present, a correlation pair with an intersection is indicated
            
            pairs_df_vn <- list()
            for (i in 1:length(unique(phenoData[,1]))) {
              for (j in 1:length(mic_list[[i]])) {
                pairs_df_vn[[unique(phenoData[,1])[i]]][j] <- paste0(colnames(group_res_list[[i]][[paste0("r_",phenotype)]])[meta_list[[i]][j]],"-",rownames(group_res_list[[i]][[paste0("r_",phenotype)]])[mic_list[[i]][j]])
              }
              # Output significant correlation pairs for each group
              write.csv(pairs_df_vn[[unique(phenoData[,1])[i]]],paste0("./results/Inter-Cor/pairwise/",unique(phenoData[,1])[i],"_sig_cor_pairs.csv"),row.names = F)
            }
            
            # Output intersections of significant correlation pairs for each group
            inter_pairs <- Reduce(intersect,pairs_df_vn)
            write.csv(inter_pairs,"./results/Inter-Cor/pairwise/inter_pairs.csv")
            
            fill_col <- rainbow(length(unique(phenoData[,1])))
            
            ## pairs_null_flag:Determines whether no element exists
            pairs_null_flag = 0
            for (i in 1:length(pairs_df_vn)) {
              if(length(pairs_df_vn[[i]]) == 0){
                pairs_null_flag = 1
              }
            }
            if(pairs_null_flag == 0){
              venn.plot <- venn.diagram(x = pairs_df_vn,NULL,col = "black",fill = fill_col,
                                        cat.cex = 2,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")
              pdf(file="./results/Inter-Cor/pairwise/vn_plot_of_sig_cor_pairs_between_groups.pdf")
              grid.draw(venn.plot)  # library(grDevices)
              dev.off()
            }
          }
        }
        return(res_list)
      }
      
      res_list.list <- list()
      # correlation for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        if(is.na(confounderData)){
          tempconfounderData <- NA
        }else{
          tempconfounderData <- as.data.frame(confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]) 
          rownames(tempconfounderData) <- rownames(tempmetaData)
          colnames(tempconfounderData) <- colnames(confounderData)
        }
        temp_list <- cate_cor_func(metaInddf = tempmetaData,micInddf = tempmicData,confdf = tempconfounderData, phenotype = i,method = corMethod)
        res_list.list <- append(res_list.list,temp_list)
      }
      
      # correlation for all samples
      temp_list <- cate_cor_func(metaInddf = metaIndData,micInddf = micIndData,confdf= confounderData, phenotype = "All",method = corMethod)
      
      res_list.list <- append(res_list.list,temp_list)
      return(res_list.list)
      
    }
    
    ### Method 4 - Generalized linear models
    if (corMethod == "glm") {
      cate_cor_func <- function(metaInddf,micInddf,confdf,phenotype){
        mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
        dim(mic_meta_cor)
        rownames (mic_meta_cor) = colnames (micInddf)
        colnames (mic_meta_cor) = colnames (metaInddf)
        mic_meta_p = mic_meta_cor
        mic_meta_p.adjust = mic_meta_cor
        for (m in colnames (metaInddf)) {
          # exist confounder
          if(is.na(confdf)){
            mic_meta_cor [ , m] <- apply (micInddf, MARGIN = 2, FUN = function (x) 
              summary(glm(metaInddf[,m] ~ x,family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (micInddf, MARGIN = 2, FUN = function (x) 
              summary(glm(metaInddf[,m] ~ x,family = gaussian))$coefficients[2,4])
          }else{
            temp_formula <- paste0("confdf$",colnames(confdf),collapse = "+")
            mic_meta_cor [ , m] <- apply (micInddf, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("metaInddf[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
            mic_meta_p [ , m] <- apply (micInddf, MARGIN = 2, FUN = function (x) 
              summary(glm(as.formula(paste0("metaInddf[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
          }
        }
        mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
          p.adjust(x,method = "BH"))
        # write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
        # write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
        # write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
        write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_cor_",phenotype,".csv"))
        write.csv(mic_meta_p, paste0("./results/Inter-Cor/pairwise/mic_meta_p_",phenotype,".csv"))
        write.csv(mic_meta_p.adjust, paste0("./results/Inter-Cor/pairwise/mic_meta_p.adjust_",phenotype,".csv"))
        res_list <- list()
        res_list[[paste0("r_",phenotype)]] = as.data.frame(mic_meta_cor)
        res_list[[paste0("p_",phenotype)]] = as.data.frame(mic_meta_p)
        res_list[[paste0("p.adjust_",phenotype)]] = as.data.frame(mic_meta_p.adjust)
        
        if(phenotype == "All"){
          ### Bar plot,box plot in all groups
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 && p<0.05 was counted
          counts <- 0 
          for (i in 1:nrow(mic_meta_cor)) {
            for (j in 1:ncol(mic_meta_cor)) {
              if(abs(mic_meta_cor[i,j]) > 0.3 && mic_meta_p.adjust[i,j] < 0.05){
                counts <- counts + 1
                # print(paste0("i:",i,"j:",j))
              }
            }
          }
          # counts <- sum(abs(mic_meta_cor > 0.5))
          r_df <- data.frame("All groups",counts)
          counts <- sum(mic_meta_p.adjust < 0.05)
          p.adjust_df <- data.frame("All groups",counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df <- rbind(r_df,p.adjust_df)
          colnames(bar_df)[1] <- "group_name"
          
          ### Box plot
          box_df <- mic_meta_cor
          box_df <- as.matrix(box_df)
          dim(box_df) <- c(ncol(box_df) * nrow(box_df),1)
          box_df <- as.data.frame(box_df)
          box_df$group <- "All groups"
          colnames(box_df) <- c("r_value","Groups")
          
          ### Bar plot,box plot and Venn Diagram between each group
          
          ### Identify groups
          temp_phenoData <- phenoData
          temp_phenoData$SampleID <- rownames(temp_phenoData)
          
          temp_metaInddf <- metaInddf
          temp_metaInddf$SampleID <- rownames(temp_metaInddf)
          
          temp_micInddf <- micInddf
          temp_micInddf$SampleID <- rownames(temp_micInddf)
          
          if(is.na(confdf)){
            temp_covdata <- NA
          }else{
            temp_covdata <- confdf
            temp_covdata$SampleID <- rownames(temp_covdata)
            temp_covdata <- merge(temp_covdata,temp_phenoData,all.x = T)
          }
          
          temp_metaInddf <- merge(temp_metaInddf,temp_phenoData,all.x = T)
          temp_micInddf <- merge(temp_micInddf,temp_phenoData,all.x = T)
          
          ### Cor between metabolites and microbes in each group
          unique(phenoData[,1]) 
          group_res_list <- list()
          for(i in 1:length(unique(phenoData[,1]))){
            group_metaInddf <- temp_metaInddf[which(temp_metaInddf[,ncol(temp_metaInddf)]== unique(phenoData[,1])[i]),]
            group_metaInddf <- group_metaInddf[,-1]
            group_metaInddf <- group_metaInddf[,-ncol(group_metaInddf)]
            
            group_micInddf <- temp_micInddf[which(temp_micInddf[,ncol(temp_micInddf)]== unique(phenoData[,1])[i]),]
            group_micInddf <- group_micInddf[,-1]
            group_micInddf <- group_micInddf[,-ncol(group_micInddf)]
            
            if(is.na(confdf)){
              group_confdf <- NA
            }else{
              group_confdf <- temp_covdata[which(temp_covdata[,ncol(temp_covdata)]== unique(phenoData[,1])[i]),]
              group_confdf <- group_confdf[,-1]
              group_confdf <- group_confdf[,-ncol(group_confdf)]
              group_confdf <- as.data.frame(group_confdf)
              rownames(group_confdf) <- rownames(group_micInddf)
            }
            
            mic_meta_cor = matrix (NA, nrow = ncol (micInddf), ncol = ncol (metaInddf))
            dim(mic_meta_cor)
            rownames (mic_meta_cor) = colnames (micInddf)
            colnames (mic_meta_cor) = colnames (metaInddf)
            mic_meta_p = mic_meta_cor
            mic_meta_p.adjust = mic_meta_cor
            for (m in colnames (group_metaInddf)) {
              # exist confounder
              if(is.na(confdf)){
                mic_meta_cor [ , m] <- apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                  summary(glm(group_metaInddf[,m] ~ x,family = gaussian))$coefficients[2,1])
                mic_meta_p [ , m] <- apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                  summary(glm(group_metaInddf[,m] ~ x,family = gaussian))$coefficients[2,4])
              }else{
                temp_formula <- paste0("confdf$",colnames(confdf),collapse = "+")
                mic_meta_cor [ , m] <- apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                  summary(glm(group_metaInddf[,m] ~ x + group_confdf[,1],family = gaussian))$coefficients[2,1])
                mic_meta_p [ , m] <- apply (group_micInddf, MARGIN = 2, FUN = function (x) 
                  summary(glm(group_metaInddf[,m] ~ x + group_confdf[,1],family = gaussian))$coefficients[2,4])
              }
            }

            mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
              p.adjust(x,method = "BH"))
            
            mic_meta_cor <- as.data.frame(mic_meta_cor)
            mic_meta_p <- as.data.frame(mic_meta_p)
            mic_meta_p.adjust <- as.data.frame(mic_meta_p.adjust)
            rownames(mic_meta_cor) = rownames(mic_meta_p) = rownames(mic_meta_p.adjust) = colnames(group_micInddf)
            colnames(mic_meta_cor) = colnames(mic_meta_p) = colnames(mic_meta_p.adjust) = colnames(group_metaInddf)
            
            ### Delete NA
            mic_meta_cor <- na.omit(mic_meta_cor)
            mic_meta_p <- na.omit(mic_meta_p)
            mic_meta_p.adjust <- na.omit(mic_meta_p.adjust)
            
            # write.csv(mic_meta_cor,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_cor.csv"))
            # write.csv(mic_meta_p,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.csv"))
            # write.csv(mic_meta_p.adjust,paste0("./results/Inter-Cor/pairwise/mic_meta_",unique(phenoData[,1])[i],"_p.adjust.csv"))
            
            group_spearman <- list()
            group_spearman[[paste0("r_",phenotype)]]  <- mic_meta_cor
            group_spearman[[paste0("p_",phenotype)]]  <- mic_meta_p
            group_spearman[[paste0("p.adjust_",phenotype)]]  <- mic_meta_p.adjust
            group_res_list[[unique(phenoData[,1])[i]]] <- group_spearman
          }
          
          ### Bar plot
          ### The number of p<0.05 and |r| >0.5 was counted
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            temp_counts <- 0
            for(j in 1: nrow(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])){
              for (k in 1:ncol(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]])) {
                if(abs(group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]][j,k]) > 0.5 && group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]][j,k] < 0.05){
                  temp_counts <- temp_counts + 1
                  # print(paste0("j:",j,"k:",k))
                }
              }
            }
            counts <- append(counts,temp_counts) 
          }
          group_name <- unique(phenoData[,1])
          r_df <- data.frame(group_name,counts)
          
          counts <- c()
          for (i in 1:length(unique(phenoData[,1]))) {
            counts <- append(counts,sum(group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]] < 0.05)) 
          }
          p.adjust_df <- data.frame(group_name,counts)
          
          r_df$counts <- as.numeric(r_df$counts)
          r_df$object <- "|r| > 0.5 and p < 0.05"
          p.adjust_df$counts <- as.numeric(p.adjust_df$counts)
          p.adjust_df$object <- "p < 0.05"
          
          bar_df2 <- rbind(r_df,p.adjust_df)
          
          bar_df <- rbind(bar_df,bar_df2)
          ### Bar plot
          bar_plot <- ggplot(bar_df,aes(group_name,counts))+
            geom_bar(stat="identity",position="dodge",fill = "#00b692")+
            labs(x = 'Groups', y = 'Counts of correlation pairs', size = 2) +
            theme(legend.position='none')+
            facet_grid(object ~ .)+
            coord_flip()
          
          ggsave(paste0("./results/Inter-Cor/pairwise/Counts_between_groups_bar_plot.pdf"),width = 7,height = 5)
          
          ### Box plot
          box_df2 <- group_res_list[[unique(phenoData[,1])[1]]][[paste0("r_",phenotype)]]
          box_df2 <- as.matrix(box_df2)
          dim(box_df2) <- c(ncol(box_df2) * nrow(box_df2),1)
          box_df2 <- as.data.frame(box_df2)
          box_df2$group <- unique(phenoData[,1])[1]
          for (i in 2:length(unique(phenoData[,1]))){
            temp_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_df <- as.matrix(temp_df)
            dim(temp_df) <- c(ncol(temp_df) * nrow(temp_df),1)
            temp_df <- as.data.frame(temp_df)
            temp_df$group <- unique(phenoData[,1])[i]
            box_df2 <- rbind(box_df2,temp_df)
          }
          colnames(box_df2) <- c("r_value","Groups")
          
          box_df <- rbind(box_df,box_df2)
          compaired <- combn(unique(box_df$Groups),2,simplify = F)
          
          box_plot <- ggplot(box_df,aes(Groups,r_value,fill = Groups)) + 
            geom_boxplot(width = 0.6,colour = "black") + 
            theme(plot.title=element_text(size = 15),
                  axis.text.x=element_text(size=10,angle=0),
                  axis.text.y=element_text(size=10),
                  axis.title.x=element_text(size = 15),
                  axis.title.y=element_text(size = 15)) + 
            labs(x='Groups', y= 'R value') + 
            geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
          ggsave(paste0("./results/Inter-Cor/pairwise/R_value_between_groups_box_plot.pdf"),width = 7,height = 5)
          
          ### Venn Diagram
          ### Substances with |r| > venn.rthreshold && p < venn.pthreshold were screened, record the subscript 
          meta_list <- list()
          mic_list <- list()
          for (i in 1:length(unique(phenoData[,1]))){
            temp_r_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("r_",phenotype)]]
            temp_p.adjust_df <- group_res_list[[unique(phenoData[,1])[i]]][[paste0("p.adjust_",phenotype)]]
            meta_df <- c()
            mic_df <- c()
            for (j in 1:nrow(temp_r_df)){
              for (k in 1:ncol(temp_r_df)){
                if(abs(temp_r_df[j,k]) >= venn.rthreshold && temp_p.adjust_df[j,k] < venn.pthreshold){
                  meta_df <- append(meta_df,k)
                  mic_df <- append(mic_df,j)
                }
              }
            }
            meta_list[[unique(phenoData[,1])[i]]] = meta_df
            mic_list[[unique(phenoData[,1])[i]]] = mic_df
          }
          if(length(meta_list) != length(unique(phenoData[,1])) || length(mic_list) != length(unique(phenoData[,1]))){
            print(paste0("There was no pair of metabolites and microbes with |r| >",venn.rthreshold," and p-value < ",venn.pthreshold," between different groups!"))
          }else{
            # Meta and MIC subscripts are combined to represent pairs names, 
            # and if the same pair names are present, a correlation pair with an intersection is indicated
            
            pairs_df_vn <- list()
            for (i in 1:length(unique(phenoData[,1]))) {
              for (j in 1:length(mic_list[[i]])) {
                pairs_df_vn[[unique(phenoData[,1])[i]]][j] <- paste0(colnames(group_res_list[[i]][[paste0("r_",phenotype)]])[meta_list[[i]][j]],"-",rownames(group_res_list[[i]][[paste0("r_",phenotype)]])[mic_list[[i]][j]])
              }
              # Output significant correlation pairs for each group
              write.csv(pairs_df_vn[[unique(phenoData[,1])[i]]],paste0("./results/Inter-Cor/pairwise/",unique(phenoData[,1])[i],"_sig_cor_pairs.csv"),row.names = F)
            }
            
            # Output intersections of significant correlation pairs for each group
            inter_pairs <- Reduce(intersect,pairs_df_vn)
            write.csv(inter_pairs,"./results/Inter-Cor/pairwise/inter_pairs.csv")
            
            fill_col <- rainbow(length(unique(phenoData[,1])))
            
            ## pairs_null_flag:Determines whether no element exists
            pairs_null_flag = 0
            for (i in 1:length(pairs_df_vn)) {
              if(length(pairs_df_vn[[i]]) == 0){
                pairs_null_flag = 1
              }
            }
            if(pairs_null_flag == 0){
              venn.plot <- venn.diagram(x = pairs_df_vn,NULL,col = "black",fill = fill_col,
                                        cat.cex = 2,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")
              pdf(file="./results/Inter-Cor/pairwise/vn_plot_of_sig_cor_pairs_between_groups.pdf")
              grid.draw(venn.plot)  # library(grDevices)
              dev.off()
            }
          }
        }
        return(res_list)
      }
      res_list.list <- list()
      # glm for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micIndData[rownames(phenoData)[which(phenoData[,1] == i)],]
        if(is.na(confounderData)){
          tempconfounderData <- NA
        }else{
          tempconfounderData <- as.data.frame(confounderData[rownames(phenoData)[which(phenoData[,1] == i)],]) 
          rownames(tempconfounderData) <- rownames(tempmetaData)
          colnames(tempconfounderData) <- colnames(confounderData)
        }
        temp_list <- cate_cor_func(metaInddf = tempmetaData,micInddf = tempmicData,confdf = tempconfounderData, phenotype = i)
        res_list.list <- append(res_list.list,temp_list)
      }
      
      # glm for all samples
      temp_list <- cate_cor_func(metaInddf = metaIndData,micInddf = micIndData,confdf= confounderData, phenotype = "All")
      
      res_list.list <- append(res_list.list,temp_list)
      return(res_list.list)
    }
  }
  
  ############ No phenotypes ############ 
  else{
    ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
    if(corMethod == "gramm"){
      naivegramm <- naivegramm(metaIndData,micIndData,covdata = confounderData)
      write.csv(naivegramm[[paste0("r_",phenotype)]],"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(naivegramm[[paste0("p_",phenotype)]], "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      mic_meta_p.adjust = apply(naivegramm[[paste0("p_",phenotype)]], 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      naivegramm[[paste0("p.adjust_",phenotype)]] = as.data.frame(mic_meta_p.adjust)
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      write.csv(naivegramm[["type"]], "./results/Inter-Cor/pairwise/mic_meta_cor_type.csv")
      
      res_list <- naivegramm
      return(res_list)
    }
    
    ### Method 2 - spearman,pearson,kendall
    if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        mic_meta_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
        mic_meta_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      
      res_list <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      return(res_list)
    }
    
    ### Method 3 - partial spearman,partial pearson,partial kendall
    if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
      corMethod = strsplit(corMethod," ")[[1]][2]
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        mic_meta_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          ppcor::pcor.test (x, metaIndData[,m], confounderData,method = corMethod)$estimate)
        mic_meta_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          ppcor::pcor.test (x, metaIndData[,m], confounderData,method = corMethod)$p.value)
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      
      res_list <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      return(res_list)
    }
    
    ### Method 4 - Generalized linear models
    if (corMethod == "glm") {
      mic_meta_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (metaIndData))
      dim(mic_meta_cor)
      rownames (mic_meta_cor) = colnames (micIndData)
      colnames (mic_meta_cor) = colnames (metaIndData)
      mic_meta_p = mic_meta_cor
      mic_meta_p.adjust = mic_meta_cor
      for (m in colnames (metaIndData)) {
        # exist confounder
        if(is.na(confdf)){
          mic_meta_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
        }else{
          temp_formula <- paste0("confdf$",colnames(confdf),collapse = "+")
          mic_meta_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
          mic_meta_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
            summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
        }
      }
      mic_meta_p.adjust = apply(mic_meta_p, 2, FUN = function(x)
        p.adjust(x,method = "BH"))
      write.csv(mic_meta_cor,"./results/Inter-Cor/pairwise/mic_meta_cor.csv")
      write.csv(mic_meta_p, "./results/Inter-Cor/pairwise/mic_meta_p.csv")
      write.csv(mic_meta_p.adjust, "./results/Inter-Cor/pairwise/mic_meta_p.adjust.csv")
      
      res_list <- list(r = as.data.frame(mic_meta_cor), p = as.data.frame(mic_meta_p), p.adjust = as.data.frame(mic_meta_p.adjust))
      return(res_list)
    }
  }
}

