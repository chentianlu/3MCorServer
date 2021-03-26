# function - microbesPrepro #
########################################################################
# File: microbesPrepro.R
# Aim : Preprocessing for microbiome data
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
##    missPro --- Flag of whether to complete the missing value                     ##
##                Default: TRUE                                                     ##
##    missMethod --- Method for completing missing values                           ##
##                Default: mean                                                     ##
##    rarePro --- Flag of whether to rarefy                                         ##
##                Default: TRUE                                                     ##
##    scalingPro --- Flag of whether to scale                                       ##
##                Default: TRUE                                                     ##
##    transPro --- Flag of whether to transform                                     ##
##                Default: TRUE                                                     ##
##    kValue --- The number of nearest neighbours to use in knn imputation          ##
##                Default: 3                                                        ##
##    phenoData --- A dataframe of phenotype data                                   ##
##                Default: NA  (rows: samples,columns: phenotype)                   ##
##    phenoDataType --- Phenotype data type                                         ##
##                Default: continuous                                               ##
######################################################################################
## Output:                                                                          ## 
##     microbes_preproRes --- The dataframe of microbes after preprocessing         ##
##                       (rows: samples,columns: microbes)                          ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------

microbesPrepro <- function(micData, missPro = TRUE,missMethod = "mean",rarePro = T, scalingPro = T, transPro = T, kValue = 3,phenoData = NA,phenoDataType = "continuous")
{
  if(!file.exists("./results/prepro/")){
    dir.create("./results/prepro/",recursive = T)
  }
  
  microbes_preproRes <- micData
  
  microbes <- micData
  
  if(missPro){
    ######################################
    ############ Step 1 - knn ############
    ######################################
    if(missMethod == "knn"){

      microbes_preproRes <- DMwR::knnImputation(microbes_preproRes,k = kValue,scale = TRUE)
    }else if(missMethod == "QRILC"){
      obj.QRILC = imputeLCMD::impute.QRILC(microbes_preproRes)
      microbes_preproRes = obj.QRILC[[1]]
    }else if(missMethod == "mean"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = mean)
      }
    }else if(missMethod == "min"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = min)
      }
    }else if(missMethod == "median"){
      for (i in colnames(microbes_preproRes)) {
        microbes_preproRes[,i] <- Hmisc::impute(microbes_preproRes[,i], fun = median)
      }
    } 
  }
  
  if(phenoDataType != "categorical"){
    # If the value of 0 is more than 70%, the microbe is deleted 
    index <- c()
    for (i in 1:ncol(microbes_preproRes)) {
      sum_count <- sum(microbes_preproRes[,i] == 0)
      if(sum_count > (nrow(microbes_preproRes)*0.7) ){
        index <- append(index,i)
      }
    }
    if(!is.null(index)){
      microbes_preproRes <- microbes_preproRes[,-index]
    }
  }else if(phenoDataType == "categorical"){
    if(is.na(phenoData)){
      stop("Please input your phenoData!")
    }else{
      index <- c()
      for (i in unique(phenoData[,1])) {
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        # If the value of 0 is more than 70% in each group, the microbe is deleted
        for (j in 1:ncol(tempmicData)) {
          sum_count <- sum(tempmicData[,j] == 0)
          if(sum_count > (nrow(tempmicData)*0.7) ){
            index <- append(index,j)
          }
        }
      }
      if(!is.null(index)){
        index <- unique(index)
        microbes_preproRes <- microbes_preproRes[,-index]
      }
    }
  }
  
  if(rarePro){
    
    ######################################
    ######## Step 1 - rarefaction ########
    ######################################
    
    microbes_preproRes <- apply(microbes,2,as.integer)
    microbes_preproRes <- t(microbes_preproRes)
    colnames(microbes_preproRes) <- rownames(microbes)
    microbes_number <- nrow(microbes_preproRes) * ncol(microbes_preproRes)
    
    # Create a pretend taxonomy table
    taxmat <- matrix(sample(letters, microbes_number, replace = TRUE), nrow = nrow(microbes_preproRes),
                     ncol = 7)
    
    rownames(taxmat) <- rownames(microbes_preproRes)
    colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family",
                          "Genus", "Species")
    OTU = phyloseq::otu_table(microbes_preproRes, taxa_are_rows = TRUE)
    TAX = phyloseq::tax_table(taxmat)
    physeq = phyloseq(OTU, TAX)
    ## set.seed(100)
    raretest <- phyloseq::rarefy_even_depth(physeq,rngseed=100)
    microbes_preproRes <- raretest@otu_table@.Data
    microbes_preproRes <- t(microbes_preproRes)
    colnames(microbes_preproRes) <- colnames(microbes)
    rownames(microbes_preproRes) <- rownames(microbes)
    
    ##  deal with 0 value
    microbes_preproRes[microbes_preproRes == 0] <- 1
  }
  
  if(scalingPro){
    ######################################
    ####### Step 2 - normalization #######
    ######################################

    microbes_preproRes<- apply(microbes_preproRes,2,function(each_col){
      col_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col / col_sum * 30000
      return(each_col)
    })
  }
    
  if(transPro){
    ######################################
    ############ Step 3 - clr ############
    ######################################

    microbes_preproRes <- apply(microbes_preproRes,2,function(each_col){
      log_col <- psych::geometric.mean(each_col)
      each_col <- abs(log(each_col) / log(log_col))
      return(each_col)
    })
  } 
    
  microbes_preproRes <- as.data.frame(microbes_preproRes)
  write.csv(microbes_preproRes,"./results/prepro/microbes_prepro_res.csv")
  return(microbes_preproRes)

}
