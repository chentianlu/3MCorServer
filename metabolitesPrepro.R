# function - metabolitesPrepro #
########################################################################
# File: metabolitesPrepro.R
# Aim : Preprocessing for metabolome data
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
##    metaData --- A dataframe of metabolites                                       ##
##                       (rows: samples,columns: metabolites)                       ##
##    missPro --- Flag of whether to complete the missing value                     ##
##                Default: TRUE                                                     ##
##    missMethod --- Method for completing missing values                           ##
##                Default: QRILC                                                    ##
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
##    metabolites_preproRes --- The dataframe of metabolites after preprocessing    ##
##                       (rows: samples,columns: metabolites)                       ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
metabolitesPrepro <- function(metaData, missPro = TRUE, missMethod = "QRILC", scalingPro = TRUE, transPro = TRUE, kValue = 3,phenoData = NA,phenoDataType = "continuous")  
{
  if(!file.exists("./results/prepro")){
    dir.create("./results/prepro",recursive = T)
  }

  metabolites_preproRes <- metaData

  ######################################
  ############ missing value ###########
  ######################################
  if(missPro){
    if(missMethod == "knn"){
      # metabolites_preproRes[metabolites_preproRes == 0] <- NA
      metabolites_preproRes <- DMwR::knnImputation(metabolites_preproRes,k = kValue,scale = F)
    }else if(missMethod == "QRILC"){
      obj.QRILC = imputeLCMD::impute.QRILC(metabolites_preproRes)
      metabolites_preproRes = obj.QRILC[[1]]
    }else if(missMethod == "mean"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = mean)
      }
    }else if(missMethod == "min"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = min)
      }
    }else if(missMethod == "median"){
      for (i in colnames(metabolites_preproRes)) {
        metabolites_preproRes[,i] <- Hmisc::impute(metabolites_preproRes[,i], fun = median)
      }
    } 
  }
  
  if(phenoDataType != "categorical"){
    # If the value of 0 is more than 70%, the metabolite is deleted
    index <- c()
    for (i in 1:ncol(metabolites_preproRes)) {
      sum_count <- sum(metabolites_preproRes[,i] == 0)
      if(sum_count > (nrow(metabolites_preproRes)*0.7) ){
        index <- append(index,i)
      }
    }
    if(!is.null(index)){
      metabolites_preproRes <- metabolites_preproRes[,-index]
    }
  }else if(phenoDataType == "categorical"){
    if(is.na(phenoData)){
      stop("Please input your phenoData!")
    }else{
      index <- c()
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        # If the value of 0 is more than 70% in each group, the metabolite is deleted
        for (j in 1:ncol(tempmetaData)) {
          sum_count <- sum(tempmetaData[,j] == 0)
          if(sum_count > (nrow(tempmetaData)*0.7) ){
            index <- append(index,j)
          }
        }
      }
      if(!is.null(index)){
        index <- unique(index)
        metabolites_preproRes <- metabolites_preproRes[,-index]
      }
    }
  }
  
  ######################################
  ############ normalization ###########
  ######################################
  if(scalingPro){
    metabolites_preproRes<- apply(metabolites_preproRes,2,function(each_col){
      col_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col / col_sum * 30000
      return(each_col)
    })
  }
  
  ######################################
  ########## transformation ############
  ######################################
  if(transPro){
    metabolites_preproRes <- log(metabolites_preproRes + 1)
  }
  
  metabolites_preproRes <- as.data.frame(metabolites_preproRes) 
  write.csv(metabolites_preproRes,"./results/prepro/metabolites_prepro_res.csv")
  return(metabolites_preproRes)
}
