# function - naivegramm #
########################################################################
# File: naivegramm.R
# Aim : Generalized Correlation Analysis for Metabolome and Microbiome (GRaMM), for inter-correlation pairs discovery 
#       among metabolome and microbiome.
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
##    micData --- A dataframe of microbes                                           ##
##                       (rows: samples,columns: microbes)                          ##  
##    covdata --- A dataframe of confounder                                         ##
##                       (rows:samples,columns:confounder)                          ##
##    r --- Correlation coefficient threshold                                       ##
##                Default: TRUE                                                     ##
##    alpha --- FDR threshold                                                       ##
##                Default: 0.99                                                     ##
##    pheno --- Phenotype name                                                      ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    res --- A list included coefficient, p value and type of correlation          ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
naivegramm<-function(metaData, micData, covdata, r = 0.3,alpha = 0.99,pheno="All")
{
  metaData <- as.data.frame(t(metaData))
  micData <- as.data.frame(t(micData))
  covdata <- as.data.frame(t(covdata))
  ## Exist confounders
  if(is.data.frame(covdata) && is.na(covdata[1,1]) != T){
    options(warn = -1)
    
    ## Generate result table
    r_result <- p_result <- r2_result <- type_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    rownames(p_result) <- rownames(r_result) <- rownames(type_result) <- rownames(metaData)
    colnames(p_result) <- colnames(r_result) <- colnames(type_result) <- rownames(micData)
    
    ## lm function was used to obtain the adjusted data and fit the linear model include linear regression and multiple linear regression
    for(i in seq_len(nrow(metaData))){
      for(j in seq_len(nrow(micData))){
        options(warn = -1)
        x1 <- t(metaData[i,])
        y1 <- t(rbind(micData[j,], covdata))
        lmx <- lm(x1 ~ y1)
        x11 <- y1[,1] - lmx$coefficients[3] * y1[,2]
        x2 <- t(rbind(metaData[i,], covdata))
        y2 <- t(micData[j,])
        lmy <- lm(y2 ~ x2)
        y22 <- x2[,1] - lmy$coefficients[3] * x2[,2]
        summ <- summary(lmx)
        pvalue <- summ$coefficients[2,4]
        rvalue <- lmx$coefficients[2] * sd(y1[,1]) / sd(x1)
        r2value <- rvalue^2
        
        ## linear: multiple linear regression
        if(pvalue < alpha & rvalue > r & !is.na(pvalue) & !is.na(rvalue)){
          p_result[i,j] <- pvalue
          r_result[i,j] <- rvalue
          r2_result[i,j] <- r2value
          type_result[i,j] <- "linear"
        }
        
        ## nonlinear: MIC
        else{
          micxy <- minerva::mine(y22,x11)
          micr <- micxy$MIC
          micp <- 0
          for(k in seq_len(101))
          {
            #bootstrap
            options(warn = -1)
            bootx <- matrix(sample(x11, replace = TRUE))
            booty <- matrix(sample(y22, replace = TRUE))
            if(sd(bootx) == 0){
              bootx <- matrix(sample(x11, replace = FALSE))
            }
            if(sd(booty) == 0){
              booty <- matrix(sample(y22, replace = FALSE))
            }
            tmp <- minerva::mine(booty, bootx)
            MICtp <- tmp$MIC
            tempM <- ifelse(micr <= MICtp,1,0)
            micp <- tempM + micp
          }
          micp <- micp / 101
          prsmic <- cor.test(y22, x11, method="pearson",use = "pairwise.complete.obs")
          prsr <- prsmic$estimate
          ltmic <- 1-micr + prsr^2
          p_result[i,j] <- micp
          r_result[i,j] <- micr
          #r2_result[i,j]<-ltmic
          type_result[i,j] <- "nonlinear"
        }

      }

    }
  }
  
  ## No confounders
  if(is.na(covdata[1,1]) == T){
    options(warn = -1)
    r_result <- p_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    r2_result <- type_result <- matrix(0, nrow = nrow(metaData),ncol = nrow(micData))
    rownames(p_result) <- rownames(r_result) <- rownames(type_result) <- rownames(metaData)
    colnames(p_result) <- colnames(r_result) <- colnames(type_result) <- rownames(micData)

    ## lm function was used to obtain the adjusted data and fit the linear model include linear regression and multiple linear regression
    for(i in seq_len(nrow(metaData))){
      for(j in seq_len(nrow(micData))){
        x1 <- t(metaData[i,])
        y1 <- t(micData[j,])
        lmx <- lm(x1 ~ y1)
        summ <- summary(lmx)
        pvalue <- summ$coefficients[2,4]
        rvalue <- lmx$coefficients[2] * sd(y1) / sd(x1)
        
        ## linear: linear regression
        if(pvalue < alpha & rvalue > r){
          p_result[i,j] <- pvalue
          r_result[i,j] <- rvalue
          r2value <- rvalue^2
          r2_result[i,j] <- r2value
          type_result[i,j] <- "linear"
        }
        
        ## nonlinear: MIC
        else{
          warnings('off')
          micxy <- minerva::mine(y1,x1)
          micr <- micxy$MIC
          micp <- 0
          for(k in seq_len(101))
          {
            bootx <- matrix(sample(x1, replace = FALSE))
            booty <- matrix(sample(y1, replace = FALSE))
            tmp <- minerva::mine(bootx, booty)
            MICtp <- tmp$MIC
            tempM <- ifelse(micr <= MICtp,1,0)
            micp <- tempM + micp
          }
          micp <- micp/101
          prsmic <- cor.test(y1,x1,method="pearson",use = "pairwise.complete.obs")
          prsr <- prsmic$estimate
          ltmic <- 1 - micr + prsr^2
          p_result[i,j] <- micp
          r_result[i,j] <- micr
          type_result[i,j] <- "nonlinear"
        }
      }
    }
  }
  res <- list()
  res[[paste0("r_",pheno)]] <- as.data.frame(t(r_result))
  res[[paste0("p_",pheno)]] <- as.data.frame(t(p_result))
  res[[paste0("type_",pheno)]] <- as.data.frame(t(type_result))
  return(res)
}

