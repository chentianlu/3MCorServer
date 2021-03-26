# function - globalCor #
########################################################################
# File: globalCor.R
# Aim : The global correlation between metabolome and microbiome
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
##    corMethod --- Method of the global correlation                                ##
##                  Default: CCA                                                    ##
##    phenoData --- A dataframe of phenotype data                                   ##
##                Default: NA  (rows: samples,columns: phenotype)                   ##
##    phenoDataType --- Phenotype data type                                         ##
##                Default: continuous                                               ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    Graph results in PDF format                                                   ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
globalCor <- function(metaData,micData,corMethod = "CCA",phenoData = NA,phenoDataType = "continuous")
{
  if(!file.exists("./results/Inter-Cor/global")){
    dir.create("./results/Inter-Cor/global",recursive = T)
  }

  ### CCA
  if(corMethod == "CCA"){
    ##The significance test of canonical correlation coefficient
    # corcoef.test <- function(r, n, p, q, alpha = 0.05){
    #   m <- length(r)
    #   Q <- rep(0, m)
    #   lambda <- 1
    #   for (k in m:1)
    #   {
    #     lambda <- lambda * (1 - r[k]^2);
    #     Q[k] <- -log(lambda)
    #   }
    #   s <- 0
    #   i <- m
    #   temp_chi <- c()
    #   for (k in 1:m)
    #   {
    #     Q[k] <- (n - k + 1 - 1 / 2 * (p + q + 3) + s) * Q[k]
    #     chi <- 1 - pchisq(Q[k], (p - k + 1) *( q - k + 1))
    #     temp_chi <- append(temp_chi,chi)
    #     if (chi > alpha)
    #     {
    #       i <- k - 1
    #       break
    #     }
    #     s <- s + 1 / r[k]^2
    #   }
    #   if(i == 0){
    #     i = 1
    #   }
    #   res <- list(r = cca$cor[i],p = temp_chi[i],i = i)
    #   return(res)
    # }
    
    cca_func <- function(metadf,micdf,phenotype){
      # Deletes rows that sum to 0
      meta_temp_i <- c()
      for (i in 1:nrow(metadf)) {
        if(sum(metadf[i,]) == 0){
          meta_temp_i <- append(meta_temp_i,i)
        }
      }
      if(!is.null(meta_temp_i)){
        metadf <- metadf[-meta_temp_i,]
        micdf <- micdf[-meta_temp_i,]
      }
      
      # Deletes rows that sum to 0
      mic_temp_i <- c()
      for (i in 1:nrow(micdf)) {
        if(sum(micdf[i,]) == 0){
          mic_temp_i <- append(mic_temp_i,i)
        }
      }
      if(!is.null(mic_temp_i)){
        metadf <- metadf[-mic_temp_i,]
        micdf <- micdf[-mic_temp_i,]
      }

      CCAdata1_Z <- metadf
      CCAdata2_Z <- micdf
      
      # # Cancor is used for calculation and output of canonical correlation analysis
      # cca <- cancor(CCAdata1_Z,CCAdata2_Z)

      # # canonical correlation coefficient results
      # # n represents the number of training data samples, P represents the number of indicators in the first group, and Q represents the number of indicators in the second group
      # corcoef_test <- corcoef.test(r = cca$cor,n = nrow(metadf),p = ncol(metadf),q = ncol(micdf))
      
      cca1 <- cca(CCAdata2_Z,CCAdata1_Z)
      cca2 <- cca(CCAdata1_Z,CCAdata2_Z)
      cca_res1 <- summary(cca1, scaling = 1)
      cca_res2 <- summary(cca2, scaling = 1)
      corcoef_test <- cor.test(cca_res1$sites[,1],cca_res2$sites[,1])
      
      ### biplot in CCA1 and CCA2 
      metacoef <- data.frame(CCA1 = cca_res2$species[,1],CCA2 = cca_res2$species[,2])
      rownames(metacoef) <- rownames(cca_res2$species)
      metacoef$color <- "#6300ea"
      metacoef$type <- "Metabolites"
      metacoef$sq <- sqrt(metacoef[,1]^2 + metacoef[,2]^2)
      ### Top 10 Metabolites and Microbes
      if(nrow(metacoef) > 10 ){
        metacoef =  metacoef[order(metacoef$sq,decreasing = T),][1:10,] #按第一列递增排序
      }
      
      miccoef <- data.frame(CCA1 = cca_res1$species[,1],CCA2 = cca_res1$species[,2])
      rownames(miccoef) <- rownames(cca_res1$species)
      miccoef$color <- "#00b692"  
      miccoef$type <- "Microbes"
      miccoef$sq <- sqrt(miccoef[,1]^2 + miccoef[,2]^2)
      if(nrow(miccoef) > 10 ){
        miccoef =  miccoef[order(miccoef$sq,decreasing = T),][1:10,] #按第一列递增排序
      }
      
      biplot_df <- rbind(metacoef,miccoef)
      colnames(metacoef) <- c("CCA1","CCA2","color","Type","sq")
      colnames(miccoef) <- c("CCA1","CCA2","color","Type","sq")
      colnames(biplot_df) <- c("CCA1","CCA2","color","Type","sq")

      cols=c('Metabolites' = "#6300ea", 'Microbes' = "#00b692")
      shapes=c('Metabolites' = 17,'Microbes' = 16)
        
      p <- ggplot(biplot_df, aes(CCA1, CCA2))+
        geom_point(data = metacoef, size = 4, aes(colour = "Metabolites", shape = "Metabolites") )+
        geom_point(data = miccoef, size = 4, aes(colour = "Microbes", shape = "Microbes")) +
        scale_colour_manual(name="Type",values = cols) +
        scale_shape_manual(name="Type",values = shapes) +
        ggtitle(paste0("r:",signif(corcoef_test$estimate,3)," p:",signif(corcoef_test$p.value,3)))+
        xlab("CCA1")+ylab("CCA2") +
        # geom_segment(data = biplot_df,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),color= biplot_df$color,
        #              arrow = arrow(length = unit(0.2, "cm")),size = 1.2) +
        # geom_text(data = biplot_df,aes(x = CCA1, y = CCA2,label = rownames(biplot_df)),size = 3) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) 
      
      ggsave(paste0("./results/Inter-Cor/global/CCA_plot_",phenotype,".pdf"),width = 7,height = 5)
      
      ### importance of metabolites and microbes in CCA1
      meta_importance <- data.frame(importance = cca_res2$species[,1])
      mic_importance <- data.frame(importance = cca_res1$species[,1])
      write.csv(meta_importance,paste0("./results/Inter-Cor/global/CCA1_metabolites_importance_",phenotype,".csv"))
      write.csv(mic_importance,paste0("./results/Inter-Cor/global/CCA1_microbes_importance_",phenotype,".csv"))
      
      ### bar plot 
      meta_plot_data <- meta_importance
      mic_plot_data <- mic_importance
      meta_name <- rownames(meta_plot_data)
      mic_name <- rownames(mic_plot_data)
      meta_df <- data.frame('name' = meta_name,'importance' = meta_plot_data$importance)
      mic_df <- data.frame('name' = mic_name,'importance' = mic_plot_data$importance)
      
      ### Top 10 metabolites or microbes
      meta_df$importance <- abs(meta_df$importance)
      meta_df <- meta_df[order(meta_df$importance,decreasing = F),]
      meta_df$importance <- as.numeric(meta_df$importance)
      
      if(nrow(meta_plot_data) > 10){
        meta_df <- meta_df[(nrow(meta_df)-9):nrow(meta_df),]
      }
      
      meta_df$name <- factor(meta_df$name,levels = meta_df$name)
      
      mic_df$importance <- abs(mic_df$importance)
      mic_df <- mic_df[order(mic_df$importance,decreasing = F),]
      mic_df$importance <- as.numeric(mic_df$importance)
      if(nrow(mic_plot_data) > 10){
        mic_df <- mic_df[(nrow(mic_df)-9):nrow(mic_df),]
      }
      mic_df$name <- factor(mic_df$name,levels = mic_df$name)
      
      ### plot
      meta_plot <- ggplot(data = meta_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#6300ea") + 
        labs(title = paste0("importance of metabolites in CCA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(meta_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/CCA1_metabolites_importanceRank_",phenotype,".pdf"),width = 7,height = 5 )
      
      mic_plot <- ggplot(data = mic_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#00b692") + 
        labs(title = paste0("importance of microbes in CCA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),  
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(mic_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/CCA1_microbes_importanceRank_",phenotype,".pdf"),width = 7,height = 5)
    }
    
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # cca for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        cca_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # cca for all samples
      cca_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      cca_func(metadf = metaData,micdf = micData,phenotype = "All")
    }

  }

  
  ### CIA
  else if (corMethod == "CIA"){
    cia_func <- function(metadf,micdf,phenotype){
      dudi1 <- dudi.pca(metadf, scale = T, scan = FALSE, nf = 6)
      dudi2 <- dudi.pca(micdf, scale = T, scan = FALSE, nf = 6)
      
      #Check to see if the two-step sort has the same row weight, which can only be used for fusion if TRUE
      all.equal(dudi1$lw, dudi2$lw)
      #PCA-PCA fusion of CoIA
      coin <- coinertia(dudi1, dudi2, scan = FALSE, nf = 5)
      # summary(coin) 
      
      #The displacement test determined the significance of CoIA axis and showed that it was significant
      rand_test <- randtest(coin, nrepet = 999)
      
      {
        # summary_coin <- summary(coin)
        # #Eigenvalues of each axis
        # coin$eig
        # #Contribution of eigenvalue of each axis
        # coin$eig / sum(coin$eig)
        # #dudi1 object coordinates
        # coin$co
        # #dudi2 object coordinates
        # coin$li
      }
      
      metaco <- as.data.frame(coin$co[,1:2])
      metaco$color <- "#6300ea"
      metaco$type <- "Metabolites"
      metaco$sq <- sqrt(metaco[,1]^2 + metaco[,2]^2)
      ### Top 10 Metabolites and Microbes
      if(nrow(metaco) > 10 ){
        metaco =  metaco[order(metaco$sq,decreasing = T),][1:10,] #按第一列递增排序
      }
      colnames(metaco) <- c("Axis1","Axis2","color","type","sq")
      
      micli <- as.data.frame(coin$li[,1:2])
      micli$color <- "#00b692"
      micli$type <- "Microbes"
      micli$sq <- sqrt(micli[,1]^2 + micli[,2]^2)
      if(nrow(micli) > 10 ){
        micli =  micli[order(micli$sq,decreasing = T),][1:10,] #按第一列递增排序
      }
      
      biplot_df <- rbind(metaco,micli)
      colnames(biplot_df) <- c("CIA1","CIA2","color","Type","sq")
      colnames(metaco) <- c("CIA1","CIA2","color","Type","sq")
      colnames(micli) <- c("CIA1","CIA2","color","Type","sq")
      
      cols=c('Metabolites' = "#6300ea", 'Microbes' = "#00b692")
      shapes=c('Metabolites' = 17,'Microbes' = 16)
      
      p <- ggplot(biplot_df, aes(CIA1, CIA2))+
        geom_point(data = metaco, size = 4, aes(colour = "Metabolites", shape = "Metabolites") )+
        geom_point(data = micli, size = 4, aes(colour = "Microbes", shape = "Microbes")) +
        scale_colour_manual(name="Type",values = cols) +
        scale_shape_manual(name="Type",values = shapes) +
        ggtitle(paste0("RV:",signif(coin$RV,3)," p:",signif(rand_test$pvalue,3)))+
        xlab("CIA1")+ylab("CIA2") +
        # geom_segment(data = biplot_df,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),color= biplot_df$color,
        #              arrow = arrow(length = unit(0.2, "cm")),size = 1.2) +
        # geom_text(data = biplot_df,aes(x = CCA1, y = CCA2,label = rownames(biplot_df)),size = 3) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) 
      
      ggsave(paste0("./results/Inter-Cor/global/CIA_plot_",phenotype,".pdf"),width = 7,height = 5)
      #CoIA 分析协惯量矩阵特征根的分解情况
      # summary_coin$EigDec
      
      ### importance of metabolites and microbes in Axis1
      meta_importance <- as.data.frame(coin$co[,1])
      colnames(meta_importance) <- "importance"
      rownames(meta_importance) <- rownames(coin$co)
      mic_importance <- as.data.frame(coin$li[,1])
      colnames(mic_importance) <- "importance"
      rownames(mic_importance) <- rownames(coin$li)
      write.csv(meta_importance,paste0("./results/Inter-Cor/global/CIA1_metabolites_importance_",phenotype,".csv"))
      write.csv(mic_importance,paste0("./results/Inter-Cor/global/CIA1_microbes_importance_",phenotype,".csv"))
      
      ### bar plot 
      meta_plot_data <- meta_importance
      mic_plot_data <- mic_importance
      meta_name <- rownames(meta_plot_data)
      mic_name <- rownames(mic_plot_data)
      meta_df <- data.frame('name' = meta_name,'importance' = meta_plot_data$importance)
      mic_df <- data.frame('name' = mic_name,'importance' = mic_plot_data$importance)
      
      ### Top 10 metabolites or microbes
      meta_df$importance <- abs(meta_df$importance)
      meta_df <- meta_df[order(meta_df$importance,decreasing = F),]
      meta_df$importance <- as.numeric(meta_df$importance)
      
      if(nrow(meta_plot_data) > 10){
        meta_df <- meta_df[(nrow(meta_df)-9):nrow(meta_df),]
      }
      
      meta_df$name <- factor(meta_df$name,levels = meta_df$name)
      
      mic_df$importance <- abs(mic_df$importance)
      mic_df <- mic_df[order(mic_df$importance,decreasing = F),]
      mic_df$importance <- as.numeric(mic_df$importance)
      if(nrow(mic_plot_data) > 10){
        mic_df <- mic_df[(nrow(mic_df)-9):nrow(mic_df),]
      }
      mic_df$name <- factor(mic_df$name,levels = mic_df$name)
      
      ### plot
      
      meta_plot <- ggplot(data = meta_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#6300ea") + 
        labs(title = paste0("importance of metabolites in CIA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),  
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(meta_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/CIA1_metabolites_importanceRank_",phenotype,".pdf"),width = 7,height = 5 )
      
      mic_plot <- ggplot(data = mic_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#00b692") + 
        labs(title = paste0("importance of microbes in CIA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) + 
        scale_y_continuous(limits = c(0,as.integer(max(mic_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/CIA1_microbes_importanceRank_",phenotype,".pdf"),width = 7,height = 5)
    }
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # cia for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        cia_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # cia for all samples
      cia_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      cia_func(metaData,micData,phenotype = "All")
    }
  }
  
  ### PA
  else if (corMethod == "PA"){
    pa_func <- function(metadf,micdf,phenotype){
      #The PCA of the environment variable needs to be standardized
      meta_pca <- rda(metadf, scale = TRUE)
      #Hellinger pre-transformation of species data
      mic_hel <- decostand(micdf, method = 'hellinger')
      
      #PCA is performed on transformed species data without standardization
      mic_pca <- rda(mic_hel, scale = FALSE)
      
      #Procrustes analysis
      #Quadrature sorting coordinates in two PCA were extracted, both taking i-scale as an example
      site_meta <- summary(meta_pca, scaling = 1)$site
      site_mic <- summary(mic_pca, scaling = 1)$site
      
      #Perform Procrustes analysis
      #Take symmetry analysis for example(symmetric = TRUE)
      proc <- procrustes(X = meta_pca, Y = mic_pca, symmetric = TRUE)
      summary(proc)
      
      head(proc$Yrot)  #The Y coordinates of Procrustes are analyzed
      head(proc$X)  #The X coordinates of Procrustes are analyzed
      proc$ss  #Deviation squared and M2 statistics
      proc$rotation  #The coordinate position of the axis of rotation can be obtained through this value
      
      #PROTEST test
      #Take 999 permutations
      #Note: The symmetric Procrustes analysis is performed in Protest (), and the allocation swap of X and Y does not affect the calculation of the M2 statistic
      set.seed(123)
      prot <- protest(X = meta_pca, Y = mic_pca, permutations = how(nperm = 999))
      prot
      
      #Extraction of important statistics
      names(prot)
      prot$signif  #p
      prot$ss  #Deviation squared and M2 statistics
      
      Yrot <- data.frame(proc$Yrot)
      Yrot$type <- "Microbes"
      X <- data.frame(proc$X)
      X$type <- "Metabolites"
      colnames(Yrot) <- colnames(X)
      Y <- rbind(X, Yrot)
      arrow_data <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
      
      v1 <- Y[Y$type == "Metabolites",1:2]
      v2 <- Y[Y$type == "Microbes",1:2]
      
      cols=c('Metabolites' = "#6300ea", 'Microbes' = "#00b692")
      shapes=c('Metabolites' = 17,'Microbes' = 16)
      p <- ggplot(Y, aes(PC1, PC2)) +
        geom_segment(data = arrow_data,aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')),
                     color = '#696969', size = 0.8) +
        geom_point(data = v1, size = 4, aes(colour = "Metabolites", shape = "Metabolites") )+
        geom_point(data = v2, size = 4, aes(colour = "Microbes", shape = "Microbes")) +
        scale_colour_manual(name="Type",values = cols) +
        scale_shape_manual(name="Type",values = shapes) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
        labs(x = 'PC1', y = 'PC2', color = '') +
        geom_vline(xintercept = 0, color = 'black', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'black', linetype = 2, size = 0.3) +
        ggtitle(paste0("M^2:",signif(prot$ss,3),"  p:",signif(prot$signif,3)))+
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) 
      ggsave(paste0('./results/Inter-Cor/global/PA_plot_',phenotype,'.pdf'), p, width = 7, height = 5)
      
      
      ### importance of metabolites and microbes in Axis1
      meta_importance <- as.data.frame(summary(meta_pca, scaling = 1)$species[,1])
      colnames(meta_importance) <- "importance"
      rownames(meta_importance) <- names(meta_pca$colsum)
      mic_importance <- as.data.frame(summary(mic_pca, scaling = 1)$species[,1])
      colnames(mic_importance) <- "importance"
      rownames(mic_importance) <- names(mic_pca$colsum)
      write.csv(meta_importance,paste0("./results/Inter-Cor/global/PA1_metabolites_importance_",phenotype,".csv"))
      write.csv(mic_importance,paste0("./results/Inter-Cor/global/PA1_microbes_importance_",phenotype,".csv"))
      
      ### bar plot 
      meta_plot_data <- meta_importance
      mic_plot_data <- mic_importance
      meta_name <- rownames(meta_plot_data)
      mic_name <- rownames(mic_plot_data)
      meta_df <- data.frame('name' = meta_name,'importance' = meta_plot_data$importance)
      mic_df <- data.frame('name' = mic_name,'importance' = mic_plot_data$importance)
      
      ### Top 10 metabolites or microbes
      meta_df$importance <- abs(meta_df$importance)
      meta_df <- meta_df[order(meta_df$importance,decreasing = F),]
      meta_df$importance <- as.numeric(meta_df$importance)
      
      if(nrow(meta_plot_data) > 10){
        meta_df <- meta_df[(nrow(meta_df)-9):nrow(meta_df),]
      }
      
      meta_df$name <- factor(meta_df$name,levels = meta_df$name)
      
      mic_df$importance <- abs(mic_df$importance)
      mic_df <- mic_df[order(mic_df$importance,decreasing = F),]
      mic_df$importance <- as.numeric(mic_df$importance)
      if(nrow(mic_plot_data) > 10){
        mic_df <- mic_df[(nrow(mic_df)-9):nrow(mic_df),]
      }
      mic_df$name <- factor(mic_df$name,levels = mic_df$name)
      
      ### plot
      
      meta_plot <- ggplot(data = meta_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#6300ea") + 
        labs(title = paste0("importance of metabolites in PA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(meta_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/PA1_metabolites_importanceRank_",phenotype,".pdf"),width = 7,height = 5 )
      
      mic_plot <- ggplot(data = mic_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#00b692") + 
        labs(title = paste0("importance of microbes in PA1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(mic_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/PA1_microbes_importanceRank_",phenotype,".pdf"),width = 7,height = 5)
    }
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # pa for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        pa_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # pa for all samples
      pa_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      pa_func(metaData,micData,phenotype = "All")
    }
  }
  
  ### O2PLS
  else if (corMethod == "O2PLS"){
    o2pls_func <- function(metadf,micdf,phenotype){
      X = metadf
      Y = micdf
      
      # X = scale(X, scale = T)
      # Y = scale(Y, scale = T)
      # try(gplots::heatmap.2(cor(X,Y), Rowv=F,Colv=F, col=gplots::bluered(100)),
      # silent = TRUE)
      
      set.seed(1221L)
      crossval_o2m(X, Y, 1:4, 1:2, 1:2, nr_folds = 10)
      fit0 = o2m(X = X, Y = Y, n = 2, nx = 2, ny = 2)
      ### Joint X importance：fit0$Tt
      ### 
      ### Joint X loadings: fit0$W.
      ### 
      ### Joint Y importance：fit0$U
      ### 
      ### Joint Y loadings: fit0$C.
      
      ### Plot fit curve
      meta_joint1 <- as.vector(fit0$Tt[,1])
      mic_joint1 <- as.vector(fit0$U[,1])
      r <- cor.test(meta_joint1,mic_joint1,method = "pearson")
      
      joint1 <- as.data.frame(cbind(meta_joint1,mic_joint1))
      # colnames(joint1) <- c("Metablotes Joint PC1","Microbes Joint PC1")
      p <- ggplot(data = joint1, aes(x = meta_joint1, y = mic_joint1)) +
        geom_point(color = "#d7191c",size = 3) +
        geom_smooth(method = "lm",color = "#00b692") +
        ggtitle(paste0("r:",signif(r$estimate,3),"  p:",signif(r$p.value,3)))+
        xlab("Metabolites joint PC1") + ylab("Microbes joint PC1") +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) 
      ggsave(paste0("./results/Inter-Cor/global/O2PLS_plot_",phenotype,".pdf"),width = 7,height = 5)
      
      ### importance of metabolites and microbes in Axis1
      meta_importance <- as.data.frame(as.vector(fit0$W.[,1]))
      colnames(meta_importance) <- "importance"
      rownames(meta_importance) <- names(fit0$W.[,1])
      mic_importance <- as.data.frame(as.vector(fit0$C.[,1]))
      colnames(mic_importance) <- "importance"
      rownames(mic_importance) <- names(fit0$C.[,1])
      write.csv(meta_importance,paste0("./results/Inter-Cor/global/O2PL1_metabolites_importance_",phenotype,".csv"))
      write.csv(mic_importance,paste0("./results/Inter-Cor/global/O2PL1_microbes_importance_",phenotype,".csv"))
      
      ### bar plot 
      meta_plot_data <- meta_importance
      mic_plot_data <- mic_importance
      meta_name <- rownames(meta_plot_data)
      mic_name <- rownames(mic_plot_data)
      meta_df <- data.frame('name' = meta_name,'importance' = meta_plot_data$importance)
      mic_df <- data.frame('name' = mic_name,'importance' = mic_plot_data$importance)
      
      ### Top 10 metabolites or microbes
      meta_df$importance <- abs(meta_df$importance)
      meta_df <- meta_df[order(meta_df$importance,decreasing = F),]
      meta_df$importance <- as.numeric(meta_df$importance)
      
      if(nrow(meta_plot_data) > 10){
        meta_df <- meta_df[(nrow(meta_df)-9):nrow(meta_df),]
      }
      
      meta_df$name <- factor(meta_df$name,levels = meta_df$name)
      
      mic_df$importance <- abs(mic_df$importance)
      mic_df <- mic_df[order(mic_df$importance,decreasing = F),]
      mic_df$importance <- as.numeric(mic_df$importance)
      if(nrow(mic_plot_data) > 10){
        mic_df <- mic_df[(nrow(mic_df)-9):nrow(mic_df),]
      }
      mic_df$name <- factor(mic_df$name,levels = mic_df$name)
      
      ### plot
      
      meta_plot <- ggplot(data = meta_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#6300ea") + 
        labs(title = paste0("importance of metabolites in O2PLS1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(meta_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/O2PLS1_metabolites_importanceRank_",phenotype,".pdf"),width = 7,height = 5 )
      
      mic_plot <- ggplot(data = mic_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#00b692") + 
        labs(title = paste0("importance of microbes in O2PLS1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(mic_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/O2PLS1_microbes_importanceRank_",phenotype,".pdf"),width = 7,height = 5)
    }
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # o2pls for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        o2pls_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # o2pls for all samples
      o2pls_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      o2pls_func(metaData,micData,phenotype = "All")
    }
  }
  
  ### Mantel test
  else if (corMethod == "Mantel"){
    mt_func <- function(metadf,micdf,phenotype){
      veg.dist <- vegdist(metadf) # Bray-Curtis
      env.dist <- vegdist(micdf, "euclid")
      mantel_res <- mantel(veg.dist, env.dist, method="spear")
      
      res <- as.data.frame(cbind(mantel_res$statistic,mantel_res$signif)) 
      colnames(res) <- c("R","p")
      write.csv(res,paste0("./results/Inter-Cor/global/Mantel_test_",phenotype,".csv"),row.names = F)
    }
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # Mantel test for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        mt_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # Mantel test for all samples
      mt_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      mt_func(metaData,micData,phenotype = "All")
    }
  }
  
  ### sPLS
  else if (corMethod == "sPLS"){
    spls_func <- function(metadf,micdf,phenotype){
      spls_res <- spls(metadf,micdf) 
      
      # R and p
      meta_com1 <- spls_res[["variates"]][["X"]]
      mic_com1 <- spls_res[["variates"]][["Y"]]
      r <- cor.test(meta_com1,mic_com1,method = "pearson")
      
      meta_axis_c1_c2 <- spls_res[["variates"]][["X"]]
      mic_axis_c1_c2 <- spls_res[["variates"]][["Y"]]
      
      # arrow data 
      arrow_data <- cbind(data.frame(mic_axis_c1_c2),data.frame(meta_axis_c1_c2))
      colnames(arrow_data) <- c("mic_c1","mic_c2","meta_c1","meta_c2")
      
      v1 <- as.data.frame(spls_res[["variates"]][["X"]]) 
      v1$Type <- "Metabolites"
      v2 <- as.data.frame(spls_res[["variates"]][["Y"]])
      v2$Type <- "Microbes"
      
      # plot data
      plot_data <- rbind(v1,v2)
      
      cols=c('Metabolites' = "#6300ea", 'Microbes' = "#00b692")
      shapes=c('Metabolites' = 17,'Microbes' = 16)
      p <- ggplot(plot_data, aes(comp1, comp2)) +
        geom_segment(data = arrow_data,aes(x = mic_c1, y = mic_c2, xend = meta_c1, yend = meta_c2), arrow = arrow(length = unit(0.1, 'cm')),
                     color = '#696969', size = 1.2) +
        geom_point(data = v1, size = 4, aes(colour = "Metabolites", shape = "Metabolites") )+
        geom_point(data = v2, size = 4, aes(colour = "Microbes", shape = "Microbes")) +
        scale_colour_manual(name="Type",values = cols) +
        scale_shape_manual(name="Type",values = shapes) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
        labs(x = 'PC1', y = 'PC2', color = '') +
        geom_vline(xintercept = 0, color = 'black', linetype = 2, size = 0.3) +
        geom_hline(yintercept = 0, color = 'black', linetype = 2, size = 0.3) +
        ggtitle(paste0("r:",signif(r$estimate,3),"  p:",signif(r$p.value,3)))+
        theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),   
              plot.background = element_rect(fill = "transparent",colour = NA)) 
      ggsave(paste0("./results/Inter-Cor/global/sPLS_plot_",phenotype,".pdf"), p, width = 7, height = 5)
      
      ### importance of metabolites and microbes in Axis1
      meta_importance <- as.data.frame(spls_res[["loadings"]][["X"]][,1])
      colnames(meta_importance) <- "importance"
      rownames(meta_importance) <- colnames(metadf)
      mic_importance <- as.data.frame(spls_res[["loadings"]][["Y"]][,1])
      colnames(mic_importance) <- "importance"
      rownames(mic_importance) <- colnames(micdf)
      write.csv(meta_importance,paste0("./results/Inter-Cor/global/sPLS1_metabolites_importance_",phenotype,".csv"))
      write.csv(mic_importance,paste0("./results/Inter-Cor/global/sPLS1_microbes_importance_",phenotype,".csv"))
      
      ### bar plot 
      meta_plot_data <- meta_importance
      mic_plot_data <- mic_importance
      meta_name <- rownames(meta_plot_data)
      mic_name <- rownames(mic_plot_data)
      meta_df <- data.frame('name' = meta_name,'importance' = meta_plot_data$importance)
      mic_df <- data.frame('name' = mic_name,'importance' = mic_plot_data$importance)
      
      ### Top 10 metabolites or microbes
      meta_df$importance <- abs(meta_df$importance)
      meta_df <- meta_df[order(meta_df$importance,decreasing = F),]
      meta_df$importance <- as.numeric(meta_df$importance)
      
      if(nrow(meta_plot_data) > 10){
        meta_df <- meta_df[(nrow(meta_df)-9):nrow(meta_df),]
      }
      
      meta_df$name <- factor(meta_df$name,levels = meta_df$name)
      
      mic_df$importance <- abs(mic_df$importance)
      mic_df <- mic_df[order(mic_df$importance,decreasing = F),]
      mic_df$importance <- as.numeric(mic_df$importance)
      if(nrow(mic_plot_data) > 10){
        mic_df <- mic_df[(nrow(mic_df)-9):nrow(mic_df),]
      }
      mic_df$name <- factor(mic_df$name,levels = mic_df$name)
      
      ### plot
      
      meta_plot <- ggplot(data = meta_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#6300ea") + 
        labs(title = paste0("importance of metabolites in sPLS1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = 'black'),   
              panel.grid.minor = element_blank(),   
              panel.grid.major = element_blank(),  
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(meta_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/sPLS1_metabolites_importanceRank_",phenotype,".pdf"),width = 7,height = 5 )
      
      mic_plot <- ggplot(data = mic_df, aes(x = name, y = importance)) +
        geom_bar(stat = "identity", fill = "#00b692") + 
        labs(title = paste0("importance of microbes in sPLS1")  ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "transparent",colour = 'black'), 
              panel.grid.minor = element_blank(),  
              panel.grid.major = element_blank(),  
              plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_y_continuous(limits = c(0,as.integer(max(mic_df$importance) + 1)),expand = c(0,0)) +
        coord_flip() 
      ggsave(paste0("./results/Inter-Cor/global/sPLS1_microbes_importanceRank_",phenotype,".pdf"),width = 7,height = 5)
    }
    ## categorical phenotypes
    if(!is.na(phenoData) && phenoDataType == "categorical"){
      # sPLS test for each group
      for (i in unique(phenoData[,1])) {
        tempmetaData <- metaData[rownames(phenoData)[which(phenoData[,1] == i)],]
        tempmicData <- micData[rownames(phenoData)[which(phenoData[,1] == i)],]
        spls_func(metadf = tempmetaData,micdf = tempmicData,phenotype = i)
      }
      # sPLS test for all samples
      spls_func(metadf = metaData,micdf = micData,phenotype = "All")
    }
    ## continuous or no phenotypes
    else{
      spls_func(metaData,micData,phenotype = "All")
    }
  }
    
  else{
    stop("No such method!")
  }
}

