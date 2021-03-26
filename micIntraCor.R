# function - micIntraCor #
########################################################################
# File: micIntraCor.R
# Aim : Intra Correlation analysis for microbes
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
##    micIndData --- A dataframe of microbes                                        ##
##                       (rows:samples,columns:microbes)                            ##
##    corMethod --- The method of correlation analysis                              ##
##                Default: spearman                                                 ##
##    count --- SparCC flag                                                         ##
##                Default: TRUE                                                     ##
##    confounderData --- A dataframe of confounder                                  ##
##                       (rows:samples,columns:confounder)                          ##
##                Default: NA                                                       ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    mic_intra_result --- A list included coefficient, p value and p.adjust value  ##
##                          of intra correlation analysis results for microbes      ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
micIntraCor <- function(micIndData,corMethod = "spearman",count = TRUE,confounderData = NA)
{
  
  if(!file.exists("./results/Intra-Cor/microbiome")){
    dir.create("./results/Intra-Cor/microbiome",recursive = T)
  }
  mic_intra_result <-list()
  
  ### Method 1 - Generalized coRrelation analysis for micbolome and Microbiome (GRaMM)
  if(corMethod == "gramm"){
    ### microbes to microbes
    naivegramm_res <- naivegramm(micIndData,micIndData,covdata = confounderData)
    write.csv(naivegramm_res[["r_All"]],"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
    write.csv(naivegramm_res[["p_All"]], "./results/Intra-Cor/microbiome/mic_intra_p.csv")
    naivegramm_res[["p.adjust_All"]] = apply(naivegramm_res[["p_All"]], 2, FUN = function(x)
      p.adjust(x,method = "BH"))
    naivegramm_res[["p.adjust_All"]] = as.data.frame(naivegramm_res[["p.adjust_All"]]) 
    write.csv(naivegramm_res[["p.adjust_All"]], "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
    write.csv(naivegramm_res[["type_All"]], "./results/Intra-Cor/microbiome/mic_intra_cor_type.csv")
    mic_intra_cor = naivegramm_res[["r_All"]]
    mic_intra_p = naivegramm_res[["p_All"]]
    mic_intra_p.adjust = naivegramm_res[["p.adjust_All"]]
    mic_intra_result <- list(r = as.data.frame(mic_intra_cor), p = as.data.frame(mic_intra_p), p.adjust = as.data.frame(mic_intra_p.adjust))
  }
  
  ### Method 2 - Spearman,Pearson,Kendall
  if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
    ### microbes to microbes
    mic_intra_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (micIndData))
    dim(mic_intra_cor)
    rownames (mic_intra_cor) = colnames (micIndData)
    colnames (mic_intra_cor) = colnames (micIndData)
    mic_intra_p = mic_intra_cor
    mic_intra_p.adjust = mic_intra_cor
    for (m in colnames (micIndData)) {
      mic_intra_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
        cor.test (x, micIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
      mic_intra_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
        cor.test (x, micIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
    }
    mic_intra_p.adjust = apply(mic_intra_p, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    write.csv(mic_intra_cor,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
    write.csv(mic_intra_p, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
    write.csv(mic_intra_p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
    mic_intra_result <- list(r = as.data.frame(mic_intra_cor), p = as.data.frame(mic_intra_p), p.adjust = as.data.frame(mic_intra_p.adjust))
    
  }
  
  ### Method 3 - Partial correlation
  if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
    corMethod = strsplit(corMethod," ")[[1]][2]
    ### microbes to microbes
    mic_intra_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (micIndData))
    dim(mic_intra_cor)
    rownames (mic_intra_cor) = colnames (micIndData)
    colnames (mic_intra_cor) = colnames (micIndData)
    mic_intra_p = mic_intra_cor
    mic_intra_p.adjust = mic_intra_cor
    for (m in colnames (micIndData)) {
      # exist confounder
      if(is.na(confounderData)){
        stop("Lack of confounders!")
      }else{
        mic_intra_cor [ , m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          pcor.test (x, micIndData[,m], confounderData,method = corMethod)$estimate)
        mic_intra_p[,m] = apply (micIndData, MARGIN = 2, FUN = function (x) 
          pcor.test (x, micIndData[,m], confounderData,method = corMethod)$p.value)
      }
    }
    mic_intra_p.adjust = apply(mic_intra_p, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    mic_intra_p[which(is.na(mic_intra_p))]<-0
    mic_intra_p.adjust[which(is.na(mic_intra_p.adjust))]<-0
    write.csv(mic_intra_cor,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
    write.csv(mic_intra_p, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
    write.csv(mic_intra_p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
    mic_intra_result <- list(r = as.data.frame(mic_intra_cor), p = as.data.frame(mic_intra_p), p.adjust = as.data.frame(mic_intra_p.adjust))
  }
  
  
  ### Method 4 - SparCC
  if(corMethod == "SparCC"){
    if(count == T){
      mic_intra_result <- SparCC.both(x = micIndData)
      mic_intra_result$cor <- as.data.frame(mic_intra_result$cor)
      colnames(mic_intra_result$cor) <- colnames(micIndData)
      rownames(mic_intra_result$cor) <- colnames(micIndData)
      mic_intra_result$p <- as.data.frame(mic_intra_result$p)
      colnames(mic_intra_result$p) <- colnames(micIndData)
      rownames(mic_intra_result$p) <- colnames(micIndData)
      mic_intra_result$p.adjust = as.data.frame(apply(mic_intra_result$p, 2, FUN = function(x){
        p.adjust(x,method = "BH")})) 
      colnames(mic_intra_result$p.adjust) <- colnames(micIndData)
      rownames(mic_intra_result$p.adjust) <- colnames(micIndData)
      names(mic_intra_result) <- c("r","p","p.adjust")
      write.csv(mic_intra_result$cor,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
      write.csv(mic_intra_result$p, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
      write.csv(mic_intra_result$p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
      mic_intra_cor <- mic_intra_result$cor
      mic_intra_p <- mic_intra_result$p
      mic_intra_p.adjust <- mic_intra_result$p.adjust
      mic_intra_result <- list(r = as.data.frame(mic_intra_result$cor), p = as.data.frame(mic_intra_result$p), p.adjust = as.data.frame(mic_intra_result$p.adjust))
    }
    else{
      mic_intra_result <- SparCC.frac(x = micIndData)
      mic_intra_result$cor <- as.data.frame(mic_intra_result$cor)
      colnames(mic_intra_result$cor) <- colnames(micIndData)
      rownames(mic_intra_result$cor) <- colnames(micIndData)
      mic_intra_result$p <- as.data.frame(mic_intra_result$p)
      colnames(mic_intra_result$p) <- colnames(micIndData)
      rownames(mic_intra_result$p) <- colnames(micIndData)
      mic_intra_result$p.adjust = as.data.frame(apply(mic_intra_result$p, 2, FUN = function(x){
        p.adjust(x,method = "BH")})) 
      colnames(mic_intra_result$p.adjust) <- colnames(micIndData)
      rownames(mic_intra_result$p.adjust) <- colnames(micIndData)
      names(mic_intra_result) <- c("r","p","p.adjust")
      write.csv(mic_intra_result$cor,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
      write.csv(mic_intra_result$p, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
      write.csv(mic_intra_result$p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
      mic_intra_cor <- mic_intra_result$cor
      mic_intra_p <- mic_intra_result$p
      mic_intra_p.adjust <- mic_intra_result$p.adjust
      mic_intra_result <- list(r = as.data.frame(mic_intra_result$cor), p = as.data.frame(mic_intra_result$p), p.adjust = as.data.frame(mic_intra_result$p.adjust))
    }
  }
  
  ### Method 5 - CCLasso
  if(corMethod == "CCLasso"){
    micIndData <- micIndData + 0.001
    mic_intra_result <- cclasso(x = micIndData)
    mic_intra_result$cor_w <- as.data.frame(mic_intra_result$cor_w)
    colnames(mic_intra_result$cor_w) <- colnames(micIndData)
    rownames(mic_intra_result$cor_w) <- colnames(micIndData)
    mic_intra_result$p_vals <- as.data.frame(mic_intra_result$p_vals)
    colnames(mic_intra_result$p_vals) <- colnames(micIndData)
    rownames(mic_intra_result$p_vals) <- colnames(micIndData)
    mic_intra_result$p.adjust = as.data.frame(apply(mic_intra_result$p_vals, 2, FUN = function(x){
      p.adjust(x,method = "BH")})) 
    colnames(mic_intra_result$p.adjust) <- colnames(micIndData)
    rownames(mic_intra_result$p.adjust) <- colnames(micIndData)
    write.csv(mic_intra_result$cor_w,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
    write.csv(mic_intra_result$p_vals, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
    write.csv(mic_intra_result$p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
    mic_intra_cor <- mic_intra_result$cor_w
    mic_intra_p <- mic_intra_result$p_vals
    mic_intra_p.adjust <- mic_intra_result$p.adjust
    mic_intra_result <- list(r = as.data.frame(mic_intra_result$cor_w), p = as.data.frame(mic_intra_result$p_vals), p.adjust = as.data.frame(mic_intra_result$p.adjust))
  }
  
  ### Method 6 - GLM
  if(corMethod == "glm"){
    ### microbes to microbes
    mic_intra_cor = matrix (NA, nrow = ncol (micIndData), ncol = ncol (micIndData))
    dim(mic_intra_cor)
    rownames (mic_intra_cor) = colnames (micIndData)
    colnames (mic_intra_cor) = colnames (micIndData)
    mic_intra_p = mic_intra_cor
    mic_intra_p.adjust = mic_intra_cor
    for (m in colnames (micIndData)) {
      # exist confounder
      if(is.na(confounderData)){
        mic_intra_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(micIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
        mic_intra_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(micIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
      }else{
        temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
        mic_intra_cor [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x)
          summary(glm(as.formula(paste0("micIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
        mic_intra_p [ , m] <- apply (micIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(as.formula(paste0("micIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
      }
    }
    mic_intra_p.adjust = apply(mic_intra_p, 2, FUN = function(x)
      p.adjust(x,method = "BH"))
    write.csv(mic_intra_cor,"./results/Intra-Cor/microbiome/mic_intra_cor.csv")
    write.csv(mic_intra_p, "./results/Intra-Cor/microbiome/mic_intra_p.csv")
    write.csv(mic_intra_p.adjust, "./results/Intra-Cor/microbiome/mic_intra_p.adjust.csv")
    mic_intra_result <- list(r = as.data.frame(mic_intra_cor), p = as.data.frame(mic_intra_p), p.adjust = as.data.frame(mic_intra_p.adjust))
  }
  ### Method 7 - MIC
  if(corMethod == "MIC"){
    ### microbes to microbes
    mic_intra_cor <- mic_intra_p.adjust <- matrix(0, nrow = ncol(micIndData),ncol = ncol(micIndData))
    rownames(mic_intra_p.adjust) <- rownames(mic_intra_cor) <- colnames(micIndData)
    colnames(mic_intra_p.adjust) <- colnames(mic_intra_cor) <- colnames(micIndData)
    for(i in seq_len(ncol(micIndData))){
      for(j in seq_len(ncol(micIndData))){
        x1 <- micIndData[,i]
        y1 <- micIndData[,j]
        ### MIC
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
        prsmic <- cor.test(y1,x1,method="pearson")
        prsr <- prsmic$estimate
        ltmic <- 1 - micr + prsr^2
        mic_intra_p.adjust[i,j] <- micp
        mic_intra_cor[i,j] <- micr
      }
    }
    p.adjust_result = apply(mic_intra_p.adjust, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    mic_intra_result <- list(r = as.data.frame(mic_intra_cor), p = as.data.frame(mic_intra_p.adjust), p.adjust = as.data.frame(p.adjust_result))
  }
  ######### heatmap #########
  mic_intra_cor <- as.data.frame(mic_intra_cor)
  mic_intra_p.adjust <- as.data.frame(mic_intra_p.adjust) 
  
  ### create matrix for heatmap
  plotmat = mic_intra_cor

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
    stop("Unable to generate heat map due to too few significant correlation pairs")
  }else{
    row_dendrogram = FALSE
    col_dendrogram = TRUE
  }
  
  ### create matrix with significance stars
  plotmat_p = mic_intra_p.adjust ### p-values to make stars for heatmap
  # rownames (plotmat_p) <- nmic
  # colnames (plotmat_p) <- nmic
  stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
  rownames (stars) <- rownames(mic_intra_p.adjust)
  colnames (stars) <- colnames(mic_intra_p.adjust)
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
           fontsize_number = 10,
           number_color = "black",
           cluster_row = row_dendrogram,
           cluster_cols = col_dendrogram,
           # annotation_row = annotation_row,
           # annotation_colors = ann_colors,
           angle_col = "90",filename = "./results/Intra-Cor/microbiome/heatmap3M.pdf")
  
  
  ######### chord #########
  
  CairoPDF('./results/Intra-Cor/microbiome/Chord.pdf',width = 8,height = 8)
  migration <- as.matrix(mic_intra_cor) 
  # sectorcol <- rainbow(nrow(migration) * ncol(migration))
  linkcol <- c()
  for (i in 1:ncol(mic_intra_cor)) {
    for(j in 1:nrow(mic_intra_cor)){
      if(mic_intra_cor[j,i] > 0) {
        linkcol <- append(linkcol,"#d7191c")
      }
      else{
        linkcol <- append(linkcol,"#4bacc6")
      }
    }
  }
  linktran <- c()
  for (i in 1:ncol(mic_intra_p.adjust)) {
    for(j in 1:nrow(mic_intra_p.adjust)){
      temp_tran <- mic_intra_p.adjust[j,i]
      linktran <- append(linktran,temp_tran)
    }
  }
  
  circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
  if(ncol(mic_intra_cor) <= 2){
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
  
  
  ######### network #########
  ### Generate edge table
  CorrDF <- function(cormat, pmat) {
    ut = matrix (TRUE, nrow = nrow (cormat), ncol = ncol (pmat))
    data.frame(
      from = rownames(cormat)[row(cormat)[ut]],
      to = colnames(cormat)[col(cormat)[ut]],
      r =(cormat)[ut],
      p.adjust = pmat[ut]
    )
  }
  edge_df <- CorrDF(mic_intra_cor,mic_intra_p.adjust)
  if(nrow(edge_df) == 0 ){
    stop("There are no edges that satisfy these thresholds!")
  }
  nodeattrib <- data.frame(node = union(edge_df$from,edge_df$to))
  co_net <- graph_from_data_frame(edge_df,direct = F, vertices = nodeattrib) 
  V(co_net)$shape <- "circle"
  V(co_net)$color <- "#00b692"
  
  ### Set edges color -- Positive:red; Negative: blue
  edge_col <- c()
  for (i in edge_df$r)
  {
    if (i > 0 ){
      edge_col <- append(edge_col,"#d7191c")
    }else{
      edge_col <- append(edge_col,"#4bacc6")
    }
  }
  co_net <- co_net %>% set_edge_attr("color", value = edge_col)
  
  ### Set edges size according to the value of r
  edge_size <- edge_df$r
  co_net <- co_net %>% set_edge_attr("weight", value = edge_size)
  
  ### Set nodes size according to the value of degre
  nodes_size = log2(centr_degree(co_net)$res + 2) * 3
  
  # ### The nodes in the network are clustered using fast greedy
  # cfg <- cluster_fast_greedy(as.undirected(co_net))
  # cluster_modules <- sort(table(membership(cfg)),decr=T)
  
  # ### modules mark groups is opposite to mods_list_cs
  # mark_groups <- list()
  # 
  # for (i in 1:length(cluster_modules)) {
  #   mark_groups[[i]] <- mods_list_cs[[length(cluster_modules) + 1 - i]]
  # }
  
  ### Set the number of random seeds, and start from the same number of 
  ### random seeds in subsequent drawing to ensure that the shapes of the 
  ### drawing before and after are corresponding
  pdf(paste0("./results/Intra-Cor/microbiome/network3M_label.pdf"),width = 7,height = 5)
  set.seed(123)
  coords <- layout_(co_net,with_fr(niter = 9999, grid = "auto"))
  plot(co_net,vertex.label = V(co_net)$name,vertex.label.cex = edge_size,vertex.label.color = "black",
       layout = coords,vertex.size = nodes_size,
       edge.lty = 1,edge.curved=TRUE,
       margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight)
  )
  dev.off()
  
  pdf(paste0("./results/Intra-Cor/microbiome/network3M.pdf"),width = 7,height = 5)
  set.seed(123)
  coords <- layout_(co_net,with_fr(niter = 9999, grid = "auto"))
  plot(co_net,layout = coords,vertex.label = NA,
       vertex.size = nodes_size,
       edge.lty = 1,edge.curved=TRUE,
       margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight)
  )
  dev.off()
  
  return(mic_intra_result)
}
