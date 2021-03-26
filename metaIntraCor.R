# function - metaIntraCor #
########################################################################
# File: metaIntraCor.R
# Aim : Intra Correlation analysis for metabolites
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
##    metaIndData --- A dataframe of metabolites                                    ##
##                       (rows:samples,columns:metabolites)                         ##
##    corMethod --- The method of correlation analysis                              ##
##                Default: spearman                                                 ##
##    confounderData --- A dataframe of confounder                                  ##
##                       (rows:samples,columns:confounder)                          ##
##                Default: NA                                                       ##
##                                                                                  ##
######################################################################################
## Output:                                                                          ## 
##    meta_intra_result --- A list included coefficient, p value and p.adjust value ##
##                          of intra correlation analysis results for metabolites   ##
##                                                                                  ##
######################################################################################
#---------------------------------------------------------------------------------------------------------------------
metaIntraCor <- function(metaIndData,corMethod = "spearman",confounderData = NA)
{
  if(!file.exists("./results/Intra-Cor/metabolome")){
    dir.create("./results/Intra-Cor/metabolome",recursive = T)
  }
  ### Method 1 - Generalized coRrelation analysis for Metabolome and Microbiome (GRaMM)
  if(corMethod == "gramm"){
    ### metabolites to metabolites
    naivegramm_res <- naivegramm(metaIndData,metaIndData,covdata = confounderData)
    write.csv(naivegramm_res[["r_All"]],"./results/Intra-Cor/metabolome/mic_meta_cor.csv")
    write.csv(naivegramm_res[["p_All"]], "./results/Intra-Cor/metabolome/mic_meta_p.csv")
    naivegramm_res[["p.adjust_All"]] = apply(naivegramm_res[["p_All"]], 2, FUN = function(x)
      p.adjust(x,method = "BH"))
    naivegramm_res[["p.adjust_All"]] = as.data.frame(naivegramm_res[["p.adjust_All"]]) 
    write.csv(naivegramm_res[["p.adjust_All"]], "./results/Intra-Cor/metabolome/mic_meta_p.adjust.csv")
    write.csv(naivegramm_res[["type_All"]], "./results/Intra-Cor/metabolome/mic_meta_cor_type.csv")
    meta_intra_cor = naivegramm_res[["r_All"]]
    meta_intra_p = naivegramm_res[["p_All"]]
    meta_intra_p.adjust = naivegramm_res[["p.adjust_All"]]
    
    
    meta_intra_result <- list(r = as.data.frame(meta_intra_cor), p = as.data.frame(meta_intra_p), p.adjust = as.data.frame(meta_intra_p.adjust))
  }
  
  ### Method 2 - Spearman,Pearson,Kendall
  if(corMethod == "spearman" || corMethod == "pearson" || corMethod == "kendall"){
    ### metabolites to metabolites
    meta_intra_cor = matrix (NA, nrow = ncol (metaIndData), ncol = ncol (metaIndData))
    dim(meta_intra_cor)
    rownames (meta_intra_cor) = colnames (metaIndData)
    colnames (meta_intra_cor) = colnames (metaIndData)
    meta_intra_p = meta_intra_cor
    meta_intra_p.adjust = meta_intra_cor
    for (m in colnames (metaIndData)) {
      meta_intra_cor [ , m] = apply (metaIndData, MARGIN = 2, FUN = function (x) 
        cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$estimate)
      meta_intra_p[,m] = apply (metaIndData, MARGIN = 2, FUN = function (x) 
        cor.test (x, metaIndData[,m], method = corMethod, use = "pairwise.complete.obs")$p.value)
    }
    meta_intra_p.adjust = apply(meta_intra_p, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    write.csv(meta_intra_cor,"./results/Intra-Cor/metabolome/meta_intra_cor.csv")
    write.csv(meta_intra_p, "./results/Intra-Cor/metabolome/meta_intra_p.csv")
    write.csv(meta_intra_p.adjust, "./results/Intra-Cor/metabolome/meta_intra_p.adjust.csv")
    meta_intra_result <- list(r = as.data.frame(meta_intra_cor), p = as.data.frame(meta_intra_p), p.adjust = as.data.frame(meta_intra_p.adjust))
  }
  
  ### Method 3 - Partial correlation
  if(corMethod == "partial spearman" || corMethod == "partial pearson" || corMethod == "partial kendall"){
    corMethod = strsplit(corMethod," ")[[1]][2]
    ### metabolites to metabolites
    meta_intra_cor = matrix (NA, nrow = ncol (metaIndData), ncol = ncol (metaIndData))
    dim(meta_intra_cor)
    rownames (meta_intra_cor) = colnames (metaIndData)
    colnames (meta_intra_cor) = colnames (metaIndData)
    meta_intra_p = meta_intra_cor
    meta_intra_p.adjust = meta_intra_cor
    for (m in colnames (metaIndData)) {
      # exist confounder
      if(is.na(confounderData)){
        stop("Lack of confounders!")
      }else{
        meta_intra_cor [ , m] = apply (metaIndData, MARGIN = 2, FUN = function (x) 
         { pcor.test (x, metaIndData[,m], confounderData,method = corMethod)$estimate})
        meta_intra_p[,m] = apply (metaIndData, MARGIN = 2, FUN = function (x) 
          {pcor.test (x, metaIndData[,m], confounderData,method = corMethod)$p.value})
      }
    }
    meta_intra_p.adjust = apply(meta_intra_p, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    meta_intra_p[which(is.na(meta_intra_p))]<-0
    meta_intra_p.adjust[which(is.na(meta_intra_p.adjust))]<-0
    write.csv(meta_intra_cor,"./results/Intra-Cor/metabolome/meta_intra_cor.csv")
    write.csv(meta_intra_p, "./results/Intra-Cor/metabolome/meta_intra_p.csv")
    write.csv(meta_intra_p.adjust, "./results/Intra-Cor/metabolome/meta_intra_p.adjust.csv")
    meta_intra_result <- list(r = as.data.frame(meta_intra_cor), p = as.data.frame(meta_intra_p), p.adjust = as.data.frame(meta_intra_p.adjust))
  }
  
  ### Method 4 - GLM
  if(corMethod == "glm"){
    ### metabolites to metabolites
    meta_intra_cor = matrix (NA, nrow = ncol (metaIndData), ncol = ncol (metaIndData))
    dim(meta_intra_cor)
    rownames (meta_intra_cor) = colnames (metaIndData)
    colnames (meta_intra_cor) = colnames (metaIndData)
    meta_intra_p = meta_intra_cor
    meta_intra_p.adjust = meta_intra_cor
    for (m in colnames (metaIndData)) {
      # exist confounder
      if(is.na(confounderData)){
        meta_intra_cor [ , m] <- apply (metaIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,1])
        meta_intra_p [ , m] <- apply (metaIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(metaIndData[,m] ~ x,family = gaussian))$coefficients[2,4])
      }else{
        temp_formula <- paste0("confounderData$",colnames(confounderData),collapse = "+")
        meta_intra_cor [ , m] <- apply (metaIndData, MARGIN = 2, FUN = function (x)
          summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,1])
        meta_intra_p [ , m] <- apply (metaIndData, MARGIN = 2, FUN = function (x) 
          summary(glm(as.formula(paste0("metaIndData[,m] ~ x + ",temp_formula)),family = gaussian))$coefficients[2,4])
      }
    }
    meta_intra_p.adjust = apply(meta_intra_p, 2, FUN = function(x)
      p.adjust(x,method = "BH"))
    write.csv(meta_intra_cor,"./results/Intra-Cor/metabolome/mic_meta_cor.csv")
    write.csv(meta_intra_p, "./results/Intra-Cor/metabolome/mic_meta_p.csv")
    write.csv(meta_intra_p.adjust, "./results/Intra-Cor/metabolome/mic_meta_p.adjust.csv")
    meta_intra_result <- list(r = as.data.frame(meta_intra_cor), p = as.data.frame(meta_intra_p), p.adjust = as.data.frame(meta_intra_p.adjust))
  }
  
  ## Method 5 - MIC
  if(corMethod == "MIC"){
    ### metabolites to metabolites
    meta_intra_cor <- meta_intra_p.adjust <- matrix(0, nrow = ncol(metaIndData),ncol = ncol(metaIndData))
    rownames(meta_intra_p.adjust) <- rownames(meta_intra_cor) <- colnames(metaIndData)
    colnames(meta_intra_p.adjust) <- colnames(meta_intra_cor) <- colnames(metaIndData)
    for(i in seq_len(ncol(metaIndData))){
      for(j in seq_len(ncol(metaIndData))){
        x1 <- metaIndData[,i]
        y1 <- metaIndData[,j]
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
        meta_intra_p.adjust[i,j] <- micp
        meta_intra_cor[i,j] <- micr
      }
    }
    p.adjust_result = apply(meta_intra_p.adjust, 2, FUN = function(x){
      p.adjust(x,method = "BH")})
    meta_intra_result <- list(r = as.data.frame(meta_intra_cor), p = as.data.frame(meta_intra_p.adjust), p.adjust = as.data.frame(p.adjust_result))
  }
  ######### heatmap #########
  meta_intra_cor <- as.data.frame(meta_intra_cor)
  meta_intra_p.adjust <- as.data.frame(meta_intra_p.adjust) 
  
  ### create matrix for heatmap
  plotmat = meta_intra_cor
  
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
  plotmat_p = meta_intra_p.adjust ### p-values to make stars for heatmap

  stars = matrix ("", ncol = ncol (plotmat_p), nrow = nrow (plotmat_p))
  rownames (stars) <- rownames(meta_intra_p.adjust)
  colnames (stars) <- colnames(meta_intra_p.adjust)
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
           angle_col = "90",filename = "./results/Intra-Cor/metabolome/heatmap3M.pdf")
  
  
  ######### chord #########
  
  CairoPDF('./results/Intra-Cor/metabolome/Chord.pdf',width = 8,height = 8)
  migration <- as.matrix(meta_intra_cor) 
  # sectorcol <- rainbow(nrow(migration) * ncol(migration))
  linkcol <- c()
  for (i in 1:ncol(meta_intra_cor)) {
    for(j in 1:nrow(meta_intra_cor)){
      if(meta_intra_cor[j,i] > 0) {
        linkcol <- append(linkcol,"#d7191c")
      }
      else{
        linkcol <- append(linkcol,"#4bacc6")
      }
    }
  }
  linktran <- c()
  for (i in 1:ncol(meta_intra_p.adjust)) {
    for(j in 1:nrow(meta_intra_p.adjust)){
      temp_tran <- meta_intra_p.adjust[j,i]
      linktran <- append(linktran,temp_tran)
    }
  }
  
  circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
  if(ncol(meta_intra_cor) <= 2){
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
  edge_df <- CorrDF(meta_intra_cor,meta_intra_p.adjust)
  if(nrow(edge_df) == 0 ){
    stop("There are no edges that satisfy these thresholds!")
  }
  nodeattrib <- data.frame(node = union(edge_df$from,edge_df$to))
  co_net <- graph_from_data_frame(edge_df,direct = F, vertices = nodeattrib) 
  V(co_net)$shape <- "circle"
  V(co_net)$color <- "#6300ea"
  
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
  pdf(paste0("./results/Intra-Cor/metabolome/network3M_label.pdf"),width = 7,height = 5)
  set.seed(123)
  coords <- layout_(co_net,with_fr(niter = 9999, grid = "auto"))
  plot(co_net,vertex.label = V(co_net)$name,vertex.label.cex = edge_size,vertex.label.color = "black",
       layout = coords,vertex.size = nodes_size,
       edge.lty = 1,edge.curved=TRUE,
       margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight)
  )
  dev.off()
  
  pdf(paste0("./results/Intra-Cor/metabolome/network3M.pdf"),width = 7,height = 5)
  set.seed(123)
  coords <- layout_(co_net,with_fr(niter = 9999, grid = "auto"))
  plot(co_net,layout = coords,vertex.label = NA,
       vertex.size = nodes_size,
       edge.lty = 1,edge.curved=TRUE,
       margin = c(0,0,0,0),edge.width = abs(E(co_net)$weight)
  )
  dev.off()
  
  return(meta_intra_result)
}
