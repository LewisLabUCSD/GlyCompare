library(ReorderCluster)


#######################################################################################
# Visualize significant correlations: Manhatann Plot for pathway vs phenotype/glycan  #  
# ----------------------------------------------------------------------------------- #
myA03.Clustering <- function(GMProfile=NULL, BioProfile=NULL, 
                             PathwayProfile=NULL, ResultDir=NULL){
  ## get Dendogram and Heatmap
  GP <- .myFun.doCluster(M=GMProfile,      Profile.NAME="GMProfile",      ResultDir)
  BP <- .myFun.doCluster(M=BioProfile,     Profile.NAME="BioProfile",     ResultDir)
  PP <- .myFun.doCluster(M=PathwayProfile, Profile.NAME="PathwayProfile", ResultDir)
  return(list(GP=GP,BP=BP,PP=PP))
}
  
# ----------------------------------------------------------------------------------- #
# Subfunctions for doing Clustering analysis
# ----------------------------------------------------------------------------------- #
## get Clustering
.myFun.doCluster <- function(M=NULL, Profile.NAME=NULL, ResultDir=NULL){
  require(Heatplus)
  require(RColorBrewer)
  require(ggpubr)
  require(factoextra)
  
  ## Parameter setting
  switch(Profile.NAME,
         GMProfile      = {cuth.r=13; cuth.c=12; ncols=0.9},
         BioProfile     = {cuth.r=12; cuth.c=12; ncols=1.5},
         PathwayProfile = {cuth.r=15; cuth.c=12; ncols=0.65}
  )
  
  ## Dendogram for row & col
  HC  <- .myFun.doHierCluster( ProfileMatrix=M,   hc_method = "ward.D2",
                               k_cluster.row = 3, k_cluster.col = 9)
  
  ## Two-way clustering using annHeatmap
  k_cluster=3
  cols.map <- colorRampPalette(brewer.pal(10, "RdBu"))(256*ncols)
  cols.hc  <- get_palette(palette = "aaas", k=k_cluster*1)
  
  labels <- list( Row = list(side = 4, cex=1.2, nrow = 20),
                  Col = list(labels = NULL))  
  
  if(!is.matrix(M)) {M <- as.matrix(M)}
  HP =  annHeatmap2( M, scale = "none",
                     cluster=list(Col=list(cuth=cuth.c), Row=list(cuth=cuth.r), col = cols.hc),
                     dendrogram = list(Col = list(dendro = as.dendrogram(HC$col)), 
                                       Row = list(dendro = as.dendrogram(HC$row))),                   
                     col=rev(cols.map), legend = 3, labels =labels)
  
  ## Visualization
  .myPlot.Heatmap(HC=HC,HP=HP, Profile.NAME=Profile.NAME, ResultDir=ResultDir)
  
  return(list(HP=HP, HC=HC))
}

## Dendogram for row & col
.myFun.doHierCluster <- function( ProfileMatrix = NULL,  hc_method = NULL, 
                                  k_cluster.row = NULL,  k_cluster.col = NULL){
  ## Dendogram for row
  hc.row <- hcut(ProfileMatrix, k = k_cluster.row, stand = FALSE, 
                 hc_func = "hclust", hc_method = hc_method)
  dist=dist(ProfileMatrix)
  class = as.numeric( unname(hc.row$cluster) )
  res=RearrangeJoseph(hc.row,as.matrix(dist),class,TRUE)
  hcl1=res$hcl

  ## Dendogram for col
  hc.col <- hcut(t(ProfileMatrix), k = k_cluster.col, stand = FALSE, 
                 hc_func = "hclust", hc_method = hc_method)
  dist=dist(t(ProfileMatrix))
  class = as.numeric( unname(hc.col$cluster) )
  res2=RearrangeJoseph(hc.col,as.matrix(dist),class,TRUE)
  hcl2=res2$hcl
  
  return(list(row=hcl1,col=hcl2))
}


## Visualization Heatmap with Clustering Dendogram 
.myPlot.Heatmap <- function(HC=NULL, HP=NULL, 
                            Profile.NAME=NULL, ResultDir=NULL){
  hc.row <- HC$row
  hc.col <- HC$col
  
  ## 1. Plot Dendogram
  pdf (file=paste0(ResultDir,Profile.NAME,".wLabel.pdf"), width=10, height=3)
  print(fviz_dend(hc.row, k=9,  rect = TRUE, cex = .5, labels_track_height = 5))
  print(fviz_dend(hc.col, k=9,  rect = TRUE, cex = .5, labels_track_height = 5))
  print(fviz_dend(hc.col, k=9,  rect = TRUE, cex = .5,  
                  show_labels = FALSE, rect_border = "grey", rect_lty =2,
                  labels_track_height = 5))
  dev.off()
  
  ## 2. Plot Dendogram w/o Label
  pdf (file=paste0(ResultDir,Profile.NAME,".woLabel.pdf"), width=5, height=3)
  print(fviz_dend(hc.row, k=9,  rect = TRUE, cex = .5,  
                  labels_track_height = 5))
  print(fviz_dend(hc.col, k=9,  rect = TRUE, cex = .5,  show_labels = FALSE,
                  rect_border = "grey", rect_lty =2,
                  labels_track_height = 5))
  dev.off()
  
  ## 3. Write Cluster
  if(Profile.NAME==GMProfile){
    fp = paste0(ResultDir,Profile.NAME,".Cluster2.csv")
    motif.cluster <- list()
    motif.cluster$ID <- hc.col$labels[ hc.col$order ]
    motif.cluster$cluster <- hc.col$cluster[ hc.col$order ]
    motif.cluster <- as.data.frame(motif.cluster)
    write.csv(motif.cluster, file=fp, row.names = FALSE)
  }
  
  ## 4. Plot Heatmap
  width=16; height=12
  pdf (paste0(ResultDir,Profile.NAME,".Heatmap.pdf"), width=16, height=12)
  plot(HP)
  dev.off()
}

