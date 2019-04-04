cnm <- abundance_matrix$X1
abundance_matrix2 <- abundance_matrix[,-1]
abundance_matrix3 <- t(abundance_matrix2)
colnames(abundance_matrix3) <- cnm
ProfileMatrix=abundance_matrix3   
Profile.NAME="GMProfile"      
ResultDir='~/'


k_cluster.row = 5
k_cluster.col = 25

hc_method = "ward.D2"

hc.row <- hcut(ProfileMatrix, k = k_cluster.row, stand = FALSE, 
                 hc_func = "hclust", hc_method = 'complete', hc_metric='pearson')
dd <- as.dendrogram(hc.row)
dd.reorder <- reorder(dd,c(30,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
plot(dd.reorder, main = "reorder(dd, 10:1)")
dd.reorder
#dist=dist(ProfileMatrix)
#?dist
#class = as.numeric(unname(hc.row$cluster))
#res=RearrangeJoseph(hc.row,as.matrix(dist),class,TRUE)
#hcl1=res$hcl
#?RearrangeJoseph
  ## Dendogram for col
#hc.col <- hcut(t(ProfileMatrix), k = k_cluster.col, stand = FALSE, 
                 hc_func = "hclust", hc_method = 'complete', hc_metric='spearman')
dist=dist(t(ProfileMatrix))
class = as.numeric( unname(hc.col$cluster) )
res2=RearrangeJoseph(hc.col,as.matrix(dist),class,TRUE)
hcl2=res2$hcl
  ## Two-way clustering using annHeatmap
  
k_cluster = 25
  
cols.map <- colorRampPalette(brewer.pal(10, "RdBu"))(256*ncols)
cols.hc  <- get_palette(palette = "simpsons", k=k_cluster*1.1)
  
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
  