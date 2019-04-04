set.seed(123)
x <- rnorm(10)
hc <- hclust(dist(x))

dd <- as.dendrogram(hc)
dd.reorder <- reorder(dd,c(1,10,40,10,10,40,10,20,10,1))
dd
op <- par(mfrow = c(1,3))
plot(dd, main = "random dendrogram 'dd'")
plot(dd.reorder, main = "reorder(dd, 10:1)")
#plot(reorder(dd, 10:1, agglo.FUN = mean), main = "reorder(dd, 10:1, mean)")
par(op)
?reorder
