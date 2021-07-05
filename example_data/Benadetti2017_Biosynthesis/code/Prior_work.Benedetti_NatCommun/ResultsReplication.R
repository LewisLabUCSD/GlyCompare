#### Initialize ----

library(GeneNet)
library(readxl)
library(ggplot2)

setwd("./Benedetti_NatCommun/")

#### Define Functions ----

# build contingency table from two adjacency matrices
contab <- function(adja, data) {
  ct_a <- length(which(adja == 1 & data == 1)) # true positive, upper left
  ct_b <- length(which(adja == 0 & data == 1)) # false positive, upper right
  ct_c <- length(which(adja == 1 & data == 0)) # false negative, lower left
  ct_d <- length(which(adja == 0 & data == 0)) # true negative, lower right
  ct <- cbind(c(ct_a, ct_c), c(ct_b, ct_d))
  return(ct)
}

#### Load Data ----

# load preprocessed glycomics data
df <- as.data.frame(read_excel(path="SupplementalMaterial_DatasetS1_PreprocessedData.xls",sheet="Korčula2013_residuals",col_names=T))
data <- df[,2:dim(df)[2]]
rownames(data) <- df[,1]

# load prior knowledge adjacency matrix (inferred pathway in our paper)
adja <- read.csv("InferredPathway.csv", header = TRUE, sep = ";", dec = ",", row.names = 1)
adja[is.na(adja)] <- 0

#### Compute Partial Correlations ----

# compute partial correlation using GeneNet
pcor <- ggm.estimate.pcor(as.matrix(data), method = "dynamic")
# compute associated p-values
pvals <- network.test.edges(pcor)
# apply multiple testing correction using Benjamini-Hochberg method for FDR
pvals$padj <- p.adjust(pvals$pval,"BH")
# get correlation cutoff for GGM
pvals$pval[pvals$pval > 0.01] <- NaN
fdr_threshold <- abs(pvals$pcor[which.max(pvals$pval)])

#### Plot Partial Correlations at Pathway Distance 1 ----

pd1 <- pcor[adja==1]
pd1 <- pd1[!is.na(pd1)]

p <- data.frame(pcor=pd1) %>% 
  ggplot(aes(x="",y=pcor)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylim(-1,1) +
  geom_hline(yintercept = fdr_threshold, color = "gray50", linetype = 2) +
  geom_hline(yintercept = - fdr_threshold, color = "gray50", linetype = 2) +
  theme_bw() +
  xlab("Pathway Distance 1") +
  ylab("Partial Correlation Coefficients") +
  ggtitle("Distribution of Coefficients at Pathway Distance 1")
p

#### Compute Fisher's Exact Test p-value ----

# create adjacency matrix from GGM
adja_data <- pcor
adja_data[abs(adja_data) >= fdr_threshold] <- 1
adja_data[abs(adja_data) < fdr_threshold] <- 0

# only consider upper triangular matrix (since both adjacency and correlation matrices are symmetric)
adja_data[lower.tri(adja_data, diag = TRUE)] <- NA
adja[lower.tri(adja, diag = TRUE)] <- NA

# compute contingency table
contin <- contab(adja, adja_data)
# compute Fisher's exact test p-value (cfr. Figure 5B for the selected model)
fis <- fisher.test(contin)
fis$p.value