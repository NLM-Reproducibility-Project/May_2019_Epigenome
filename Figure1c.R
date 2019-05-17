# figure 1c
library(dplyr)

data <- read.csv("~/Documents/NLM_Repro/NLM_Repro_Workshop/data/DMR_final_with_level.tsv", header = T, sep = "\t")
dmr <- data[ , -c(1:4)] # dim: 1,198,131
dmr <- na.omit(dmr) # rm NAs, dim: 1,143,033

# # filter for >= 0.3 methylation level
# dmr.30 <- filter_all(dmr, all_vars(. >= 0.3))

dmr.dist <- dist(t(dmr))

dmr.pca <- cmdscale(dmr.dist, k=15)
r <- rownames(dmr.pca)
r <- sub("methylation_level_", "", r)
rownames(dmr.pca) <- r

library(gplots)
dmr.pca <- t(dmr.pca)
labels <- c("glands", "glands", rep("epithelial",5), rep("fat",3), rep("mucosa",3),
            rep("epithelial",2),"glands",rep("muscle",2),rep("glands",3),
            rep("muscle",6),rep("mucosa",5),rep("immune",4))

heatmap.2(dmr.pca[, 1:15], hclustfun=function(x) hclust(x, method="ward.D2"), 
          dendrogram = "column")
