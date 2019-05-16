library('tidyverse')
library('factoextra')
library('ggdendro')
library('magittr')
## read in dmr file
matr=read_tsv('~/Downloads/Supplementary_Information/Ziller_et_al_DMR_finding/DMR_final_with_level.tsv')

## create pca matrix
matr.pca <- prcomp(na.omit(matr[,c(6:40)]), center = TRUE,scale. = TRUE, rank. = 15)
summary(matr.pca)
fviz_eig(matr.pca)
sdat <- t(scale(t(matr[,5:40])))

pr.dis <- dist(t(sdat), method = "euclidean")
hclust(pr.dis,method = "complete") -> xx
xx <-as.dendrogram(xx)
##colored by label
xx %>% set("labels_colors", k = 10) %>% plot(main = "Default colors")
## colored by branches
xx %>% set("branches_k_color", k = 10) %>% plot(main = "Default colors")
##boring
plot(xx)