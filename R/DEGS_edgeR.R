library( "DESeq2" )
library("RColorBrewer")
library(genefilter)
library(gplots)
library(pheatmap)
library(edgeR)
library(limma)
library(pcaExplorer)
library(rgl)
library(ggfortify)
library(cluster)
library('factoextra')
library('ggdendro')
library('magittr')



keep <- rowSums(cpm(count.matrix)>1) >= 12
table(keep)
counts <- count.matrix[keep, ]

labels <- c("glands", "glands", rep("epithelial",5), rep("fat",3), rep("mucosa",3),
            rep("epithelial",2),"glands",rep("muscle",2),rep("glands",3),
            rep("muscle",6),rep("mucosa",5),rep("immune",4))
D <- data.frame(row.names = colnames(counts))
D$Patients <- c(2,3,2,3,1,2,3,1,2,3,1,2,3,1,2,1,1,3,2,2,3,1,2,3,3,1,3,1,2,3,1,3,1,2,3,1)
D$Organ <- labels

Patients <- as.factor(D$Patients)
Organ <- as.factor(D$Organ)
design <- model.matrix(~Patients+Organ)
Organ # ref epithelial
Patients # ref 1

data <- DGEList(counts=counts, group=D$Organ)

data <- estimateDisp(data, design, robust=TRUE)
fit <- glmQLFit(data, design)
colnames(fit)

qlf <- glmQLFTest(fit, coef = 4:8)
topTags(qlf)


deg <- qlf$table
deg$id <- rownames(qlf$table)

deg.05 <- filter(deg, PValue <= 0.05)
head(deg.05)

counts.sig <- counts[deg.05$id, ]


matr.pca <- prcomp(na.omit(counts.sig), center = TRUE,scale. = TRUE, rank. = 15)
summary(matr.pca)
fviz_eig(matr.pca)
sdat <- t(scale(t(counts.sig)))

pr.dis <- dist(t(sdat), method = "euclidean")
hclust(pr.dis,method = "complete") -> xx
xx <-as.dendrogram(xx)
##colored by label
xx %>% set("labels_colors", k = 10) %>% plot(main = "Transcriptome")
## colored by branches
xx %>% set("branches_k_color", k = 10) %>% plot(main = "Transcriptome")
##boring
plot(xx)



