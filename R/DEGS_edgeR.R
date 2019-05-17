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


# import transcriptomic matrix with FPKM
raw.fpkm <- read.csv("data/RNA_matrix_FPKM.tsv", header = T, sep = '\t')

# data cleaning
fpkm <- raw.fpkm[ , 10:45]
tmp <- colnames(fpkm)
tmp <- sub("_FPKM", "", tmp)rownames(fpkm) <- raw.fpkm$gene_id
colnames(fpkm) <- tmp
rm(tmp)
dim(fpkm)

# data file with chromosome coordinates
locus <- read.csv("data/location.csv", header = F, stringsAsFactors = F)
locus <- locus$V1
transcript.length <- c()                        #toal number of transcripts

# for loop end-start to get the gene length
for (l in locus) {
  a <- strsplit(l, "-")
  b <- unlist(a)
  len <- abs(as.integer(b[2]) - as.integer(b[1]))
  len <- len / 1000
  transcript.length <- c(transcript.length, len)
}


# import data with total read count mapped fragments
total.mapped <- read.csv("data/total.mapped.csv", header = T)
total.mapped <- total.mapped$Properly.Paired.RNA.seq.Reads


# function to convert FPKM to read counts
convertFPKM.to.counts <- function (frag, tran, total) {
  count.matrix <- matrix(nrow = length(tran))
  for (c in 1:ncol(frag)){
    D <- c()
    #print(D)
    for (r in 1:nrow(frag)){
      #print(frag[c,r])
      counts <- (frag[r,c] * tran[r] * total[c]) / 1000000
      counts <- as.integer(counts)
      #print(counts)
      D <- c(D, counts)
    }
    #print (length(D))
    count.matrix <- cbind(count.matrix, D)
  }
  return(count.matrix)
}

# count.matrix
count.matrix <- convertFPKM.to.counts(fpkm, transcript.length, total.mapped)
count.matrix <- count.matrix[, -1]
rownames(count.matrix) <- rownames(fpkm)
colnames(count.matrix) <- colnames(fpkm)

# filtering step for low counts
keep <- rowSums(cpm(count.matrix)>1) >= 12
table(keep)
counts <- count.matrix[keep, ]

# labels for desing matrix
labels <- c("glands", "glands", rep("epithelial",5), rep("fat",3), rep("mucosa",3),
            rep("epithelial",2),"glands",rep("muscle",2),rep("glands",3),
            rep("muscle",6),rep("mucosa",5),rep("immune",4))
D <- data.frame(row.names = colnames(counts))
D$Patients <- c(2,3,2,3,1,2,3,1,2,3,1,2,3,1,2,1,1,3,2,2,3,1,2,3,3,1,3,1,2,3,1,3,1,2,3,1)
D$Organ <- labels

Patients <- as.factor(D$Patients)
Organ <- as.factor(D$Organ)
# desig
design <- model.matrix(~Patients+Organ)
Organ # ref epithelial
Patients # ref 1

data <- DGEList(counts=counts, group=D$Organ)

# estimate dispersion (edgeR)
data <- estimateDisp(data, design, robust=TRUE)

# fitting the model
fit <- glmQLFit(data, design)
colnames(fit)

# quasi-likelihood F-test
qlf <- glmQLFTest(fit, coef = 4:8)
topTags(qlf)


deg <- qlf$table
deg$id <- rownames(qlf$table)

# DEGs genes
deg.05 <- filter(deg, PValue <= 0.05)
head(deg.05)

counts.sig <- counts[deg.05$id, ]

# dendrogram
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



