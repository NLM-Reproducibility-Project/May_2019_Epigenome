# import FPKM data
setwd("~/Documents/NLM_Repro/NLM_Repro_Workshop")
raw.fpkm <- read.csv("data/RNA_matrix_FPKM.tsv", header = T, sep = '\t')


fpkm <- raw.fpkm[ , 10:45]
rownames(fpkm) <- raw.fpkm$gene_id
# remove last 
tmp <- colnames(fpkm)
tmp <- sub("_FPKM", "", tmp)
colnames(fpkm) <- tmp
rm(tmp)
dim(fpkm)


locus <- read.csv("data/location.csv", header = F, stringsAsFactors = F)
locus <- locus$V1
transcript.length <- c()

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




names <- colnames(fpkm)

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

count.matrix <- convertFPKM.to.counts(fpkm, transcript.length, total.mapped)
count.matrix <- count.matrix[, -1]
rownames(count.matrix) <- rownames(fpkm)
colnames(count.matrix) <- colnames(fpkm)






