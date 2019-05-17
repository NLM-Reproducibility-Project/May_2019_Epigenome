library(EDASeq)
# read table contains fpkm data of donor samples
tissue_fpkm <- read.table(file = "RNA_matrix_FPKM.tsv",header = TRUE, stringsAsFactors = FALSE)

# Get coordinates of genes
chr_start <- gsub(".*:","",tissue_fpkm$locus)
chr_start <- gsub("-.*","",chr_start)
chr_end <- gsub(".*-","",tissue_fpkm$locus)

# Gene length
tissue_fpkm$length <- abs(as.numeric(chr_start) - as.numeric(chr_end) )

# Get exon length of the genes which will be used to transform FPKM back into raw count
ensembl_gene_id <- gsub("\\..*","",tissue_fpkm$gene_id)

# Are there duplicates in the ensembl gene id? 
sum(duplicated(ensembl_gene_id)) # [1] 0 # No duplicates
# Get exon length and GC content
gene_length_and_gc_content <- getGeneLengthAndGCContent(ensembl_gene_id, "hsa")

# Exon length of genes
tissue_fpkm$exon_length <- gene_length_and_gc_content[,1]

# For ensembl gene id discarded, which exon_length with be 'NA'. Use gene length as exon length
tissue_fpkm$exon_length[is.na(tissue_fpkm$exon_length)] <- tissue_fpkm$length[is.na(tissue_fpkm$exon_length)]

# discard Placenta, IMR90 and H1 data and keep only data from donor samples
# replace 'length' column with 'exon_length'
donor_tissue_fpkm <- tissue_fpkm[,c(1:7,51,9:45)]
# Get rid of "_FPKM" in the colnames
colnames(donor_tissue_fpkm) <- gsub("_FPKM","",colnames(donor_tissue_fpkm))

# Get total reads for each sample so that raw count can be calculated from FPKM
# The total_reads_mapped_each_sample.csv is part of Sup. Table 1.
total_reads_per_sample <- read.csv(file = "total_reads_mapped_each_sample.csv", stringsAsFactors = FALSE)

# After manually checking the file;
# one abbreviation maybe wrong! "BL-3" should be "BL-1"
total_reads_per_sample$Abbreviation[5] # "BL-3"
# Correct the error
total_reads_per_sample$Abbreviation[5] <- "BL-1"

# Check tissue names(abbreivations) agree between total_reads_per_sample and donor_tissue_fpkm
colnames(donor_tissue_fpkm)[-(1:9)] %in% total_reads_per_sample$Abbreviation

# NOTE: the sample abbreviation is noted as "EG_3" in tissue_fpkm dataframe
#       while the abbreviation in total_reads_per_sample is noted as "EG-3"
#       make them the same first: all use '_'!
total_reads_per_sample$Abbreviation <- gsub('-','_',total_reads_per_sample$Abbreviation)

# All sample names in donor_tissue_fpkm also in total_reads_per_sample?
all(colnames(donor_tissue_fpkm)[-(1:9)] %in% total_reads_per_sample$Abbreviation)

# Order donor_tissue_fpkm by tissue abbreviates
ordered_donor_tissue_fpkm <- cbind(donor_tissue_fpkm[,1:9], donor_tissue_fpkm[, sort(colnames(donor_tissue_fpkm)[-(1:9)])] )

# Order total_reads_per_sample by tissue abbreviates
ordered_total_reads_per_sample <- total_reads_per_sample[order(total_reads_per_sample$Abbreviation),]

# Transfer FPKM data into real raw counts 
# RPKM = (ExonMappedReads * 10^9 ) / (TotalMappedReads * ExonLength), 
# So for paired-end reads here, read_counts = 2* FPKM*(TotalMappedReads*ExonLength)/10^9

# Construct a matrix containing total reads of each sample to do calculation later
sample_total_read_matrix <- matrix(rep(ordered_total_reads_per_sample$Mapped_RNA_seq.Reads, nrow(ordered_donor_tissue_fpkm) ), 
                                   nrow=nrow(ordered_donor_tissue_fpkm), byrow=TRUE)

# calculate raw counts; use round() function to get integers
donor_raw_counts <- round(2* ordered_donor_tissue_fpkm[,-(1:9)] * (ordered_donor_tissue_fpkm$exon_length) * sample_total_read_matrix/10^9)

# Add back annotation data
donor_raw_counts_with_annotation <- cbind(ordered_donor_tissue_fpkm[,1:9], donor_raw_counts)
# Write the raw count into tsv file
write.table(donor_raw_counts_with_annotation, file="integer_raw_counts_of_genes_among_donor_tissues.tsv", sep='\t', row.names = FALSE)
