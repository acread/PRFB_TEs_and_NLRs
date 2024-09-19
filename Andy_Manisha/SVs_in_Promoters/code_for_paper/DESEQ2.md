R code
`````
library(DESeq2)
library(gplots)
library(GenomeInfoDb)
library(data.table)
#library(tidyverse)

test="asdg"
write.csv(test, "check.txt")

### Set working directory and read in matrix data ###
setwd("/home/springer/read0094/NAM_RNAseq/DESEQ")
counts <- read.delim(file="NAM_RNAseq_htseq_cat.txt", head=T, sep="\t", row.names=1)

### Set up the data
# Create dataframe with tissues as conditions 
#design = data.frame(row.names=colnames(counts), condition=c(""
design=fread("NAM_RNAseq_design.txt", header=TRUE)
condition = design$condition

# Create count data set varialble within DESeq2
cds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = design,
                              design = ~ condition)

# Normalize the data
cds <- estimateSizeFactors(cds)
sizeFactors(cds)

cds <- estimateDispersions(cds)

# Create object with normalized values
norm <- counts(cds, normalized = TRUE)

write.csv(norm, "/scratch.global/read0094/NAM_RNAseq/DESEQnorm_NAM_RNAseq_counts.txt")
