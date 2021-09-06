## Intersecting Gene and TE coordinates

While I work on this, Nicole and Noelle are asking related questions, but are \
looking at percent of gene and gene-proximal regions are annotated as TE \

First I took the most recent TE annotation and removed TEs < 100bp
```
B73TEs_subset100.bed
```

For gene data, used a primary gene list made from the most recent NAM annotation
```
Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.primary.gff3
```

I used grep commands to get gene and UTR file coordinates \
I used R "GenomicFeatures" to get an intron file\
### NOTE that this code includes the TE filtering and formatting the gene+UTR gffs as beds
```R
library("dotPlotly")
library("plotly")
library("optparse")
library("tidyverse")

setwd("/Users/read0094/Desktop/Maize")

BiocManager::install("GenomicFeatures")
library("GenomicFeatures")

x=makeTxDbFromGFF("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.primary.gff3")
all.introns = intronicParts(x)

y=as.data.frame(all.introns)
y=y %>% select(1,2,3,8)

write.table(y, file='B73_introns.bed', quote=FALSE, sep='\t', col.names = FALSE, 
            row.names = FALSE)

B73TEs=read_tsv("Zm-B73-REFERENCE-NAM-5.0.TE.txt", col_names = FALSE)
B73TEs$len=B73TEs$X5-B73TEs$X4
B73TEs_Subset50=subset(B73TEs, len > 50)
B73TEs_Subset100=subset(B73TEs, len > 100)
B73TEs_Subset1000=subset(B73TEs, len > 1000)

B73TEs_Subset50=B73TEs_Subset50 %>% select(1,4,5,9)
B73TEs_Subset100=B73TEs_Subset100 %>% select(1,4,5,9)
B73TEs_Subset1000=B73TEs_Subset1000 %>% select(1,4,5,9)


write.table(B73TEs_Subset50, file="B73TEs_Subset50.bed", quote=FALSE, sep='\t',
            col.names=FALSE, row.names=FALSE)
write.table(B73TEs_Subset100, file="B73TEs_Subset100.bed", quote=FALSE, sep='\t',
            col.names=FALSE, row.names=FALSE)
write.table(B73TEs_Subset1000, file="B73TEs_Subset1000.bed", quote=FALSE, sep='\t',
            col.names=FALSE, row.names=FALSE)


B73_UTRS=read_tsv("Zm-B73-REFERENCE-NAM-5.0-UTRs.gff3", col_names = FALSE)
B73_UTRS$cat=paste(B73_UTRS$X9,B73_UTRS$X3, sep="-")
B73_UTRS=B73_UTRS %>% select(1,4,5,10)
write.table(B73_UTRS, file="B73_UTRs.bed", quote=FALSE, sep='\t',
            col.names=FALSE, row.names=FALSE)

B73_Genes=read_tsv("Zm-B73-REFERENCE-NAM-5.0-gene.gff3", col_names= FALSE)
B73_Genes=B73_Genes %>% select(1,4,5,9)
write.table(B73_Genes, file="B73_Genes.bed", quote=FALSE, sep='\t',
            col.names=FALSE, row.names=FALSE)
```

Using bedtools to get overlaps of TEs and Genes
Here are the commands I used to get overlaps \
### Note that, as written, this returns TEs that overlap Introns+outside the gene as Intron Overlaps.... \
### I need to fix this.

```
  811  bedtools intersect -a B73TEs_Subset100.bed -b B73_Genes.bed -wo -f 1.0 > Intersect_TE100Subset_Genes_100percent.txt
  812  bedtools intersect -a B73TEs_Subset100.bed -b B73_introns.bed -wo -f 1.0 > Intersect_TE100Subset_Introns_100percent.txt
  813  bedtools intersect -a B73TEs_Subset100.bed -b B73_UTRs.bed -wo -f 1.0 > Intersect_TE100Subset_UTRs_100percent.txt
  814  bedtools intersect -a B73TEs_Subset100.bed -b B73_Genes.bed -wo > Intersect_TE100Subset_Genes_Anypercent.txt
```

