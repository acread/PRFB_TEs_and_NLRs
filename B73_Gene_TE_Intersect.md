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
### Note that, as written, this returns TEs that overlap Introns+outside the gene as Intron Overlaps.... 
### I need to fix this.

```
  811  bedtools intersect -a B73TEs_Subset100.bed -b B73_Genes.bed -wo -f 1.0 > Intersect_TE100Subset_Genes_100percent.txt
  812  bedtools intersect -a B73TEs_Subset100.bed -b B73_introns.bed -wo -f 1.0 > Intersect_TE100Subset_Introns_100percent.txt
  813  bedtools intersect -a B73TEs_Subset100.bed -b B73_UTRs.bed -wo -f 1.0 > Intersect_TE100Subset_UTRs_100percent.txt
  814  bedtools intersect -a B73TEs_Subset100.bed -b B73_Genes.bed -wo > Intersect_TE100Subset_Genes_Anypercent.txt
```

### Getting rid of overlapping features \
I created a MINIMAL version of this to use instead.
```
egrep -v "ID=long_terminal_repeat|ID=lTSD|ID=rTSD|ID=lLTR|ID=rLTR|ID=repeat_region|knob|centromeric_repeat|scaf_|low_complexity|subtelomere|rDNA_intergenic_spacer_element" Zm-B73-REFERENCE-NAM-5.0.TE.gff3 > B73_MINIMAL_TEs.gff3

#this was fed into R to generate a 100bp+ bed file
B73TEs_MINIMAL_Subset100.bed
```

Repeat bedtools intersect using this minimal TE list
```
  bedtools intersect -a B73TEs_MINIMAL_Subset100.bed -b B73_Genes.bed -wo -f 1.0 > Intersect_TE100Subset_Genes_100percent.txt
  bedtools intersect -a B73TEs_MINIMAL_Subset100.bed -b B73_introns.bed -wo -f 1.0 > Intersect_TE100Subset_Introns_100percent.txt
  bedtools intersect -a B73TEs_MINIMAL_Subset100.bed -b B73_UTRs.bed -wo -f 1.0 > Intersect_TE100Subset_UTRs_100percent.txt
  bedtools intersect -a B73TEs_MINIMAL_Subset100.bed -b B73_Genes.bed -wo > Intersect_TE100Subset_Genes_Anypercent.txt
```

This method seems to have worked
I've pulled everything together into this file:  0910_100bpTEs_OverlapGenomicFeatures

Summary below:

| 	|	N	|	NLRs (n=144)            |
|-----------|-----------------------|-----------------------------------| 
| N gene and TE annotation overlap (any %)	|	91664	|	
|     Unique Tes	|	89897	|	            |
|     Unique Genes	|	22260	|	96          |
| N gene and TE annotation 100% within Gene	|	79204	|	
|     Unique Tes	|	78772	|	            |
|     Unique Genes	|	16733	|	74          |
| N gene and TE annotation 100% within INTRON	|	72154	|	
|     Unique Tes	|	71906	|	            |
|     Unique Genes	|	14410	|	63          |
| N gene and TE annotation 100% within UTR	|	2674	|	
|     Unique Tes	|	2659	|	            | 
| Unique Genes	|	2202	|	            | 
| Unique UTR	|	2229	|	16          |


## People were surprised by how high these numbers are -- suggest using only the STRUCTURAL TEs

````
grep "Method=structural" B73TEs_MINIMAL_Subset100.bed > B73TEs_MINIMAL_Structural_Subset100.bed 

bedtools intersect -a B73TEs_MINIMAL_Structural_Subset100.bed -b B73_Genes.bed -wo -f 1.0 > Intersect_TE100Struct_Genes_100percent.txt
bedtools intersect -a B73TEs_MINIMAL_Structural_Subset100.bed -b B73_introns.bed -wo -f 1.0 > Intersect_TE100Struct_Introns_100percent.txt
bedtools intersect -a B73TEs_MINIMAL_Structural_Subset100.bed -b B73_UTRs.bed -wo -f 1.0 > Intersect_TE100Struct_UTRs_100percent.txt
bedtools intersect -a B73TEs_MINIMAL_Structural_Subset100.bed -b B73_Genes.bed -wo > Intersect_TE100Struct_Genes_Anypercent.txt
````

Output brought into R to figure out how many unique elements
````R
#Any percent overlap with gene
TE_GeneAnypercent_Overlap=read_tsv("Intersect_TE100Struct_Genes_Anypercent.txt", col_names = FALSE)
Unique_Genes=unique(TE_GeneAnypercent_Overlap$X8)
Unique_TEs=unique(TE_GeneAnypercent_Overlap$X4)

#TEs completely within a gene
TE_Gene_Overlap=read_tsv("Intersect_TE100Struct_Genes_100percent.txt", col_names = FALSE)
Unique_Genes=unique(TE_Gene_Overlap$X8)
Unique_TEs=unique(TE_Gene_Overlap$X4)

#TEs completely within introns
TE_Intron_Overlap=read_tsv("Intersect_TE100Struct_Introns_100percent.txt", col_names = FALSE)
Unique_Genes=unique(TE_Intron_Overlap$X8)
Unique_TEs=unique(TE_Intron_Overlap$X4)

#TEs completely within UTRs (5' and/or 3')
TE_UTR_Overlap=read_tsv("Intersect_TE100Struct_UTRs_100percent.txt", col_names = FALSE)
Unique_Genes=unique(TE_UTR_Overlap$X8)
Unique_TEs=unique(TE_UTR_Overlap$X4)
````
Summary:

| 	|	N	|
|-----------|-----------------------|
| gene and TE annotation overlap (any %)	|	
|     Unique Tes	|	8379	|
|     Unique Genes	|	7619	|
| gene and TE annotation 100% within Gene	|	
|     Unique Tes	|	5183	|
|     Unique Genes	|	3832	|
| gene and TE annotation 100% within INTRON	|	
|     Unique Tes	|	4746	|
|     Unique Genes	|	3498	|
| gene and TE annotation 100% within UTR	|	
|     Unique Tes	|	117	| 
|     Unique Genes	|	117         |

