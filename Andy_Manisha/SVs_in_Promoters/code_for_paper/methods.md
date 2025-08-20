## Methods

### RNAseq alignment to B73v5
The B73v5 genome was indexed using STAR – the genome and .gtf were pulled from MSI `/home/springer/shared/ns-genome/Zmays_B73v5`
`````
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=stargenomegenerate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu

module load star/2.5.3a

#mkdir /scratch.global/read0094/NAM_RNAseq/B73v5_index_STAR

STAR --runMode genomeGenerate --runThreadN 20 --genomeDir /scratch.global/read0094/NAM_RNAseq/B73v5_index_STAR --genomeFastaFiles /scratch.global/read0094/NAM_RNAseq/10.fasta --sjdbGTFfile /scratch.global/read0094/NAM_RNAseq/10.gtf
`````
RNAseq data from the NAM publication `table NAM_10tissue_RNAseq_SRA.xlsx`\
RNAseq data was trimmed with Trim_Galore and aligned to the B73v5 genome using STAR. Gene counts were generated with HT-SEQ and normalized in R using DESEQ2.
`````
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32gb
#SBATCH -J NAM_RNAseq_AR
#SBATCH --array=1-2


#newdir=${SLURM_JOB_ID}

myFile=`sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch.global/read0094/NAM_RNAseq/NAM_ERR_list.txt`


########## Modules #################

module load fastqc/0.11.5
module load cutadapt
module load star/2.5.3a
module load htseq/0.7.2
module load python3

########## Set up dirs #################
#mkdir /scratch.global/read0094/NAM_RNAseq/trimmed
#mkdir /scratch.global/read0094/NAM_RNAseq/fastqc



##Download RNAseq
#/home/springer/read0094/Software/enatools/enaBrowserTools-1.1.0/python3/enaDataGet -f fastq -d /scratch.global/read0094/NAM_RNAseq ${myFile}



##Trim_Galore
/home/springer/read0094/Software/trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir fastqc" -o trimmed/ \
 --paired /scratch.global/read0094/NAM_RNAseq/${myFile}/${myFile}_1.fastq.gz /scratch.global/read0094/NAM_RNAseq/${myFile}/${myFile}_2.fastq.gz
 
##STAR align
STAR --genomeDir /scratch.global/read0094/NAM_RNAseq/B73v5_index_STAR \
     --readFilesCommand zcat --runThreadN 8 \
     --readFilesIn /scratch.global/read0094/NAM_RNAseq/trimmed/${myFile}_1_val_1.fq.gz /scratch.global/read0094/NAM_RNAseq/trimmed/${myFile}_2_val_2.fq.gz \
     --outFileNamePrefix ${myFile}_STAR \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \ --outSAMattributes Standard

##HTseq - align to genes
mkdir /scratch.global/read0094/NAM_RNAseq/htseq_out
htseq-count -f bam -s no -t gene -i ID -m union -a 20 /scratch.global/read0094/NAM_RNAseq/${myFile}_STARAligned.sortedByCoord.out.bam \
  /scratch.global/read0094/NAM_RNAseq/10.gff > /scratch.global/read0094/NAM_RNAseq/htseq_out/htseq_${myFile}.txt

##Remove the junk
rm -R ${myFile}
rm -R *STARtmp
`````

<details>
	
  <summary>Rcode</summary>
  
`````R
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
`````     
</details>


### Identifying genic Structural Variants
### Note that the categorization strategy changed - new strat includes 2k upstream of translational start as "promoter" - the code and description below do NOT reflect this
Structural variants between B73 and each of the 25 additional NAM founder lines were identified via Anchorwave genome alignment in Munasinghe 2023 \
These SVs are defined as 50bp or greater in/dels that map to a single bp in one cultivar
I'm categorizing 'genic' SVs as those that either overlap the canonical gene model, or overlap the 1kb upstream of the TSS \
Note that I am not parsing 5'UTR - it is categorized as either intron/exon \
First I generated a sorted .bed file with coordinates of promoters, introns, exons. This was intersected with each of the B73 v NAM SV bam files \
This was repeated for each genome to generate 3 distinct output files: \
1. Tx303_All_SVs_comp_B73_details.txt
`````
V6	V7_simple	V1	V2	V3	V4	V15	V8	V9	V10	V19	V20	V21	gene_feature
Zm00001eb000050	intron	chr1	109952	113761	3810	chr1_AW_BlockID_3	chr1	110302	111609	Category_5	1307	B73_Ins_Rel_Tx303	Zm00001eb000050_intron
`````
2. Tx303All_B73_genes_with_SVs_reltoTx303.txt
`````
V6	B73_Prom_Del	B73_Intron_Del	B73_Exon_Del	B73_Prom_Ins	B73_Intron_Ins	B73_Exon_Ins	sum
Zm00001eb000010	0	0	0	0	0	0	0
Zm00001eb000020	0	0	0	0	0	0	0
`````
3.  Tx303single_prom_SVs_simple.txt
`````
V6	V1	V2	V3	V15	V8	V9	V10	V19	V21	orientation
Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_25	chr1	322874	328037	Category_5	B73_Ins_Rel_Tx303	+
Zm00001eb000220	chr1	1017150	1018150	chr1_AW_BlockID_121	chr1	1015278	1018118	Category_4	B73_Ins_Rel_Tx303
`````

<details>
	
  <summary>Rcode</summary>
  
`````R
##### setup #####
library("optparse")
library("tidyverse")
library("stringr")
library("data.table")
library("GenomicFeatures")
library("bedtoolsr")
library("MetBrewer")
library("ggpubr")
library("gridExtra")

################
################
# This script outputs SV data for B73 vs each NAM line - the outputs are a file
# with all genes with SVs with columns for counts per feature (PIED)
# I am doing a little bit of filtering - "filtered" means the gene is at least 80%
# alignable at the exon level, and less than 95% alignable at the gene +/- 1kb
# level (these filters can change) -- additionally the gene must not overlap
# 'missing data' or 'unalignable' blocks at all (=0)
setwd("/Users/read0094/Desktop/Maize")

ALL_AW_blocks=fread("all_AWBlocks_CategoryAssignment.tsv")

#get rid of the downstream category
#B73_PIED_bed_sort_noDown=subset(B73_PIED_bed_sort, V7 != "Downstream")

#genome="B97"
#genome="CML103"
#genome="CML228"
#genome="CML247"
#genome="CML277"
#genome="CML322"
#genome="CML333"
#genome="CML52"
#genome="CML69"
#genome="HP301"
#genome="Il14H"
#genome="Ki11"
#genome="Ki3"
#genome="Ky21"
#genome="M162W"
#genome="M37W"
#genome="Mo18W"
#genome="Ms71"
#genome="NC350"
#genome="NC358"
#genome="Oh43"
#genome="Oh7B"
#genome="P39"
#genome="Tx303"
genome="Tzi8"

B73_PIED_bed_sort=fread("B73_PIED_bed_sort.txt")

#setwd("/Users/read0094/Desktop/Maize")
dir.create(paste(genome,"_noDown",sep=""))
#setwd(paste("/Users/read0094/Desktop/Maize/",genome,"_noDown",sep=""))


NAM_AW_blocks=subset(ALL_AW_blocks, Lineage_Comp==genome)
NAM_AW_blocks$ASM_chr=paste("chr",NAM_AW_blocks$ASM_chr, sep="")

NAM_AW_blocks_B73INS=subset(NAM_AW_blocks, Block_Type=="structural_insertion_inB73")
NAM_AW_blocks_NAMINS=subset(NAM_AW_blocks, Block_Type==paste("structural_insertion_in",genome,sep=""))

#####INSERTIONS
All_B73_SV_Inserts=bedtoolsr::bt.intersect(a=B73_PIED_bed_sort,b=NAM_AW_blocks_B73INS, wao=TRUE)
All_B73_SV_Inserts=subset(All_B73_SV_Inserts, V9 != -1)
All_B73_SV_Inserts=dplyr::filter(All_B73_SV_Inserts, !grepl(",",V6))
All_B73_SV_Inserts$V7_simple=gsub('_.*', '', All_B73_SV_Inserts$V7)
All_B73_SV_Inserts_a=All_B73_SV_Inserts %>% dplyr::select(6,21,1,2,3,4,15,8,9,10,19,20)
All_B73_SV_Inserts_a$V21=paste("B73_Ins_Rel_",genome,sep="")
All_B73_SV_Inserts_a$gene_feature=paste(All_B73_SV_Inserts_a$V6, All_B73_SV_Inserts_a$V7_simple, sep="_")


#####DELETIONS
All_B73_SV_Deletions=bedtoolsr::bt.intersect(a=B73_PIED_bed_sort,b=NAM_AW_blocks_NAMINS, wao=TRUE)
All_B73_SV_Deletions=subset(All_B73_SV_Deletions, V9 != -1)
All_B73_SV_Deletions=dplyr::filter(All_B73_SV_Deletions, !grepl(",",V6))
All_B73_SV_Deletions$V7_simple=gsub('_.*', '', All_B73_SV_Deletions$V7)
All_B73_SV_Deletions_a=All_B73_SV_Deletions %>% dplyr::select(6,21,1,2,3,4,15,8,9,10,19,20)
All_B73_SV_Deletions_a$V21=paste("B73_Del_Rel_",genome,sep="")
All_B73_SV_Deletions_a$gene_feature=paste(All_B73_SV_Deletions_a$V6, All_B73_SV_Deletions_a$V7_simple, sep="_")

#Combine the Inserts and Deletion SVs
All_B73_SVs=dplyr::bind_rows(All_B73_SV_Inserts_a, All_B73_SV_Deletions_a)

write_tsv(All_B73_SVs, paste("/Users/read0094/Desktop/Maize/",genome,"_noDown/",genome,"_All_SVs_comp_B73_details.txt",sep=""))

#subset by P,I,E feature? 
Prom_B73_SV_Inserts=subset(All_B73_SV_Inserts_a, V7_simple=="Promoter")
Prom_B73_SV_Inserts_simple=Prom_B73_SV_Inserts %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Prom_Ins=n())
Intron_B73_SV_Inserts=subset(All_B73_SV_Inserts_a, V7_simple=="intron")
Intron_B73_SV_Inserts_simple=Intron_B73_SV_Inserts %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Intron_Ins=n())
Exon_B73_SV_Inserts=subset(All_B73_SV_Inserts_a, V7_simple=="exon")
Exon_B73_SV_Inserts_simple=Exon_B73_SV_Inserts %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Exon_Ins=n())

#probably better to use a full gene list here as the thing to join to
B73_genes=fread("B73_full_gene_list.txt")
All_B73_Inserts_join=dplyr::full_join(B73_genes,Prom_B73_SV_Inserts_simple)
All_B73_Inserts_join=dplyr::full_join(All_B73_Inserts_join, Intron_B73_SV_Inserts_simple)
All_B73_Inserts_join=dplyr::full_join(All_B73_Inserts_join, Exon_B73_SV_Inserts_simple)

#subset by P,I,E,D feature? 
Prom_B73_SV_Deletions=subset(All_B73_SV_Deletions_a, V7_simple=="Promoter")
Prom_B73_SV_Deletions_simple=Prom_B73_SV_Deletions %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Prom_Del=n())
Intron_B73_SV_Deletions=subset(All_B73_SV_Deletions_a, V7_simple=="intron")
Intron_B73_SV_Deletions_simple=Intron_B73_SV_Deletions %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Intron_Del=n())
Exon_B73_SV_Deletions=subset(All_B73_SV_Deletions_a, V7_simple=="exon")
Exon_B73_SV_Deletions_simple=Exon_B73_SV_Deletions %>% dplyr::group_by(V6) %>% dplyr::summarize(B73_Exon_Del=n())

#probably better to use a full gene list here as the thing to join to
All_B73_Deletions_join=dplyr::full_join(B73_genes,Prom_B73_SV_Deletions_simple)
All_B73_Deletions_join=dplyr::full_join(All_B73_Deletions_join, Intron_B73_SV_Deletions_simple)
All_B73_Deletions_join=dplyr::full_join(All_B73_Deletions_join, Exon_B73_SV_Deletions_simple)

All_B73_SVs_count=full_join(All_B73_Deletions_join,All_B73_Inserts_join,by="V6")

#replace NA with zeros and sum rows
All_B73_SVs_count[is.na(All_B73_SVs_count)] = 0
All_B73_SVs_count$sum=rowSums(All_B73_SVs_count[,2:7])
write_tsv(All_B73_SVs_count, paste("/Users/read0094/Desktop/Maize/",genome,"_noDown/",genome,"All_B73_genes_with_SVs_relto",genome,".txt",sep=""))

All_B73_single_SVs=subset(All_B73_SVs_count, sum == 1)
colSums(All_B73_single_SVs[,-1])

rm(list=setdiff(ls(), "ALL_AW_blocks"))

#B97
#B73_Prom_Del B73_Intron_Del   B73_Exon_Del   B73_Prom_Ins B73_Intron_Ins   B73_Exon_Ins 
#1553            731            661           2620            681            411 
#sum 
#6657 
`````
</details>

	
### Expression of Promoter SV genes across tissues
Identify Promoter SVs that overlap TEs - Look at gene expression within a tissue, how many are up? down?
	
### Expression Plasticity Across Tissues
For each gene across tissues, calculate the mean and variance \
Group genes by presence/absence of TE-SV and T-test to compare \
	






