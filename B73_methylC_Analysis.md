### Analyzing 3 methylseq data files mostly for Zach

Zach is working on a meta-analysis of gene expression changes in maize following heat stress. \
He is interested in taking a look at metaplots around genes that respond to heat and has asked if I could apply the \
Springer-lab pipeline to help him out... so I'm doing that.  I also think it's a great opportunity for me to dust \
off these skills and make a version that takes single gene IDs as inputs to generate faceted plots 

### Step 0 - download the data
`````
#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH --tmp=4g
#SBATCH --job-name=download
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR534/006/ERR5347656/ERR5347656_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR534/006/ERR5347656/ERR5347656_1.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR534/004/ERR5347654/ERR5347654_2.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR534/005/ERR5347655/ERR5347655_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR534/005/ERR5347655/ERR5347655_2.fastq.gz
`````

### Step 1 - trim reads
`````
#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=maize_mC_trim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu


########## Modules #################

module load fastqc/0.11.5
module load cutadapt/1.8.1

########## Set up dirs #################

#make trimmed folder
trimmedfolder=analysis/trimmed
mkdir -p $trimmedfolder

fastqcfolder=analysis/fastqc
mkdir -p $fastqcfolder

#uncompress reads because trim_galore throws the error `gzip: stdout: Broken pipe` if I input .gz files
gunzip reads/ERR5347654_1.fastq.gz
gunzip reads/ERR5347654_2.fastq.gz
gunzip reads/ERR5347655_1.fastq.gz
gunzip reads/ERR5347655_2.fastq.gz
gunzip reads/ERR5347656_1.fastq.gz
gunzip reads/ERR5347656_2.fastq.gz


########## Run #################
/home/springer/read0094/Software/trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder --paired reads/ERR5347654_1.fastq reads/ERR5347654_2.fastq
/home/springer/read0094/Software/trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder --paired reads/ERR5347655_1.fastq reads/ERR5347655_2.fastq
/home/springer/read0094/Software/trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir $fastqcfolder" -o $trimmedfolder --paired reads/ERR5347656_1.fastq reads/ERR5347656_2.fastq

#compress original reads again
gzip reads/ERR5347654_1.fastq
gzip reads/ERR5347654_2.fastq
gzip reads/ERR5347655_1.fastq
gzip reads/ERR5347655_2.fastq
gzip reads/ERR5347656_1.fastq
gzip reads/ERR5347656_2.fastq

fi

echo Done trimming
`````

### Step 2 - align to B73v5 with bsmap

`````
#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=step2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu


########## Set up dirs #################
cd analysis
mkdir -p bsmapped

########## Run #################

/home/springer/read0094/Software/bsmap-2.74/bsmap \
-a trimmed/ERR5347654_1_val_1.fq \
-b trimmed/ERR5347654_2_val_2.fq \
-d /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta \
-o bsmapped/ERR5347654.sam \
-v 5 \
-r 0 \
-p 8 \
-q 20 \
-A $adapter_seq

bsmap \
-a trimmed/ERR5347655_1_val_1.fq \
-b trimmed/ERR5347655_2_val_2.fq \
-d /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta \
-o bsmapped/ERR5347655.sam \
-v 5 \
-r 0 \
-p 8 \
-q 20 \
-A $adapter_seq

bsmap \
-a trimmed/ERR5347656_1_val_1.fq \
-b trimmed/ERR5347656_2_val_2.fq \
-d /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta \
-o bsmapped/ERR5347656.sam \
-v 5 \
-r 0 \
-p 8 \
-q 20 \
-A $adapter_seq
`````

### Step 3 - indexing and sorting sam/bams

`````
#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=2021_02_05_04_filter
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu

########## Modules #################

#bsmap requires samtools < 1.0.0
# module load samtools/0.1.18 # no longer an available module
module load samtools

########## Set up dirs #################
cd analysis/bsmapped

########## Run #################

# make bam
samtools view -b ERR5347654.sam -o ERR5347654.bam
samtools view -b ERR5347655.sam -o ERR5347655.bam
samtools view -b ERR5347656.sam -o ERR5347656.bam

# sort by read name (needed for fixsam)
samtools sort -n ERR5347654.bam -o ERR5347654_nameSrt.bam
samtools sort -n ERR5347655.bam -o ERR5347655_nameSrt.bam
samtools sort -n ERR5347656.bam -o ERR5347656_nameSrt.bam

# fix mate pairs
samtools fixmate ERR5347654_nameSrt.bam ERR5347654_nameSrt_fixed.bam
samtools fixmate ERR5347655_nameSrt.bam ERR5347655_nameSrt_fixed.bam
samtools fixmate ERR5347656_nameSrt.bam ERR5347656_nameSrt_fixed.bam

# co-ordinate sort
samtools sort ERR5347654_nameSrt_fixed.bam -o ERR5347654_sorted.bam
samtools sort ERR5347655_nameSrt_fixed.bam -o ERR5347655_sorted.bam
samtools sort ERR5347656_nameSrt_fixed.bam -o ERR5347656_sorted.bam

# index
samtools index ERR5347654_sorted.bam
samtools index ERR5347655_sorted.bam
samtools index ERR5347656_sorted.bam

# remove intermediate files
rm ERR5347654.sam ERR5347654_nameSrt.bam ERR5347654_nameSrt_fixed.bam ERR5347654.bam
rm ERR5347655.sam ERR5347655_nameSrt.bam ERR5347655_nameSrt_fixed.bam ERR5347655.bam
rm ERR5347656.sam ERR5347656_nameSrt.bam ERR5347656_nameSrt_fixed.bam ERR5347656.bam

echo Done fixing
`````
### NOTE: the samples can (should?) be merged after step 3, but I was having memory problems on MSI - this would be very helpful

### Step 4 - mark duplicates, keep proper pairs, clip overlapping reads

`````
#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=step4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu


########## Modules #################

#module load python2/2.7.8
module load java
#module load bedtools
module load samtools
module load bamtools

########## Set up dirs #################

cd analysis
mkdir -p bsmapped_filtered

########## Run #################

#ERR5347654
        samtools collate -o bsmapped/ERR5347654_collate.bam bsmapped/ERR5347654_sorted.bam
        samtools fixmate -m bsmapped/ERR5347654_collate.bam bsmapped/ERR5347654_sorted_fixmate.bam
        samtools sort -o bsmapped/ERR5347654_sorted_collate.bam bsmapped/ERR5347654_sorted_fixmate.bam
        samtools markdup bsmapped/ERR5347654_sorted_collate.bam bsmapped_filtered/ERR347654_sorted_MarkDup.bam


        # keep properly paired reads using bamtools package
        # note that some reads marked as properly paired by bsmap actually map to different chromosomes
        bamtools filter \
        -isMapped true \
        -isPaired true \
        -isProperPair true \
        -in bsmapped_filtered/ERR347654_sorted_MarkDup.bam \
        -out bsmapped_filtered/ERR5347654_sorted_MarkDup_pairs.bam

        # clip overlapping reads using bamUtils package
        /home/springer/read0094/Software/bamUtil/bin/bam clipOverlap \
        --in bsmapped_filtered/ERR5347654_sorted_MarkDup_pairs.bam \
        --out bsmapped_filtered/ERR5347654_sorted_MarkDup_pairs_clipOverlap.bam \
        --stats

        #index bam
        # index
        samtools index bsmapped_filtered/ERR5347654_sorted_MarkDup_pairs_clipOverlap.bam
`````

### Step 5 - Sumarize methylation 

This is a pretty complex step - it takes the input in the bsmapped_filtered folder and generates several outputs: \
It makes a bed file and splits this into files for each mC context \
it makes bedgraph files and converts to bigwig for genome browser viewing (not sure I need this for this, but keeping it) \
Pete's version calculates a conversion rate based on the chloroplast sequence - requires a 'Pt' chromosome - I'm dropping this

**NOTE: It looks like Pete has switched to MethylDackel for this -- https://github.com/dpryan79/MethylDackel -- Look into this before running**

`````
#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=4
#SBATCH --mem=40g
#SBATCH --tmp=10g
#SBATCH --job-name=step5_54
#SBATCH --mail-type=ALL
#SBATCH --mail-user=read0094@umn.edu

########## Modules #################

#module load python2/2.7.8
module load python2/2.7.12_anaconda4.2
#module load java
module load bedtools


#make adaligned folder bsmaped
cd analysis
mkdir -p BSMAPratio
mkdir -p TempOut
#mkdir -p OnTargetCoverage
mkdir -p ConversionRate

########## Run #################
        # The required input is all in the folder bsmapped_filtered from the prior step
        ########################
        # extract methylation information using bsmap tool methratio.py
        python /home/springer/read0094/Software/bsmap-2.74/methratio.py \
        -o BSMAPratio/ERR347654_methratio.txt \
        -d /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta \
        -u \
        -z \
        -r bsmapped_filtered/ERR347654_sorted_MarkDup.bam

        # Well 2nt resolution makes a proper zero-based coordinate BED file... (this is what Qing did).
        awk_make_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        print $1, $2-1, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/))
                        print $1, $2-1, $2, $3, "CHG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/))
                        print $1, $2-1, $2, $3, "CHH", $5, $6, $7, $8, $9, $10, $11, $12;
else
                        print $1, $2-1, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                '

        awk_make_subcontext_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        print $1, $2-1, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^CTG../ ) || ($3=="+" &&  $4~/^..CAG/))
                        print $1, $2-1, $2, $3, "CAG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^CGG../ ) || ($3=="+" &&  $4~/^..CCG/))
                        print $1, $2-1, $2, $3, "CCG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GAC../ ) || ($3=="+" &&  $4~/^..CTG/))
                        print $1, $2-1, $2, $3, "CTG", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TTG../ ) || ($3=="+" &&  $4~/^..CAA/))
                        print $1, $2-1, $2, $3, "CAA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GTG../ ) || ($3=="+" &&  $4~/^..CAC/))
                        print $1, $2-1, $2, $3, "CAC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^ATG../ ) || ($3=="+" &&  $4~/^..CAT/))
                        print $1, $2-1, $2, $3, "CAT", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TGG../ ) || ($3=="+" &&  $4~/^..CCA/))
                        print $1, $2-1, $2, $3, "CCA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GGG../ ) || ($3=="+" &&  $4~/^..CCC/))
                        print $1, $2-1, $2, $3, "CCC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^AGG../ ) || ($3=="+" &&  $4~/^..CCT/))
                        print $1, $2-1, $2, $3, "CCT", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^TAG../ ) || ($3=="+" &&  $4~/^..CTA/))
                        print $1, $2-1, $2, $3, "CTA", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^GAG../ ) || ($3=="+" &&  $4~/^..CTC/))
                        print $1, $2-1, $2, $3, "CTC", $5, $6, $7, $8, $9, $10, $11, $12;
                      else if(($3=="-" && $4~/^AAG../ ) || ($3=="+" &&  $4~/^..CTT/))
                        print $1, $2-1, $2, $3, "CTT", $5, $6, $7, $8, $9, $10, $11, $12;
                else
                        print $1, $2-1, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                '

        awk -F$"\\t" "$awk_make_bed" \
        "BSMAPratio/ERR347654_methratio.txt" > "BSMAPratio/ERR347654_BSMAP_out.txt"

        awk -F$"\\t" "$awk_make_subcontext_bed" \
        "BSMAPratio/ERR347654_methratio.txt" > "BSMAPratio/ERR347654_BSMAP_out_subcontext.txt"

        ########################
        #For genome browser

        # begGraph ratio files for tdfs
        awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $8/$9*100, $5
        }
        '

        # split bedGraph by contex
        awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_"$5".bedGraph"
        }
        '

        # split bedGraph by sub-contex
        # only difference in the output file suffix (so that we make CG files for context and sub-contex: sainty check - they should be the same)
        awk_make_bedGraph_subcontext='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_subcontext_"$5".bedGraph"
        }
        '

        #pipe bedGraph to split by context (use dash to read from sdtin)
        # per context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/ERR347654_BSMAP_out.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

        # per sub-context
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/ERR347654_BSMAP_out_subcontext.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_subcontext" -

        #Make bigWigs per context
        /home/springer/shared/share_software/bedGraphToBigWig "BSMAPratio/ERR347654_BSMAP_out_CG.bedGraph" /home/springer/read0094/SetariaStuff/chrLength.txt \
        "BSMAPratio/ERR347654_BSMAP_out_CG.bigWig"
        /home/springer/shared/share_software/bedGraphToBigWig "BSMAPratio/ERR347654_BSMAP_out_CHG.bedGraph" /home/springer/read0094/SetariaStuff/chrLength.txt \
        "BSMAPratio/ERR347654_BSMAP_out_CHG.bigWig"
        /home/springer/shared/share_software/bedGraphToBigWig "BSMAPratio/ERR347654_BSMAP_out_CHH.bedGraph" /home/springer/read0094/SetariaStuff/chrLength.txt \
        "BSMAPratio/ERR347654_BSMAP_out_CHH.bigWig"

        #remove bedGraph it is large and not really required
        # keep bigWigs
        rm -rv BSMAPratio/ERR347654*.bedGraph

`````

### Step 6 - Get average methylation across tiles - NOTE: this is step 07 for Pete C.

`````
Fill in R and .sh scripts for this
`````

### Visualizing tile-based mC - FileZilla'd ERR347654_merge_gene_metaplot.df.txt onto my machine

`````
library(ggplot2)
library(fields)
library(janitor)

setwd("/Users/read0094/Desktop/Maize/ForZach/")

#data for all tiles
q1 <- read.table(file="ERR347654_merge_gene_metaplot.df.txt", header=T, sep="\t")

#subsetting for tiles near genes
data.sub1=subset(q1,q1$distance < 2100 & q1$distance > -2000)
write.csv(data.sub1,"data.sub1.txt")

#PETEs code to average for all genes - I didn't do this but Zach can

#get 100 bins across the relative gene distance and get stats on your data column
#stats.test1.cg=stats.bin(data.sub1$distance,data.sub1$cg,N=42)
#stats.test1.chg=stats.bin(data.sub1$distance,data.sub1$chg,N=42)
#stats.test1.chh=stats.bin(data.sub1$distance,data.sub1$chh,N=42)

#p.1.cg=cbind(matrix(stats.test1.cg$centers,ncol=1),stats.test1.cg$stats["mean",])
#p.1.chg=cbind(matrix(stats.test1.chg$centers,ncol=1),stats.test1.chg$stats["mean",])
#p.1.chh=cbind(matrix(stats.test1.chh$centers,ncol=1),stats.test1.chh$stats["mean",])

#NOTE: genes can be split by some new column (NLR vs nonNLR, or Heat responsive vs not and \
#graphed seperately prior to this graphing step -- also can convert to ggplot

#pdf("All_B73_genes.pdf")
#par(mfrow=c(2,1))
#par(mar = rep(2,4))
#plot(x=NULL, y=NULL,xlim=c(-1100,1100),ylim=c(0,1),xlab="",ylab='read count',main='CG')
#lines(p.1.cg,col=1,lwd=1)
#lines(p.1.chg,col=2,lwd=1)
#lines(p.1.chh,col=3,lwd=1)
#xline(0,lty=2,col='black')
#xline(100,lty=2,col='black')
#legend("topright",c('CG','CHG', 'CHH'),col=c(1,2,3),lty=1,lwd=2,cex=0.7)
#dev.off()

#HERE's how to graph mC for single genes
data_table=fread("data.sub1.txt")

gene="Zm00001eb045540"
gene1=subset(data_table, Gene==gene)

#by using relative_distance, the gene body is normalized to positions 0-1000,
#1000-2000 is the 1kb downstream of the gene
ggplot(gene1, aes(x=relative_distance, y=cg, color="red"))+
  geom_point(alpha=0.5, color="red")+
  geom_line(color="red")+
  geom_vline(xintercept=0, linetype="dashed",
             color="black", alpha=0.5)+
  geom_vline(xintercept=1000, linetype="dashed",
             color="black", alpha=0.5)+ 
  ggtitle(paste(gene,"_cg",sep=""))

ggplot(gene1, aes(x=relative_distance, y=chg, color="blue"))+
  geom_point(alpha=0.5, color="blue")+
  geom_line(color="blue")+
  geom_vline(xintercept=0, linetype="dashed",
             color="black", alpha=0.5)+
  geom_vline(xintercept=1000, linetype="dashed",
             color="black", alpha=0.5)+ 
  ggtitle(paste(gene,"_chg",sep=""))

ggplot(gene1, aes(x=relative_distance, y=chh, color="black"))+
  geom_point(alpha=0.5, color="black")+
  geom_line(color="black")+
  geom_vline(xintercept=0, linetype="dashed",
             color="black", alpha=0.5)+
  geom_vline(xintercept=1000, linetype="dashed",
             color="black", alpha=0.5)+ 
  ggtitle(paste(gene,"_chh",sep=""))

`````
### NOTE - probably makes sense to manipulate the dataframe, plot CG, CHG, CHH on the same plot List of genes can be given to the subset command and output can be faceted by Gene rather than doing one by one

Here's what the plots look like for a the gene I chose at random:
<img width="835" alt="image" src="https://github.com/acread/PRFB_TEs_and_NLRs/assets/43852873/ef9cafc4-0715-45b1-aef2-234b8e6c900f">

This tracks pretty well with the genome browser snapshot - we just start to see a few bumps in the promoter \
and there are TEs just downstream of the gene that are heavily methylated.  The bumps in the gene body don't\
completely track, but these may drop way down when the other two replicates are included

![image](https://github.com/acread/PRFB_TEs_and_NLRs/assets/43852873/bff25d33-18bb-4264-9c85-0e5d3bda5a45)

### Metaplots of merged data from the 3 reps

I ran the maize data through Pete's pipeline to generate tile-based mC levels in each context \
Zach provided a categorization file for each gene based on whether it is DE in: no datasets (Exp not DE), \
DE in 1, DE in 5 or more (Robust), or DE in all 6 (Core) -- I'm doing various splits based on these data \
to see if there are observable differences in methylation patterns at or near genes.

<details>
  <summary>R code</summary>
        
`````R
library(ggplot2)
library(data.table)
library(tidyverse)
library(janitor)
library(fields)

setwd("/Users/read0094/Desktop/Maize/ForZach/")

data_table=fread("B73_merge_data.sub1.txt")


#to subset to Zach's subset lists
Intersect_Zach=fread("Intersections.Up.Zm01.02.03.04.05.06.B73.tsv")
Robust_DEGs=subset(Intersect_Zach, intersectionDegree>=5)
Robust_DEGs_subset=Robust_DEGs %>% dplyr::select(1)
Robust_DEGs_subset=dplyr::rename(Robust_DEGs_subset, Gene=gid)
Robust_DEGs_subset$status="Robust"
Robust_DEGs_subset_more=dplyr::left_join(Robust_DEGs_subset, data_table)

Core_DEGs=subset(Intersect_Zach, intersectionDegree==6)
Core_DEGs_subset=Core_DEGs %>% dplyr::select(1)
Core_DEGs_subset=dplyr::rename(Core_DEGs_subset, Gene=gid)
Core_DEGs_subset$status="Core"
Core_DEGs_subset_more=dplyr::left_join(Core_DEGs_subset, data_table)


DE_once=subset(Intersect_Zach, intersectionDegree==1)
DE_once_subset=DE_once %>% dplyr::select(1)
DE_once_subset=dplyr::rename(DE_once_subset, Gene=gid)
DE_once_subset$status="DE_once"
DE_once_subset_more=dplyr::left_join(DE_once_subset, data_table)

##
#stats.test1.cg=stats.bin(Robust_DEGs_subset_more$distance,Robust_DEGs_subset_more$cg,N=42)
stats.Robust.cg=stats.bin(Robust_DEGs_subset_more$relative_distance,Robust_DEGs_subset_more$cg,N=42)
#stats.test1.chg=stats.bin(Robust_DEGs_subset_more$distance,Robust_DEGs_subset_more$chg,N=42)
stats.Robust.chg=stats.bin(Robust_DEGs_subset_more$relative_distance,Robust_DEGs_subset_more$chg,N=42)
#stats.test1.chh=stats.bin(Robust_DEGs_subset_more$distance,Robust_DEGs_subset_more$chh,N=42)
stats.Robust.chh=stats.bin(Robust_DEGs_subset_more$relative_distance,Robust_DEGs_subset_more$chh,N=42)

p.Robust.cg=cbind(matrix(stats.Robust.cg$centers,ncol=1),stats.Robust.cg$stats["mean",])
p.Robust.chg=cbind(matrix(stats.Robust.cg$centers,ncol=1),stats.Robust.chg$stats["mean",])
p.Robust.chh=cbind(matrix(stats.Robust.cg$centers,ncol=1),stats.Robust.chh$stats["mean",])

##
plot(x=NULL, y=NULL,xlim=c(-1100,2100),ylim=c(0,1),xlab="",ylab='read count',main='Robust')
lines(p.Robust.cg,col=1,lwd=1)
lines(p.Robust.chg,col=2,lwd=1)
lines(p.Robust.chh,col=3,lwd=1)
#lines(p.7.cg,col=7,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

###
stats.Core.cg=stats.bin(Core_DEGs_subset_more$relative_distance,Core_DEGs_subset_more$cg,N=42)
#stats.test1.chg=stats.bin(Core_DEGs_subset_more$distance,Core_DEGs_subset_more$chg,N=42)
stats.Core.chg=stats.bin(Core_DEGs_subset_more$relative_distance,Core_DEGs_subset_more$chg,N=42)
#stats.test1.chh=stats.bin(Core_DEGs_subset_more$distance,Core_DEGs_subset_more$chh,N=42)
stats.Core.chh=stats.bin(Core_DEGs_subset_more$relative_distance,Core_DEGs_subset_more$chh,N=42)

p.Core.cg=cbind(matrix(stats.Core.cg$centers,ncol=1),stats.Core.cg$stats["mean",])
p.Core.chg=cbind(matrix(stats.Core.cg$centers,ncol=1),stats.Core.chg$stats["mean",])
p.Core.chh=cbind(matrix(stats.Core.cg$centers,ncol=1),stats.Core.chh$stats["mean",])


###
#stats.DE_once.cg=stats.bin(DE_once_more$distance,DE_once_more$cg,N=42)
stats.DE_once.cg=stats.bin(DE_once_subset_more$relative_distance,DE_once_subset_more$cg,N=42)
#stats.DE_once.chg=stats.bin(DE_once_more$distance,DE_once_more$chg,N=42)
stats.DE_once.chg=stats.bin(DE_once_subset_more$relative_distance,DE_once_subset_more$chg,N=42)
#stats.DE_once.chh=stats.bin(DE_once_more$distance,DE_once_more$chh,N=42)
stats.DE_once.chh=stats.bin(DE_once_subset_more$relative_distance,DE_once_subset_more$chh,N=42)

p.DE_once.cg=cbind(matrix(stats.DE_once.cg$centers,ncol=1),stats.DE_once.cg$stats["mean",])
p.DE_once.chg=cbind(matrix(stats.DE_once.cg$centers,ncol=1),stats.DE_once.chg$stats["mean",])
p.DE_once.chh=cbind(matrix(stats.DE_once.cg$centers,ncol=1),stats.DE_once.chh$stats["mean",])

plot(x=NULL, y=NULL,xlim=c(-1100,2100),ylim=c(0,1),xlab="",ylab='read count',main='Robust v DEonce')
lines(p.Robust.cg,col="darkblue",lwd=1)
lines(p.DE_once.cg,col="blue",lwd=1)
lines(p.Robust.chg,col="darkred",lwd=1)
lines(p.DE_once.chg,col="red",lwd=1)
lines(p.Robust.chh,col="darkgreen",lwd=1)
lines(p.DE_once.chh,col="green",lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')


###Bringing in the 'Exp_notDE' and 'All_Exp' lists
All_Expressed_DEGs=fread("ZmB73.AllExpressedGenes.tsv")
#All_Expressed_DEGs_subset=All_Expressed_DEGs %>% dplyr::select(1)
All_Expressed_DEGs=dplyr::rename(All_Expressed_DEGs, Gene=gid)
All_Expressed_DEGs$status="All_Expressed"
All_Expressed_DEGs_more=dplyr::left_join(All_Expressed_DEGs, data_table)

stats.All_Expressed.cg=stats.bin(All_Expressed_DEGs_more$relative_distance,All_Expressed_DEGs_more$cg,N=42)
#stats.test1.chg=stats.bin(All_Expressed_DEGs_subset_more$distance,All_Expressed_DEGs_subset_more$chg,N=42)
stats.All_Expressed.chg=stats.bin(All_Expressed_DEGs_more$relative_distance,All_Expressed_DEGs_more$chg,N=42)
#stats.test1.chh=stats.bin(All_Expressed_DEGs_subset_more$distance,All_Expressed_DEGs_subset_more$chh,N=42)
stats.All_Expressed.chh=stats.bin(All_Expressed_DEGs_more$relative_distance,All_Expressed_DEGs_more$chh,N=42)

p.All_Expressed.cg=cbind(matrix(stats.All_Expressed.cg$centers,ncol=1),stats.All_Expressed.cg$stats["mean",])
p.All_Expressed.chg=cbind(matrix(stats.All_Expressed.cg$centers,ncol=1),stats.All_Expressed.chg$stats["mean",])
p.All_Expressed.chh=cbind(matrix(stats.All_Expressed.cg$centers,ncol=1),stats.All_Expressed.chh$stats["mean",])

##
plot(x=NULL, y=NULL,xlim=c(-1100,2100),ylim=c(0,1),xlab="",ylab='read count',main='All_Expressed')
lines(p.All_Expressed.cg,col=1,lwd=1)
lines(p.All_Expressed.chg,col=2,lwd=1)
lines(p.All_Expressed.chh,col=3,lwd=1)
#lines(p.7.cg,col=7,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

#CorevALL
plot(x=NULL, y=NULL,xlim=c(-610,1500),ylim=c(0,0.75),xlab="",ylab='methylation',main='Core v AllExp',xaxt='n')
lines(p.Core.cg,col="blue",lwd=1)
lines(p.All_Expressed.cg,col="blue",lwd=1, lty='dashed')
lines(p.Core.chg,col="red",lwd=1)
lines(p.All_Expressed.chg,col="red",lwd=1, lty='dashed')
lines(p.Core.chh,col="darkgreen",lwd=1)
lines(p.All_Expressed.chh,col="darkgreen",lwd=1, lty='dashed')
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("bottomright", legend=c("Core CG", "All CG","Core CHG", "All CHG","Core CHH", "All CHH"),
       col=c("blue", "blue","red", "red","darkgreen","darkgreen"), lty=1:2, cex=0.5,box.lty=0)

plot(x=NULL, y=NULL,xlim=c(-1100,2100),ylim=c(0,0.1),xlab="",ylab='read count',main='Core v AllExp CHH')
lines(p.Core.chh,col="darkgreen",lwd=1)
lines(p.All_Expressed.chh,col="green",lwd=1, lty='dashed')
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')

#RobustvAll
plot(x=NULL, y=NULL,xlim=c(-610,1500),ylim=c(0,0.75),xlab="",ylab='methylation',main='Robust v AllExp', xaxt='n')
lines(p.Robust.cg,col="blue",lwd=1)
lines(p.All_Expressed.cg,col="blue",lwd=1, lty='dashed')
lines(p.Robust.chg,col="red",lwd=1)
lines(p.All_Expressed.chg,col="red",lwd=1, lty='dashed')
lines(p.Robust.chh,col="darkgreen",lwd=1)
lines(p.All_Expressed.chh,col="darkgreen",lwd=1, lty='dashed')
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("bottomright", legend=c("Robust CG", "All CG","Robust CHG", "All CHG","Robust CHH", "All CHH"),
       col=c("blue", "blue","red", "red","darkgreen","darkgreen"), lty=1:2, cex=0.5,box.lty=0)


plot(x=NULL, y=NULL,xlim=c(-1100,2100),ylim=c(0,0.1),xlab="",ylab='read count',main='Robust v AllExp CHH')
lines(p.Robust.chh,col="darkgreen",lwd=1)
lines(p.All_Expressed.chh,col="green",lwd=1, lty='dashed')
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
`````
</details>

Here are a few examples of metaplots -- it does look like there are different methylation signatures for the genes that \
are generally upregulated after heat stress, but it is unclear whether this is biologically relevant.

![image](https://github.com/acread/PRFB_TEs_and_NLRs/assets/43852873/65c796fb-3be6-42b1-9c94-a353c86a0ca3)

![image](https://github.com/acread/PRFB_TEs_and_NLRs/assets/43852873/5d10cd19-8ce0-4c39-b19f-7be40ab0565f)





