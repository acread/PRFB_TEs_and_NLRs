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



