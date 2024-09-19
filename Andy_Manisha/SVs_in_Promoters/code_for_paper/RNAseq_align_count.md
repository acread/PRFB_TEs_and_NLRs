### Re-Indexing the B73v5 genome
I'm not sure why, but I seemed to be having some compatibility issues - so I repeated this with slight modifications \
and it seems to have worked - main difference is I specified the version of STAR
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

### Workflow to use a job array and submit all the NAM RNA-seq -- trying with just the first two files first to make sure things run
Here we're downloading the data from ENA, trimming, aligning to the genome with STAR, and counting using HT-SEQ - Note that HT-SEQ is \
aligning to "genes" here -- I'll want to repeat with "exons" to get exon specific numbers.
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
