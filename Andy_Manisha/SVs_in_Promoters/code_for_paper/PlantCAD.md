### PlantCAD to get gene scores for SV and no-SV genes for project with Manisha.
 
We have previously categorized each of the B73 genes compared to the equivalent genomic region in each other NAM line as either: 


     1. Containing an SV in the promoter, intron, or exon 
            -OR- 
     2. No SV in promoter, intron, or exon 
     
For genes that do have SVs, Manisha has grouped them based on the size and location of the SV.  So groups of NAM lines with similar patterns are grouped into haplogroups for each gene.
Importantly, this is not based on the annotation of individual NAM lines – it is only based on the annotation in B73 and the alignment of genomic regions with anchorwave.  In other words, it is agnostic to the annotation of the additional NAM lines.
In addition to our other analyses, we are interested to see if the genes we’ve identified with expression differences associated with SVs have distinct scores when measured with PlantCAD – a machine learning based scoring system based on the ‘rules’ for ‘real’ (non-psuedogenized) genes.  Note that this DOES require using the pan-gene ID to associate genes with each other (I believe). 
For each gene, Manisha has grouped them into haplogroups based on SVs and pulled the nucleotide sequence based on pan-gene ID.

### Here are the files I was provided to get started:

### gene_fastas 
A folder of fastas from 11491 B73 gene IDs. \
Each file includes the sequences of the gene from each of the NAM lines and each gene is put into a haplogroup.  This information is in the header for each entry –  Note that there are not 26 entries in each fasta because there are some genotypes that were filtered due to too many/few SVs \
examples: \
            >Zm00001eb188800:B73:H1 \
            >Zm00001eb188800:CML322:H3 \
            .. \
            >Zm00001eb188800:Oh7B:H1 
            
### 10_logit2zeroShot.R 
an R script required to scale the zeroshot scores 

### 2_HF_Logit_All_Tokens.py
This is the actual script that generates scores

Also a few template files for making SLURM scripts that I modified – I’m including some details of the modified versions 

Aimee uses a neat trick to run massive numbers of scripts on SciNet – this requires making a .cmds file that is basically a huge list of all the scripts you want to run.  In our case it’s the same command run for each fasta file
To run anything I need to set up a conda environment – initially I was able to point to Ana Berthel’s files, but there was maintenance on Atlas and the link broke.  I copied the environment to my 90day location – pytorch_pc

### 01_plantCAD.cmds
here’s the first line:
`````
python /project/psru_wheat_and_oat_genomics/plantCAD/2_HF_Logit_All_Tokens.py -input /project/psru_wheat_and_oat_genomics/plantCAD/gene_fastas/Zm00001eb000150_CDS.fa -outLogit /90daydata/psru_wheat_and_oat_genomics/plantCAD_output/Zm00001eb000150_CDS.npz  -model kuleshov-group/PlantCaduceus_l20 -numWorkers 0
plantCAD_AR.sh – this is the submission script for SciNet, note that it loads a few modules including parallel
`````
SLURM script
`````
#!/bin/bash -l
#SBATCH --time=60:00:00 # walltime limit (HH:MM:SS)
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=1 # 16 processor core(s) per node X 2 threads per core
#SBATCH --mem=50G # maximum memory per node
#SBATCH --partition=gpu-a100-mig7 # standard node(s)
#SBATCH --gres=gpu:a100_1g.10gb:1
#SBATCH --job-name="plantCAD_AR"
#SBATCH --mail-user=read0094@umn.edu # email address
#SBATCH --mail-type=ALL
#SBATCH -A psru_wheat_and_oat_genomics
 
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#module load miniconda/4.12.0
module load miniconda3
ml parallel
eval "$(conda shell.bash hook)"
 
source activate /90daydata/shared/ana.berthel/pytorch_pc
 
parallel -j 2 --joblog input${SLURM_ARRAY_TASK_ID}_joblog.txt --workdir $PWD < /90daydata/psru_wheat_and_oat_genomics/01_plantCAD.cmds
`````

The output of this script is stored in a new folder ‘plantCAD_output’ – the files are binary .npz files.  These are used as input into the R script to generate a corrected zeroShot score.

### 02_zeroShot.cmds
another long list of scripts – repeated for each gene ID

`````
Rscript /project/psru_wheat_and_oat_genomics/plantCAD/10_logit2zeroShot.R /90daydata/psru_wheat_and_oat_genomics/plantCAD_output/Zm00001eb378700_CDS.npz /project/psru_wheat_and_oat_genomics/plantCAD/gene_fastas/Zm00001eb378700_CDS.fa /90daydata/psru_wheat_and_oat_genomics/zeroShot_output/Zm00001eb378700_CDS.txt
`````

Ideally this is run as another Parallel array, however we couldn’t get it to run so I ran it as a series of array jobs.  Here’s an example for running from entry 4001 to 6000
A few notes:
1. SciNet does not like array numbers over 10000, it will fail.  I had to break out the final entries into a separate array.
2. Some of the input files are very large and they did not complete within the time given, I’ve had to go back and rescue these by submitting them with longer times

### 02_array_4k_6k.sh

`````
#!/bin/bash
#SBATCH --job-name=zeroshot
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb
#SBATCH --array=4001-6000
#SBATCH --nodes=1
#SBATCH --account=psru_wheat_and_oat_genomics
#SBATCH --mail-user=read0094@umn.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=array_slurm.out
 
ml miniconda3
source activate mafft_env
 
command=$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" /90daydata/psru_wheat_and_oat_genomics/02_zeroShot.cmds)
 
eval "$command"
`````
 
Output from this script are in zeroShot_output – the result is a .txt file for each geneID that includes the zeroShotScore and scaledzeroShotScore for each gene_NAM line – I’m not sure why column one is duplicated…  But we want either col1 or 2 combined with col 4 for all genes in a giant dataframe
`````
seqID   zeroShotScore   scaledZeroShotScore
Zm00001eb217320_CDS:Zm00001eb217320:B73:H1      Zm00001eb217320_CDS:Zm00001eb217320:B73:H1      2.85492179026573        -0.980242096259087
Zm00001eb217320_CDS:Zm00001eb217320:B97:H1      Zm00001eb217320_CDS:Zm00001eb217320:B97:H1      2.84504283317516        -1.99957264672426
Zm00001eb217320_CDS:Zm00001eb217320:CML103:H2   Zm00001eb217320_CDS:Zm00001eb217320:CML103:H2   2.86177092077827        -0.273535111750739
Zm00001eb217320_CDS:Zm00001eb217320:CML228:H2   Zm00001eb217320_CDS:Zm00001eb217320:CML228:H2   2.87064749989836        0.642368074970571
Zm00001eb217320_CDS:Zm00001eb217320:CML247:H1   Zm00001eb217320_CDS:Zm00001eb217320:CML247:H1   2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:CML277:H2   Zm00001eb217320_CDS:Zm00001eb217320:CML277:H2   2.86831217730242        0.401404818888352
Zm00001eb217320_CDS:Zm00001eb217320:CML322:H1   Zm00001eb217320_CDS:Zm00001eb217320:CML322:H1   2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:CML333:H2   Zm00001eb217320_CDS:Zm00001eb217320:CML333:H2   2.87024814719807        0.601162064650662
Zm00001eb217320_CDS:Zm00001eb217320:CML69:H2    Zm00001eb217320_CDS:Zm00001eb217320:CML69:H2    2.87639840912633        1.23575839089934
Zm00001eb217320_CDS:Zm00001eb217320:Il14H:H2    Zm00001eb217320_CDS:Zm00001eb217320:Il14H:H2    2.86177092077827        -0.273535111750739
Zm00001eb217320_CDS:Zm00001eb217320:Ki11:H1     Zm00001eb217320_CDS:Zm00001eb217320:Ki11:H1     2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:Ki3:H2      Zm00001eb217320_CDS:Zm00001eb217320:Ki3:H2      2.8706809756143 0.645822166292068
Zm00001eb217320_CDS:Zm00001eb217320:Ky21:H2     Zm00001eb217320_CDS:Zm00001eb217320:Ky21:H2     2.87164322521412        0.745109004495154
Zm00001eb217320_CDS:Zm00001eb217320:M162W:H2    Zm00001eb217320_CDS:Zm00001eb217320:M162W:H2    2.86982731020978        0.557739262527487
Zm00001eb217320_CDS:Zm00001eb217320:M37W:H2     Zm00001eb217320_CDS:Zm00001eb217320:M37W:H2     2.87064749989836        0.642368074970571
Zm00001eb217320_CDS:Zm00001eb217320:Mo18W:H1    Zm00001eb217320_CDS:Zm00001eb217320:Mo18W:H1    2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:Ms71:H2     Zm00001eb217320_CDS:Zm00001eb217320:Ms71:H2     2.87024814719807        0.601162064650662
Zm00001eb217320_CDS:Zm00001eb217320:NC350:H1    Zm00001eb217320_CDS:Zm00001eb217320:NC350:H1    2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:NC358:H1    Zm00001eb217320_CDS:Zm00001eb217320:NC358:H1    2.85344577218343        -1.13254059401838
Zm00001eb217320_CDS:Zm00001eb217320:Oh43:H2     Zm00001eb217320_CDS:Zm00001eb217320:Oh43:H2     2.87164322521412        0.745109004495154
Zm00001eb217320_CDS:Zm00001eb217320:Oh7B:H2     Zm00001eb217320_CDS:Zm00001eb217320:Oh7B:H2     2.87442956214766        1.0326088219569
Zm00001eb217320_CDS:Zm00001eb217320:P39:H2      Zm00001eb217320_CDS:Zm00001eb217320:P39:H2      2.87639840912633        1.23575839089934
Zm00001eb217320_CDS:Zm00001eb217320:Tzi8:H2     Zm00001eb217320_CDS:Zm00001eb217320:Tzi8:H2     2.87639840912633        1.23575839089934
`````
 
 
 

