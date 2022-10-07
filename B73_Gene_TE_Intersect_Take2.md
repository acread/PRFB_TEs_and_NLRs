
### Breakdown of Strategy: 
1. a curated list of B73 structural TEs was used (no helitrons) 
2. the canonical transcripts from B73 were used to generate gene and intron files 
3. I determined any TE that overlaps coding sequence and deleted it for downstream analysis 
4. Remaining TEs were intersected with Intron locations to generate a list of TEs that overlap Introns 
    Note that there are introns that overlap both coding and intron - these are classed as coding here 
5. Remaining TEs were seperately intersected with 2kb upstream (promoter) and 2kb down files 
    I chose 2kb this time based on a recent Brassica paper https://onlinelibrary.wiley.com/doi/10.1111/pbi.13807
    
    This Brassica paper has the following figure - I am recreating this for B73\

<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184729988-cfd0d635-2884-4654-b404-335fb267a11d.png?raw=true" alt="image"/>
</p>

### Modifications to Strategy:
1. I used Manisha's curated structural and homology TE list instead of Structural only (though I did a grouped and seperate analysis)
2. Canonical transcripts were still used to generate the Intron files
3. I did Intron specific overlaps, however for the specific below analysis I'm just saying "does a TE overlap a canonical gene"
4. I created files with 500bp blocks up and downstream of each gene - 4 blocks each, so looking at 2kb up and down \
Note that I did not correct for any genes that are within 4kb of eachother, so there may be some genomic regions that get counted more than once.

To generate the below figure I calculated how many genes contained at least one of each TE element (Struct/Hom) \
There was not size or count cut-offs here - so a gene is counted the same if it has one tiny element \
or several quite large elements

![image](https://user-images.githubusercontent.com/43852873/194641175-37930262-efc4-474f-afca-b0cca8f033a2.png)

### The take-home from this is that the differences between NLR and Non-NLR genes is relatively small, EXCEPT \
### that NLR promoters seem to be enriched for DTH (PIF-Harbingers) and NLR gene-bodies enriched in LTRs, particulary\
### Copia elements

## R Code
<details>
  </br>

  `````R
  library("optparse")
library("tidyverse")
library("stringr")
library("data.table")
library("GenomicFeatures")
library("bedtoolsr")
library("MetBrewer")
#library("patchwork")
library("ggpubr")

setwd("/Users/read0094/Desktop/Maize")
setwd("/Users/read0094/Desktop/Maize/Gene_NLR_TE_Analyses/")

#####Step 3 intersect intron/TE####
#requires Manisha's minimal structural TE files
file3="B73_Introns_822.txt"
file4="FunTest/screened_Zm-B73-REFERENCE-NAM-5.0.TE.bed"

#Canonical gene bed
B73_Canonical=fread("B73_Genes.bed", header=F, sep="\t")
B73_Canonical$V5="Gene"

Introns=fread(file3, header=T, sep="\t")
Introns$gene_id=Introns$tx_name
Introns$gene_id=gsub('_.*', '', Introns$gene_id)
Introns$Gene_Intron=paste(Introns$gene_id,Introns$intron_num,sep=":")
Introns=Introns %>% dplyr::select(1,2,3,10,5)
Introns$Intron_Len=Introns$end-Introns$start
#get rid of the absurdly long introns
Introns=subset(Introns, Intron_Len <= 10000)
IntronLen=Introns %>% dplyr::select(4,6)

TEs=fread(file4, header=T, sep="\t")
TEs=TEs %>% dplyr::select(1,2,3,10,9)
StructTEs=TEs %>% filter(!grepl('=homo', attributes))

#Interserct Intron and Genes vs TEs
ANY_TE_IN_INTRON=bedtoolsr::bt.intersect(Introns,TEs,wo=T,F=1.0)
ANY_TE_OVERLAP_GENE=bedtoolsr::bt.intersect(B73_Canonical,TEs,wo=T)

STRUCT_TE_IN_INTRON=bedtoolsr::bt.intersect(Introns, StructTEs,wo=T,F=1.0)
STRUCT_TE_OVERLAP_GENE=bedtoolsr::bt.intersect(B73_Canonical,StructTEs,wo=T)

##### To Get 500bp windows flanking genes ####
PosGenes=subset(join1, V3=="gene")
PosGenes=subset(PosGenes, V7=="+")
NegGenes=subset(join1, V3=="gene")
NegGenes=subset(NegGenes, V7=="-")

AllGenes_bed=subset(join1, V3=='gene')
AllGenes_bed=AllGenes_bed %>% dplyr::select(1,4,5,10)

write_tsv(AllGenes_bed,"B73_canonical_bed_fromR.txt")

#500bp up down
PosGenes500bpup=PosGenes
PosGenes500bpup$V5=PosGenes500bpup$V4-1
PosGenes500bpup$V4=PosGenes500bpup$V4-500

PosGenes500bpdown=PosGenes
PosGenes500bpdown$V4=PosGenes500bpdown$V5+1
PosGenes500bpdown$V5=PosGenes500bpdown$V4+500

NegGenes500bpup=NegGenes
NegGenes500bpup$V4=NegGenes500bpup$V5+1
NegGenes500bpup$V5=NegGenes500bpup$V5+500

NegGenes500bpdown=NegGenes
NegGenes500bpdown$V5=NegGenes500bpdown$V4-1
NegGenes500bpdown$V4=NegGenes500bpdown$V4-500

Allgenes500bpup=bind_rows(PosGenes500bpup,NegGenes500bpup)
All500bpup_bed=Allgenes500bpup %>% dplyr::select(1,4,5,10)
All500bpup_bed=All500bpup_bed %>% filter(!grepl('scaf', V1))
All500bpup_bed$Class="500bpup"


Allgenes500bpdown=bind_rows(PosGenes500bpdown,NegGenes500bpdown)
All500bpdown_bed=Allgenes500bpdown %>% dplyr::select(1,4,5,10)
All500bpdown_bed=All500bpdown_bed %>% filter(!grepl('scaf', V1))
All500bpdown_bed$Class="500bpdown"

write_tsv(All500bpup_bed,"B73_canonical_bed_500up_fromR.txt")
write_tsv(All500bpdown_bed,"B73_canonical_bed_500downfromR.txt")


#rm(NegGenes500bpup,NegGenes500bpdown,PosGenes500bpup,PosGenes500bpdown,Allgenes500bupup,Allgenes500bpdown)

#1000bp up down
PosGenes1000bpup=PosGenes
PosGenes1000bpup$V5=PosGenes1000bpup$V4-501
PosGenes1000bpup$V4=PosGenes1000bpup$V4-1000

PosGenes1000bpdown=PosGenes
PosGenes1000bpdown$V4=PosGenes1000bpdown$V5+501
PosGenes1000bpdown$V5=PosGenes1000bpdown$V4+500

NegGenes1000bpup=NegGenes
NegGenes1000bpup$V4=NegGenes1000bpup$V5+501
NegGenes1000bpup$V5=NegGenes1000bpup$V5+500

NegGenes1000bpdown=NegGenes
NegGenes1000bpdown$V5=NegGenes1000bpdown$V4-501
NegGenes1000bpdown$V4=NegGenes1000bpdown$V4-500

Allgenes1000bpup=bind_rows(PosGenes1000bpup,NegGenes1000bpup)
All1000bpup_bed=Allgenes1000bpup %>% dplyr::select(1,4,5,10)
All1000bpup_bed=All1000bpup_bed %>% filter(!grepl('scaf', V1))
All1000bpup_bed$Class="1000bpup"


Allgenes1000bpdown=bind_rows(PosGenes1000bpdown,NegGenes1000bpdown)
All1000bpdown_bed=Allgenes1000bpdown %>% dplyr::select(1,4,5,10)
All1000bpdown_bed=All1000bpdown_bed %>% filter(!grepl('scaf', V1))
All1000bpdown_bed$Class="1000bpdown"

write_tsv(All1000bpup_bed,"B73_canonical_bed_1000up_fromR.txt")
write_tsv(All1000bpdown_bed,"B73_canonical_bed_1000downfromR.txt")

rm(NegGenes500bpup,NegGenes500bpdown,PosGenes500bpup,PosGenes500bpdown,Allgenes500bupup,Allgenes500bpdown)

#rm(NegGenes1000bpup,NegGenes1000bpdown,PosGenes1000bpup,PosGenes1000bpdown,Allgenes500bupup,Allgenes1000bpdown)

#1500bp up down
PosGenes1500bpup=PosGenes
PosGenes1500bpup$V5=PosGenes1500bpup$V4-1001
PosGenes1500bpup$V4=PosGenes1500bpup$V5-500

PosGenes1500bpdown=PosGenes
PosGenes1500bpdown$V4=PosGenes1500bpdown$V5+1001
PosGenes1500bpdown$V5=PosGenes1500bpdown$V4+500

NegGenes1500bpup=NegGenes
NegGenes1500bpup$V4=NegGenes1500bpup$V5+1001
NegGenes1500bpup$V5=NegGenes1500bpup$V5+500

NegGenes1500bpdown=NegGenes
NegGenes1500bpdown$V5=NegGenes1500bpdown$V4-1001
NegGenes1500bpdown$V4=NegGenes1500bpdown$V4-500

Allgenes1500bpup=bind_rows(PosGenes1500bpup,NegGenes1500bpup)
All1500bpup_bed=Allgenes1500bpup %>% dplyr::select(1,4,5,10)
All1500bpup_bed=All1500bpup_bed %>% filter(!grepl('scaf', V1))
All1500bpup_bed$Class="1500bpup"


Allgenes1500bpdown=bind_rows(PosGenes1500bpdown,NegGenes1500bpdown)
All1500bpdown_bed=Allgenes1500bpdown %>% dplyr::select(1,4,5,10)
All1500bpdown_bed=All1500bpdown_bed %>% filter(!grepl('scaf', V1))
All1500bpdown_bed$Class="1500bpdown"

write_tsv(All1500bpup_bed,"B73_canonical_bed_1500up_fromR.txt")
write_tsv(All1500bpdown_bed,"B73_canonical_bed_1500downfromR.txt")

rm(NegGenes1000bpup,NegGenes1000bpdown,PosGenes1000bpup,PosGenes1000bpdown,Allgenes500bupup,Allgenes1000bpdown)

#2000bp up down
PosGenes2000bpup=PosGenes
PosGenes2000bpup$V5=PosGenes2000bpup$V4-1501
PosGenes2000bpup$V4=PosGenes2000bpup$V5-500

PosGenes2000bpdown=PosGenes
PosGenes2000bpdown$V4=PosGenes2000bpdown$V5+1501
PosGenes2000bpdown$V5=PosGenes2000bpdown$V4+500

NegGenes2000bpup=NegGenes
NegGenes2000bpup$V4=NegGenes2000bpup$V5+1501
NegGenes2000bpup$V5=NegGenes2000bpup$V4+500

NegGenes2000bpdown=NegGenes
NegGenes2000bpdown$V5=NegGenes2000bpdown$V4-1501
NegGenes2000bpdown$V4=NegGenes2000bpdown$V4-500

Allgenes2000bpup=bind_rows(PosGenes2000bpup,NegGenes2000bpup)
All2000bpup_bed=Allgenes2000bpup %>% dplyr::select(1,4,5,10)
All2000bpup_bed=All2000bpup_bed %>% filter(!grepl('scaf', V1))
All2000bpup_bed$Class="2000bpup"


Allgenes2000bpdown=bind_rows(PosGenes2000bpdown,NegGenes2000bpdown)
All2000bpdown_bed=Allgenes2000bpdown %>% dplyr::select(1,4,5,10)
All2000bpdown_bed=All2000bpdown_bed %>% filter(!grepl('scaf', V1))
All2000bpdown_bed$Class="2000bpdown"

write_tsv(All2000bpup_bed,"B73_canonical_bed_2000up_fromR.txt")
write_tsv(All2000bpdown_bed,"B73_canonical_bed_2000downfromR.txt")

rm(NegGenes2000bpup,NegGenes2000bpdown,PosGenes2000bpup,PosGenes2000bpdown,Allgenes2000bupup,Allgenes2000bpdown)
rm(NegGenes1500bpup,NegGenes1500bpdown,PosGenes1500bpup,PosGenes1500bpdown,Allgenes1500bupup,Allgenes1500bpdown)

#####
#All TEs
ANY_TE_500bp_UP=bedtoolsr::bt.intersect(All500bpup_bed,TEs,wo=T)
ANY_TE_1000bp_UP=bedtoolsr::bt.intersect(All1000bpup_bed,TEs,wo=T)
ANY_TE_1500bp_UP=bedtoolsr::bt.intersect(All1500bpup_bed,TEs,wo=T)
ANY_TE_2000bp_UP=bedtoolsr::bt.intersect(All2000bpup_bed,TEs,wo=T)

ANY_TE_500bp_DOWN=bedtoolsr::bt.intersect(All500bpdown_bed,TEs,wo=T)
ANY_TE_1000bp_DOWN=bedtoolsr::bt.intersect(All1000bpdown_bed,TEs,wo=T)
ANY_TE_1500bp_DOWN=bedtoolsr::bt.intersect(All1500bpdown_bed,TEs,wo=T)
ANY_TE_2000bp_DOWN=bedtoolsr::bt.intersect(All2000bpdown_bed,TEs,wo=T)

BIND_ANY_TEs=dplyr::bind_rows(ANY_TE_500bp_UP,ANY_TE_1000bp_UP,
                                 ANY_TE_1500bp_UP,ANY_TE_2000bp_UP,
                                 ANY_TE_500bp_DOWN, ANY_TE_1000bp_DOWN,
                                 ANY_TE_1500bp_DOWN, ANY_TE_2000bp_DOWN,
                                 ANY_TE_OVERLAP_GENE)

#StructTEs
STRUCT_TE_500bp_UP=bedtoolsr::bt.intersect(All500bpup_bed,StructTEs,wo=T)
STRUCT_TE_1000bp_UP=bedtoolsr::bt.intersect(All1000bpup_bed,StructTEs,wo=T)
STRUCT_TE_1500bp_UP=bedtoolsr::bt.intersect(All1500bpup_bed,StructTEs,wo=T)
STRUCT_TE_2000bp_UP=bedtoolsr::bt.intersect(All2000bpup_bed,StructTEs,wo=T)

STRUCT_TE_500bp_DOWN=bedtoolsr::bt.intersect(All500bpdown_bed,StructTEs,wo=T)
STRUCT_TE_1000bp_DOWN=bedtoolsr::bt.intersect(All1000bpdown_bed,StructTEs,wo=T)
STRUCT_TE_1500bp_DOWN=bedtoolsr::bt.intersect(All1500bpdown_bed,StructTEs,wo=T)
STRUCT_TE_2000bp_DOWN=bedtoolsr::bt.intersect(All2000bpdown_bed,StructTEs,wo=T)

BIND_STRUCT_TEs=dplyr::bind_rows(STRUCT_TE_500bp_UP,STRUCT_TE_1000bp_UP,
                       STRUCT_TE_1500bp_UP,STRUCT_TE_2000bp_UP,
                       STRUCT_TE_500bp_DOWN, STRUCT_TE_1000bp_DOWN,
                       STRUCT_TE_1500bp_DOWN, STRUCT_TE_2000bp_DOWN,
                       STRUCT_TE_OVERLAP_GENE)

#Add NLR data
B73_NLRs=fread("/Users/read0094/Desktop/Maize/Gene_NLR_TE_Analyses/B73_NLRs_2022.txt",
               header=F)
B73_NLRs$NLR="NLR"
B73_NLRs=B73_NLRs %>% dplyr::rename(V4=V1)

#Combine ANY TEs
BIND_ANY_TEs1=left_join(BIND_ANY_TEs,B73_NLRs)
rm(BIND_ANY_TEs)
BIND_ANY_TEs1$NLR=BIND_ANY_TEs1$NLR %>% replace_na("NO")
BIND_ANY_TEs1$TEtype=BIND_ANY_TEs1$V9
BIND_ANY_TEs1$TEtype=gsub("/DNA","",BIND_ANY_TEs1$TEtype)
BIND_ANY_TEs1$TEtype=gsub("/MITE","",BIND_ANY_TEs1$TEtype)
BIND_ANY_TEs1$TEtype=gsub("/L1","",BIND_ANY_TEs1$TEtype)
BIND_ANY_TEs1$TEtype=gsub("/RTE","",BIND_ANY_TEs1$TEtype)
BIND_ANY_TEs1$TEtype=gsub("LINE/unknown","LINE",BIND_ANY_TEs1$TEtype)

level_order=c("2000bpup","1500bpup","1000bpup","500bpup",
              "Gene",
              "500bpdown","1000bpdown","1500bpdown","2000bpdown")

#I think I should do counts and output into a new, small dataframe

###The Manisha approach
#Manisha helped me do some grouping and summarizing 
ANYx=BIND_ANY_TEs1 %>% group_by(V5,TEtype,NLR) %>% summarize(n=n()) %>%
  ungroup() %>% mutate(prop=case_when(
    NLR=='NLR' ~ n/176*100,
    TRUE ~ n/39035*100
  ))
#This gives the proportion of all genes that contain each type
#of TE (ANY only) vs the number of NLRs with each TE type -- note that TE 
#number is much much lower
BIND_ANY_TEs1 %>% group_by(V5,TEtype,NLR) %>% summarize(n=n()) %>%
  ungroup() %>% mutate(prop=case_when(
    NLR=='NLR' ~ n/176*100,
    TRUE ~ n/39035*100
  )) %>% mutate(level_order=factor(V5,levels=level_order)) %>% 
  ggplot(aes(x=level_order,y=prop,fill=NLR))+
  geom_col(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=met.brewer("Egypt"))+
  theme_bw()+
  facet_wrap(~TEtype)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Proportion of B73 genes that encode Hom/Struct TE elements")



#Combine STRUCTURAL TEs
BIND_STRUCT_TEs1=left_join(BIND_STRUCT_TEs,B73_NLRs)
rm(BIND_STRUCT_TEs)
BIND_STRUCT_TEs1$NLR=BIND_STRUCT_TEs1$NLR %>% replace_na("NO")
BIND_STRUCT_TEs1$TEtype=BIND_STRUCT_TEs1$V9
BIND_STRUCT_TEs1$TEtype=gsub("/DNA","",BIND_STRUCT_TEs1$TEtype)
BIND_STRUCT_TEs1$TEtype=gsub("/MITE","",BIND_STRUCT_TEs1$TEtype)

level_order=c("2000bpup","1500bpup","1000bpup","500bpup",
              "Gene",
              "500bpdown","1000bpdown","1500bpdown","2000bpdown")

#I think I should do counts and output into a new, small dataframe

###The Manisha approach
#Manisha helped me do some grouping and summarizing 
x=BIND_STRUCT_TEs1 %>% group_by(V5,TEtype,NLR) %>% summarize(n=n()) %>%
  ungroup() %>% mutate(prop=case_when(
    NLR=='NLR' ~ n/176*100,
    TRUE ~ n/39035*100
  ))
#This gives the proportion of all genes that contain each type
#of TE (STruct only) vs the number of NLRs with each TE type -- note that TE 
#number is much much lower
BIND_STRUCT_TEs1 %>% group_by(V5,TEtype,NLR) %>% summarize(n=n()) %>%
ungroup() %>% mutate(prop=case_when(
    NLR=='NLR' ~ n/176*100,
    TRUE ~ n/39035*100
  )) %>% mutate(level_order=factor(V5,levels=level_order)) %>% 
  ggplot(aes(x=level_order,y=prop,fill=NLR))+
  geom_col(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=met.brewer("Egypt"))+
  theme_bw()+
  facet_wrap(~TEtype)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Proportion of B73 genes that encode Structural TE elements")


### Create a table for each TE type

ANY_DTA=subset(ANY_TE_OVERLAP_GENE, V9=="DTA/DNA" | V9=="DTA/MITE")
ANY_DTC=subset(ANY_TE_OVERLAP_GENE, V9=="DTC/DNA" | V9=="DTC/MITE")
ANY_DTH=subset(ANY_TE_OVERLAP_GENE, V9=="DTH/DNA" | V9=="DTH/MITE")
ANY_DTT=subset(ANY_TE_OVERLAP_GENE, V9=="DTT/DNA" | V9=="DTT/MITE")
ANY_LINE=subset(ANY_TE_OVERLAP_GENE, V9=="LINE/unknown" | V9=="LINE/RTE" | V9=="LINE/L1")
ANY_Copia=subset(ANY_TE_OVERLAP_GENE, V9=="LTR/Copia")
ANY_Ty3=subset(ANY_TE_OVERLAP_GENE, V9=="LTR/Ty3")
ANY_LTRunknown=subset(ANY_TE_OVERLAP_GENE, V9=="LTR/unknown")

#write_tsv(ANY_Copia,"B73_Genic_Copia_LTRs.txt")

ANY_DTA_count <- ANY_DTA %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTA_count = ANY_DTA_count %>% dplyr::rename("DTAcount"="count")
ANY_DTA_countNLR=left_join(ANY_DTA_count,B73_NLRs)

ANY_DTC_count <- ANY_DTC %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTC_count = ANY_DTC_count %>% dplyr::rename("DTCcount"="count")
ANY_DTC_countNLR=left_join(ANY_DTC_count,B73_NLRs)

ANY_DTT_count <- ANY_DTT %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTT_count = ANY_DTT_count %>% dplyr::rename("DTTcount"="count")
ANY_DTT_countNLR=left_join(ANY_DTT_count,B73_NLRs)

ANY_DTH_count <- ANY_DTH %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTH_count = ANY_DTH_count %>% dplyr::rename("DTHcount"="count")
ANY_DTH_countNLR=left_join(ANY_DTH_count,B73_NLRs)

ANY_LINE_count <- ANY_LINE %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_LINE_count = ANY_LINE_count %>% dplyr::rename("LINEcount"="count")
ANY_LINE_countNLR=left_join(ANY_LINE_count,B73_NLRs)

ANY_Copia_count <- ANY_Copia %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_Copia_count = ANY_Copia_count %>% dplyr::rename("Copiacount"="count")
ANY_Copia_countNLR=left_join(ANY_Copia_count,B73_NLRs)

ANY_Ty3_count <- ANY_Ty3 %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_Ty3_count = ANY_Ty3_count %>% dplyr::rename("Ty3count"="count")
ANY_Ty3_countNLR=left_join(ANY_Ty3_count,B73_NLRs)

ANY_LTRunknown_count <- ANY_LTRunknown %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_LTRunknown_count = ANY_LTRunknown_count %>% dplyr::rename("LTRunknowncount"="count")
ANY_LTRunknown_countNLR=left_join(ANY_LTRunknown_count,B73_NLRs)

plot=fread("GeneswTEcount_HS.txt", header=T)

ggplot(plot, aes(x=V1,y=percent,fill=V4))+
  geom_col(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=met.brewer("Egypt"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Proportion of B73 genes that encode Struct/Hom TE elements")

##### Repeat for flanking regions #####
ANY_DTA_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="DTA/DNA" | V9=="DTA/MITE")
ANY_DTC_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="DTC/DNA" | V9=="DTC/MITE")
ANY_DTH_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="DTH/DNA" | V9=="DTH/MITE")
ANY_DTT_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="DTT/DNA" | V9=="DTT/MITE")
ANY_LINE_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="LINE/unknown" | V9=="LINE/RTE" | V9=="LINE/L1")
ANY_Copia_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="LTR/Copia")
ANY_Ty3_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="LTR/Ty3")
ANY_LTRunknown_1500bp_DOWN=subset(ANY_TE_1500bp_DOWN, V9=="LTR/unknown")

#write_tsv(ANY_Copia,"B73_Genic_Copia_LTRs.txt")

ANY_DTA_1500bp_DOWN_count <- ANY_DTA_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTA_1500bp_DOWN_count = ANY_DTA_1500bp_DOWN_count %>% dplyr::rename("DTAcount"="count")
ANY_DTA_1500bp_DOWN_countNLR=left_join(ANY_DTA_1500bp_DOWN_count,B73_NLRs)

ANY_DTC_1500bp_DOWN_count <- ANY_DTC_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTC_1500bp_DOWN_count = ANY_DTC_1500bp_DOWN_count %>% dplyr::rename("DTCcount"="count")
ANY_DTC_1500bp_DOWN_countNLR=left_join(ANY_DTC_1500bp_DOWN_count,B73_NLRs)

ANY_DTT_1500bp_DOWN_count <- ANY_DTT_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTT_1500bp_DOWN_count = ANY_DTT_1500bp_DOWN_count %>% dplyr::rename("DTTcount"="count")
ANY_DTT_1500bp_DOWN_countNLR=left_join(ANY_DTT_1500bp_DOWN_count,B73_NLRs)

ANY_DTH_1500bp_DOWN_count <- ANY_DTH_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_DTH_1500bp_DOWN_count = ANY_DTH_1500bp_DOWN_count %>% dplyr::rename("DTHcount"="count")
ANY_DTH_1500bp_DOWN_countNLR=left_join(ANY_DTH_1500bp_DOWN_count,B73_NLRs)

ANY_LINE_1500bp_DOWN_count <- ANY_LINE_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_LINE_1500bp_DOWN_count = ANY_LINE_1500bp_DOWN_count %>% dplyr::rename("LINEcount"="count")
ANY_LINE_1500bp_DOWN_countNLR=left_join(ANY_LINE_1500bp_DOWN_count,B73_NLRs)

ANY_Copia_1500bp_DOWN_count <- ANY_Copia_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_Copia_1500bp_DOWN_count = ANY_Copia_1500bp_DOWN_count %>% dplyr::rename("Copiacount"="count")
ANY_Copia_1500bp_DOWN_countNLR=left_join(ANY_Copia_1500bp_DOWN_count,B73_NLRs)

ANY_Ty3_1500bp_DOWN_count <- ANY_Ty3_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_Ty3_1500bp_DOWN_count = ANY_Ty3_1500bp_DOWN_count %>% dplyr::rename("Ty3count"="count")
ANY_Ty3_1500bp_DOWN_countNLR=left_join(ANY_Ty3_1500bp_DOWN_count,B73_NLRs)

ANY_LTRunknown_1500bp_DOWN_count <- ANY_LTRunknown_1500bp_DOWN %>%
  group_by(V4) %>%
  summarize(count = n()) %>%
  ungroup()
ANY_LTRunknown_1500bp_DOWN_count = ANY_LTRunknown_1500bp_DOWN_count %>% dplyr::rename("LTRunknowncount"="count")
ANY_LTRunknown_1500bp_DOWN_countNLR=left_join(ANY_LTRunknown_1500bp_DOWN_count,B73_NLRs)

nrow(ANY_DTA_1500bp_DOWN_countNLR)
nrow(ANY_DTC_1500bp_DOWN_countNLR)
nrow(ANY_DTH_1500bp_DOWN_countNLR)
nrow(ANY_DTT_1500bp_DOWN_countNLR)
nrow(ANY_LINE_1500bp_DOWN_countNLR)
nrow(ANY_LTRunknown_1500bp_DOWN_countNLR)
nrow(ANY_Ty3_1500bp_DOWN_countNLR)
nrow(ANY_Copia_1500bp_DOWN_countNLR)


table(ANY_DTA_1500bp_DOWN_countNLR$NLR)
table(ANY_DTC_1500bp_DOWN_countNLR$NLR)
table(ANY_DTH_1500bp_DOWN_countNLR$NLR)
table(ANY_DTT_1500bp_DOWN_countNLR$NLR)
table(ANY_LINE_1500bp_DOWN_countNLR$NLR)
table(ANY_LTRunknown_1500bp_DOWN_countNLR$NLR)
table(ANY_Ty3_1500bp_DOWN_countNLR$NLR)
table(ANY_Copia_1500bp_DOWN_countNLR$NLR)

#GeneswTEcount_inandNear_HS.txt

#####
plot=fread("GeneswTEcount_inandNear_HS.txt", header=T)
plot$V5=factor(plot$V5, levels=c("2kbUP","1.5kbUP","1kbUP","500bpUP",
                                 "GENE","500bpDOWN","1kbDOWN","1.5kbDOWN","2kbDOWN"))

#level_order1=c("2kbUP","1.5kbUP","1kbUP","500bpUP",
#              "GENE",
#              "500bpDOWN","1kbDOWN","1.5kbDOWN","2kbDOWN")

#probably a good idea to set the order of teh facet wrap
ggplot(plot, aes(x=V1,y=percent,fill=V4))+
  geom_col(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=met.brewer("Egypt"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size=4))+
  ggtitle("Proportion of B73 genes that encode Struct/Hom TE elements")+
  facet_wrap(~V5, nrow=1)
  `````

</details>

### Here is the B73 data


### ALL GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731583-ab2c567e-7b18-4ba6-a4eb-60b6953b7f9b.png?raw=true" alt="image"/>
</p>

### NLR GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731759-578c3e9c-f3d4-448f-a730-2b7cb427bfcd.png?raw=true" alt="image"/>
</p>


## Reorganizing the data
This is interesting but it's challenging because the datasets are such different sizes - Manisha suggested pulling 2 more subsets of genes \
and seeing what their patterns look like

<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/185230446-e4eaadb2-e142-4369-9f0a-9cd997697f2c.png?raw=true" alt="image"/>
</p>


