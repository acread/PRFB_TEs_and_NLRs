## Calculating Tau - a measure of plasticity/tissue specificity - for all genes in each NAM line



`````R
library(here)
library(arrow)
library(BioNERO)
library(tidyverse)
#library(ComplexHeatmap)*
#library(patchwork)*
#library(clusterProfiler)*
library(Biostrings)
#library(planttfhunter)*
  
  calculate_tau <- function(x) {
    
    tau <- NA
    if(all(!is.na(x)) & min(x, na.rm = TRUE) >= 0) {
      tau <- 0
      if(max(x) != 0) {
        x <- (1-(x / max(x, na.rm = TRUE)))
        tau <- sum(x, na.rm = TRUE)
        tau <- tau / (length(x) - 1)
      }
    }
    
    return(tau)
  }
  
  setwd("/Users/read0094/Desktop/Maize")
  
  #read in all the normalized expression
  All_Expression=fread("norm_NAM_RNAseq_counts.txt")
  Gene_Names=All_Expression %>% dplyr::select(1)
  Gene_Names=as.data.frame(Gene_Names)
  
  #set genotype
  #genotype="B97"
  #genotype="B73"
  #genotype="CML103"
  #genotype="CML228"
  genotype="CML247"
  #genotype="CML277"
  #genotype="CML322"
  #genotype="CML333"
  #genotype="CML52"
  #genotype="CML69"
  #genotype="HP301"
  #genotype="Il14H"
  #genotype="Ki11"
  #genotype="Ki3"
  #genotype="Ky21"
  #genotype="M162W"
  #genotype="M37W"
  #genotype="Mo18W"
  #genotype="Ms71"
  #genotype="NC350"
  #genotype="NC358"
  #genotype="Oh43"
  #genotype="Oh7B"
  #genotype="P39"
  #genotype="Tx303"
  #genotype="Tzi8"
  
  #Just B73
  Genotype_Expression=dplyr::select(All_Expression, matches(genotype))
  
  #This seems crazy but I'm extracting each tissue type and calculating a median
  Genotype_root=Genotype_Expression %>% dplyr::select(matches("root"))
  Genotype_root=Genotype_root %>% rowwise() %>%
    mutate(root_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_root=Genotype_root[,(ncol(Genotype_root))]
  
  Genotype_shoot=Genotype_Expression %>% dplyr::select(matches("shoot"))
  Genotype_shoot=Genotype_shoot %>% rowwise() %>%
    mutate(shoot_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_shoot=Genotype_shoot[,(ncol(Genotype_shoot))]
  
  Genotype_leaf_base=Genotype_Expression %>% dplyr::select(matches("leaf.base"))
  Genotype_leaf_base=Genotype_leaf_base %>% rowwise() %>%
    mutate(leaf_base_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_leaf_base=Genotype_leaf_base[,(ncol(Genotype_leaf_base))]
  
  Genotype_leaf=Genotype_Expression %>% dplyr::select(matches("leaf_"))
  Genotype_leaf=Genotype_leaf %>% rowwise() %>%
    mutate(leaf_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_leaf=Genotype_leaf[,(ncol(Genotype_leaf))]
  
  Genotype_leaf_tip=Genotype_Expression %>% dplyr::select(matches("leaf.tip"))
  Genotype_leaf_tip=Genotype_leaf_tip %>% rowwise() %>%
    mutate(leaf_tip_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_leaf_tip=Genotype_leaf_tip[,(ncol(Genotype_leaf_tip))]
  
  Genotype_tassel=Genotype_Expression %>% dplyr::select(matches("tassel"))
  Genotype_tassel=Genotype_tassel %>% rowwise() %>%
    mutate(tassel_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_tassel=Genotype_tassel[,(ncol(Genotype_tassel))]
  
  Genotype_ear=Genotype_Expression %>% dplyr::select(matches("ear"))
  Genotype_ear=Genotype_ear %>% rowwise() %>%
    mutate(ear_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_ear=Genotype_ear[,(ncol(Genotype_ear))]
  
  Genotype_anther=Genotype_Expression %>% dplyr::select(matches("anther"))
  Genotype_anther=Genotype_anther %>% rowwise() %>%
    mutate(anther_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_anther=Genotype_anther[,(ncol(Genotype_anther))]
  
  Genotype_endosperm=Genotype_Expression %>% dplyr::select(matches("endosperm"))
  Genotype_endosperm=Genotype_endosperm %>% rowwise() %>%
    mutate(endosperm_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_endosperm=Genotype_endosperm[,(ncol(Genotype_endosperm))]
  
  Genotype_embryo=Genotype_Expression %>% dplyr::select(matches("embryo"))
  Genotype_embryo=Genotype_embryo %>% rowwise() %>%
    mutate(embryo_med=median(c_across(where(is.numeric)),na.rm=TRUE))
  Genotype_embryo=Genotype_embryo[,(ncol(Genotype_embryo))]
  
  All_median=dplyr::bind_cols(Genotype_anther,
                              Genotype_ear,
                              Genotype_endosperm,
                              #Genotype_embryo,
                              Genotype_leaf,
                              Genotype_leaf_tip,
                              Genotype_leaf_base,
                              Genotype_root,
                              Genotype_shoot,
                              Genotype_tassel)
  All_median=as.data.frame(All_median)
  #All_median$gene_ID=Gene_Names[,1]
  rownames(All_median)=Gene_Names[,1]
  
  rm(list=ls()[grep("Genotype",ls())])
  
  #remove genes with median count <1 in all tissues
  Unexpressed_genes=as.data.frame(apply(All_median, 1 ,max ,na.rm=TRUE))
  
  Unexpressed_genes=subset(Unexpressed_genes, `apply(All_median, 1, max, na.rm = TRUE)`<1)
  
  #adding geneID for anti-join
  Unexpressed_genes$geneID=rownames(Unexpressed_genes)
  All_median$geneID=rownames(All_median)
  
  Expressed_median=dplyr::anti_join(All_median,Unexpressed_genes)
  #getrid of that geneID column
  Expressed_median=Expressed_median %>% dplyr::select(1:9)
  
  
  tau=apply(log2(Expressed_median+1),1,calculate_tau)
  
  genes_tau=data.frame(
    geneID=names(tau),
    Tau=as.numeric(tau)
  )
  
  #Add GeneID back to Expressed median...
  Expressed_median=dplyr::anti_join(All_median,Unexpressed_genes) 
  #This datframe is supposed be pivoted/melted with a new column called 'Median'
  classified_genes=dplyr::left_join(Expressed_median, genes_tau)
  classified_genes=classified_genes %>% mutate(Median_greater1 = rowSums(classified_genes > 1))
  classified_genes=classified_genes %>% mutate(Median_greater5 = rowSums(classified_genes > 5))
  classified_genes$Median_greater5=classified_genes$Median_greater5-1
  
  classified_genes=classified_genes %>% mutate(
    Classification=case_when(
      Median_greater1==0 ~"Null",
      Median_greater5==0 ~"Weak",
      Median_greater5>=1 & Tau <0.85 ~ "Broad",
      Median_greater5>=1 & Tau >0.85 ~ "Specific"
    )
  )
 
  write_tsv(classified_genes, paste(genotype,"_tau.txt",sep=""))
`````

Based on the out file from above I am looking at classification across the NAM lines to characterize genes as: \
Broad 22452 = Classified as Broad in all 26 lines \
Specific 81 = Classified as Specific in all 26 lines \
Weak  183 = Classified as NA or Weak across all 26 lines \
Variable 9716 = does not meet above criteria (in other words, the gene falls into multiple categories)

Of the Variable class - 1975 of the genes are only a combo of Broad and Specific \
A total of 3921 genes are Weak in 2 or fewer genotypes

