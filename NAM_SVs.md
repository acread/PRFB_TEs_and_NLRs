




For each gene in B73, does it have an SV in or near (+/- 1kb) the gene?, if yes, parse out those that have only one SV in either B73 or \
B97, and this SV overlaps the promoter, and only the promoter. Some additional filters to get rid of genes that are too similar (no SVs) \
these were kept as a control set, and filter to get rid of genes that are too different based on exon-only comparison, and kick out genes \
that overlap 'unalignable regions'
***Note: Go back and create a category of multi-SV genes***

This leaves us with: \
       739 B73 genes with deletions in promoters relative to B97 \
       1256 B73 genes with insertions in promoters relative to B97 \
  and \
       11761 B73 genes that have no SVs (promoter or elsewhere)

The first thing I've done is look at the distribution of expression of these genes. \
I am pulling B73 expression from an existing TPM table (B73_full-tpm.tsv), only using one set of leaf tissues for initial analysis \
"B73_V11_middle".  ***Question: were these re-run for the NAM paper? it would be nice to have similar datasets*** \
For B97, ZM helped me map the B97 data back to B73, and again I'm pulling a leaf treatment only "leaf" - 4 reps \
I calculated the mean TPM for each gene and visualized with a violin plot -- note that the distribution of TPMs includes some very \
large values that I am not including here.  I also have not calculated a ratio of B73/B97 yet (though this is easy and on my list) \

### Here I am subsetting to only include genes with TPM between 0 and 500 - it looks like the "No SV" genes are biased toward lower expression, but I need to watch out because the n is very different (see above)

![image](https://user-images.githubusercontent.com/43852873/212158984-42b2362e-2a34-4b53-8d2f-2a2fb102c8c0.png)


## B73 leaf expression vs. B97 leaf expression 
Here I'm making a scatterplot of TPM -- I used a cut-off of 500TPM in either genotype.  I think I should try a ratio instead of the raw TPM
````R
B73_v_B97_SV_noSV_Expression_sub500=subset(B73_v_B97_SV_noSV_Expression, mean_leaf_B73 < 500 & mean_leaf_B97 < 500)
#Scatter plot of B73 v B97 expression colored by insertion type
ggplot(B73_v_B97_SV_noSV_Expression_sub500, aes(x=mean_leaf_B73, y=mean_leaf_B97, color=V21))+
  geom_point()+
  #scale_fill_manual(values=met.brewer("Kandinsky"))+
  geom_smooth(method=lm)+
  facet_wrap(~ V21)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
````

![image](https://user-images.githubusercontent.com/43852873/212142212-56ddd3f3-99a4-4cdc-9bc3-90c1d29362da.png)

### Here I'm showing fold change - where negative values are higher exp in B97
````
B73_v_B97_SV_noSV_Expression$ratio=(B73_v_B97_SV_noSV_Expression$mean_leaf_B73-B73_v_B97_SV_noSV_Expression$mean_leaf_B97)/B73_v_B97_SV_noSV_Expression$mean_leaf_B73
B73_v_B97_SV_noSV_Expression_B73ins=subset(B73_v_B97_SV_noSV_Expression, V21=="B73_Ins_Rel_B97")
B73_v_B97_SV_noSV_Expression_B73del=subset(B73_v_B97_SV_noSV_Expression, V21=="B73_Del_Rel_B97")
B73_v_B97_SV_noSV_Expression_noSV=subset(B73_v_B97_SV_noSV_Expression, V21=="No_SV_gene")
#subset the large noSV dataset
B73_v_B97_SV_noSV_Expression_noSV_subset=B73_v_B97_SV_noSV_Expression_noSV[sample(nrow(B73_v_B97_SV_noSV_Expression_noSV),1500),]


a=ggplot(B73_v_B97_SV_noSV_Expression_B73ins, aes(x=V6, y=ratio))+
  geom_point(color="green")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  ylim(-10,1)
b=ggplot(B73_v_B97_SV_noSV_Expression_B73del, aes(x=V6, y=ratio))+
  geom_point(color="red")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  ylim(-10,1)
c=ggplot(B73_v_B97_SV_noSV_Expression_noSV_subset, aes(x=V6, y=ratio))+
  geom_point(color="blue")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  ylim(-10,1)

grid.arrange(a,b,c)
````

![image](https://user-images.githubusercontent.com/43852873/212154765-3127aa07-59cf-4dd5-9702-1011d23bb6d9.png)

### Density plots

### Consider doing max expression-min/min -- I lose directionality but everything goes the same way

## Next Steps
Manisha is copying the B73 v B97 anchorwave prior to her cleanup into scratch - this allows me to do a bedtools intersect with gene +/- 1kb \
vs the anchorwave file to count SNPs and potentially small indels in my Prom_SV and No_SV gene sets

Look at B73vB97 TPM (x vs y) - color by insertion type/no sv, see if No SV genes have a higher Rsquared

Continue to get different tissue expression data

Begin to parse by: distance to TSS, Manisha category, etc

Are our promoter SV genes enriched for a specific type of gene(s)?
