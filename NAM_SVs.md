




For each gene in B73, does it have an SV in or near (+/- 1kb) the gene?, if yes, parse out those that have only one SV in either B73 or \
B97, and this SV overlaps the promoter, and only the promoter. Some additional filters to get rid of genes that are too similar (no SVs) \
these were kept as a control set, and filter to get rid of genes that are too different based on exon-only comparison, and kick out genes \
that overlap 'unalignable regions'

This leaves us with: \
       739 B73 genes with deletions in promoters relative to B97 \
       1256 B73 genes with insertions in promoters relative to B97 \
  and \
       11761 B73 genes that have no SVs (promoter or elsewhere)

The first thing I've done is look at the distribution of expression of these genes. \
I am pulling B73 expression from an existing TPM table (B73_full-tpm.tsv), only using one set of leaf tissues for initial analysis \
"B73_V11_middle".  For B97, ZM helped me map the B97 data back to B73, and again I'm pulling a leaf treatment only "leaf" - 4 reps \
I calculated the mean TPM for each gene and visualized with a violin plot -- note that the distribution of TPMs includes some very \
large values that I am not including here.  I also have not calculated a ratio of B73/B97 yet (though this is easy and on my list) \

Here I am subsetting to only include genes with TPM between 0 and 500 - it looks like the "No SV" genes are biased toward lower \
expression, but I need to watch out because the n is very different (see above)

![image](https://user-images.githubusercontent.com/43852873/211930843-f541b72f-bbb7-4623-8814-86a16e94fa25.png)




## Next Steps
Manisha is copying the B73 v B97 anchorwave prior to her cleanup into scratch - this allows me to do a bedtools intersect with gene +/- 1kb \
vs the anchorwave file to count SNPs and potentially small indels in my Prom_SV and No_SV gene sets

Look at B73vB97 TPM (x vs y) - color by insertion type/no sv, see if No SV genes have a higher Rsquared

Continue to get different tissue expression data

Begin to parse by: distance to TSS, Manisha category, etc

Are our promoter SV genes enriched for a specific type of gene(s)?
