## Is structural variation associated with tissue specific expression?

Gene Expression:  RNAseq data generated as part of the NAM genome sequencing project was used.  This includes RNAseq data for 10 tissue for each \
of the 25 NAM lines and B73. \
All RNAseq data was aligned to B73v5.

### SV Genes
Structural Variation:  The Anchorwave alignment of each NAM genome aligned with B73 was used to identify structural insertions/deletions relative to B73 \
Genes with structural in/dels in the 1kb promoter were identified. \
These genes were further filtered through a few steps:
1. genes must have one and only one SV in the 1kb promoter
2. genes must not have additional SV anywhere in the gene body (intron or exon)
3. genes must be expressed at greater than 4 TPM in at least one RNAseq run
4. genes can only have two haplotypes -- with and without the SV (there are some cases were there are multiple SV types across the NAM pop - so all have one SV, but not always the same SV)
5. haplotypes must include at least 3 members (at least 3 B73-like and at least 3 alternate) 

This leaves us with 296 SV genes -- not sure how this breaks down insertions relative to B73 vs deletions relative to B73

### A Control Gene Set
A control 'No-SV' gene set was identified - these genes have no structural variation across the NAM line either in the gene body or in the 1kb promoter
