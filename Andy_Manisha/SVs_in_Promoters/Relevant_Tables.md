## Some relevant tables/files

### I'm dropping these into: /home/read0094/read0094/For_Manisha

1. NAM_TE_Exp_SV_Tablev2.txt\
Rows = Genes\
Columns = geneID, syntenic (y/n), PromoterSV?, PromoterSVfiltered?, No_SV_genes?, Structural SV in promoter, Max Expression, Is there a Struct SV

`````
geneID	syntenic	PromSV	PromSV_filtered	No_SV_genes	AnyPromTE	StructPromTE	MaxTPM	Struct_SV_TE	Any_SV_TE
Zm00001eb000010	Not_Syntenic	Not_SinglePromSV	Not_PromSV_filtered	Gene_has_PIE_SV	TE_in_Promoter	No_Struct_Prom_TE	5961.22	NA	NA
Zm00001eb000020	Not_Syntenic	Not_SinglePromSV	Not_PromSV_filtered	No_SV_gene	TE_in_Promoter	No_Struct_Prom_TE	16506.00	NA	NA
Zm00001eb000050	Not_Syntenic	Not_SinglePromSV	Not_PromSV_filtered	Gene_has_PIE_SV	TE_in_Promoter	Struct_TE_in_Promoter	382.43	NA	NA
Zm00001eb000060	Not_Syntenic	Not_SinglePromSV	Not_PromSV_filtered	Gene_has_PIE_SV	TE_in_Promoter	No_Struct_Prom_TE	3314.58	NA	NA
`````

2. norm_NAM_RNAseq_counts.txt\
Includes average expression for each RNAseq dataset
Rows = Genes\
Expression “count” in a tissue/rep = columns

(many additional columns)
`````
	B73_root_1_ERR3773807	B73_root_2_ERR3773808	B73_shoot_1_ERR3773809	B73_shoot_2_ERR3773810	B73_leaf.base_1_ERR3773811	B73_leaf.base_2_ERR3773812	B73_leaf_1_ERR3773813
Zm00001eb000010	256.35	203.62	848.01	905.83	589.96	372.07	1127.82
Zm00001eb000020	3944.06	3593.81	2341.92	2367.16	1860.58	3341.37	723.56
Zm00001eb000050	4.52	8.68	0.00	0.00	0.00	12.08	0.00
Zm00001eb000060	653.70	516.07	506.36	435.96	1060.53	844.00	780.21
`````

3. NAM_all_gene_EXP_mean_CV.txt\
Rows = Genes\
Mean, SD, and CV for expression data across all tissues (single value for all tissues)

`````
geneID	B97_mean	B97_sd	B97_cv	CML103_mean	CML103_sd	CML103_cv	CML228_mean	CML228_sd	CML228_cv	CML247_mean
Zm00001eb000010	1019.10	791.54	0.78	880.34	675.54	0.77	1163.60	788.51	0.68	827.60
Zm00001eb000020	2811.01	2913.76	1.04	3308.93	2297.04	0.69	3086.13	2472.65	0.80	2573.44
`````

4. All_gene_cultivar_combos_filtered_byhap_extendedv2.txt \
Rows = Gene_cultivar \
Columns = haplotype number, gene, Promoter_chromosome, promoter 5prime coord, promoter 3prime coord, SV_block ID, SV_chromosome, SV 5prime coord, SV 3prime coord, Category from Manisha’s data
`````
gene_lineage	haplotype_num	V6	Prom_chr	Prom_5prime	Prom_3prime	SV_block	SV_chr	SV_5prime	SV_3prime	Man_category
Zm00001eb000150_B97	1	Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_23	chr1	322874	328037	Category_5
Zm00001eb000150_CML103	1	Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_17	chr1	322874	328037	Category_5
Zm00001eb000150_CML228	1	Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_17	chr1	322874	328037	Category_5
Zm00001eb000150_CML247	1	Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_37	chr1	322874	328037	Category_5
Zm00001eb000150_CML277	1	Zm00001eb000150	chr1	327086	328086	chr1_AW_BlockID_29	chr1	322874	328037	Category_5
`````

5. NAM_SVs_details.txt \
Rows = geneID \
Columns = PIED feature, feature chromosome, feature start, feature end, feature length, AW block ID, AW block chr, AW block start, AW block end, MAM category, overlap length, Description (B73 INS/DEL relative to which NAM), gene_feature concatenated \
**Notes: This version of the table includes the 1kb downstream region -- There is some feature overlap here due to annotations that are close/overlapping**

`````
geneID	PIE	gene_PIE_chr	gene_PIE_start	gene_PIE_end	gene_PIE_len	AW_blockID	AW_chr	AW_start	AW_end	MAM_category	overlap_len	Description	gene_feature
Zm00001eb000050	intron	chr1	109952	113761	3810	chr1_AW_BlockID_9	chr1	110302	111609	Category_5	1307	B73_Ins_Rel_B97	Zm00001eb000050_intron
Zm00001eb000050	intron	chr1	109952	113761	3810	chr1_AW_BlockID_11	chr1	113391	118693	Category_4	370	B73_Ins_Rel_B97	Zm00001eb000050_intron
Zm00001eb000050	exon	chr1	113762	114382	620	chr1_AW_BlockID_11	chr1	113391	118693	Category_4	620	B73_Ins_Rel_B97	Zm00001eb000050_exon
Zm00001eb000050	Promoter	chr1	114383	115383	1000	chr1_AW_BlockID_11	chr1	113391	118693	Category_4	1000	B73_Ins_Rel_B97	Zm00001eb000050_Promoter
Zm00001eb000140	Downstream	chr1	324306	325306	1000	chr1_AW_BlockID_23	chr1	322874	328037	Category_5	1000	B73_Ins_Rel_B97	Zm00001eb000140_Downstream
`````
