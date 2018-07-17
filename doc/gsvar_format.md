# megSAP - GSvar format

In addition to the standard VCF format, megSAP creates variant lists in GSvar format in most analysis pipelines.

The GSvar format is a tab-separated text format which can be opened and filtered for example in Excel.

## Variant filtering

There are several ways to filter variant lists in GSvar format to identify rare pathogenic variants:

1. The file can be filtered in Microsoft Excel.
2. The file can be filtered on the Linux command line using the [VariantFilterAnnotations tool](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantFilterAnnotations.md) of ngs-bits.
3. The file can be filtered interactively on Windows using the GSvar tool if ngs-bits (more documentation on GSvar coming soon).

## Format details

#### Meta data header lines

Similar to VCF format, the GSvar format can contain meta data in header lines that start with `##`:

- The `##ANALYSISTYPE=` header specifies the analysis type. It can be present only once. Valid analysis types are:  
	`GERMLINE_SINGLESAMPLE`, `GERMLINE_MULTISAMPLE`,  `GERMLINE_TRIO`,  `SOMATIC_SINGLESAMPLE`,  `SOMATIC_PAIR`  

- The `##PIPELINE=` header specifies the analysis pipeline used to perform the analysis.
	
- The `##SAMPLE=` header specifies the analyzed sample name(s) and can contain additional annotations. It can be present several times.   

- The `##DESCRIPTION=` header describes columns. It can be present several times.    

- The `##FILTER=` header describes entries in the `filter` column. It can be present several times.  

 
#### Main header line and variant lines

After the meta data header line, the main header line and the variant lines follow:

	#chr	start	end	ref	obs	genotype	filter	quality	gene	variant_type	coding_and_splicing	RepeatMasker	dbSNP	1000g	ExAC	ExAC_hom	Kaviar	phyloP	Sift	MetaLR	PolyPhen2	FATHMM	CADD	OMIM	ClinVar	HGMD	COSMIC	ihdb_allsys_hom	ihdb_allsys_het	classification	classification_comment	validated	comment	gene_info
	chr1	27682481	27682481	G	A	het	off-target	QUAL=2181;DP=169;AF=0.51;MQM=60	MAP3K6	intron	MAP3K6:NM_004672.4:intron:MODIFIER:exon27/28:c.3711+36C>T:,MAP3K6:NM_001297609.1:intron:MODIFIER:exon26/27:c.3687+36C>T:		rs12569127	0.1903	0.2656	4554/3529/231	0.2588											95	396					MAP3K6 (inh=n/a pLI=0.00)
	...


[back to the start page](../README.md)







