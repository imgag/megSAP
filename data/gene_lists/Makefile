help:
	@cat Makefile
	
create_somatic_tmb_tsg_bed:
	awk -F "\t" '$$6==1 {print $$2}' ../dbs/NCG6.0/NCG6.0_oncogene.tsv | GenesToApproved | cut -f1 | sort | uniq | GenesToBed -source ensembl -mode exon | BedMerge -merge_names > somatic_tmb_tsg.bed

create_somatic_tmb_blacklist:
	echo "HLA-A\nHLA-B" | GenesToBed -source ensembl -mode exon | BedMerge -merge_names > somatic_tmb_blacklist.bed

dosage_sensitive_disease_genes:
	GenesToBed -in dosage_sensitive_disease_genes.txt -mode gene -source ensembl | BedExtend -n 20 -out dosage_sensitive_disease_genes.bed

genes:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode gene -source ensembl | BedExtend -n 20 -out genes.bed

gene_exons:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode exon -source ensembl | BedMerge -out gene_exons.bed

gene_exons_pad20:
	BedExtend -n 20 -in gene_exons.bed | BedMerge -out gene_exons_pad20.bed
