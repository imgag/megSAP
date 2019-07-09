help:
	@cat Makefile
	
dosage_sensitive_disease_genes:
	GenesToBed -in dosage_sensitive_disease_genes.txt -mode gene -source ensembl | BedExtend -n 20 -out dosage_sensitive_disease_genes.bed

genes:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode gene -source ensembl | BedExtend -n 20 -out genes.bed

gene_exons:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode exon -source ensembl | BedMerge -out gene_exons.bed

gene_exons_pad20:
	BedExtend -n 20 -in gene_exons.bed | BedMerge -out gene_exons_pad20.bed
