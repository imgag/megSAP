help:
	@cat Makefile
	
dosage_sensitive_disease_genes:
	GenesToBed -in dosage_sensitive_disease_genes.txt -mode gene -source ensembl | BedExtend -n 20 -out dosage_sensitive_disease_genes.bed

genes:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode gene -source ensembl | BedExtend -n 20 -out genes.bed

genes_exons:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode exon -source ensembl | BedExtend -n 20 | BedMerge -out genes_exons.bed
