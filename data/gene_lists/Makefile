help:
	@cat Makefile
	
dosage_sensitive_disease_genes:
	GenesToBed -in dosage_sensitive_disease_genes.txt -mode gene -source ensembl | BedExtend -n 20 -out dosage_sensitive_disease_genes.bed

genes:
	GenesToBed -in genes.txt -mode gene -source ensembl | BedExtend -n 20 -out genes.bed
