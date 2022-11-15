help:
	@cat Makefile
	
create_somatic_tmb_tsg_bed:
	awk -F "\t" '$$6==1 {print $$2}' ../dbs/NCG7.0/NCG7.0_oncogene.tsv | GenesToApproved | cut -f1 | sort | uniq | GenesToBed -source ensembl -mode exon | BedMerge -merge_names > somatic_tmb_tsg.bed

create_somatic_tmb_blacklist:
	echo "HLA-A\nHLA-B" | GenesToBed -source ensembl -mode exon | BedMerge -merge_names > somatic_tmb_blacklist.bed

genes:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode gene -source ensembl | BedExtend -n 5000 -out genes.bed

gene_exons:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode exon -source ensembl | BedMerge -out gene_exons.bed

gene_exons_pad20:
	BedExtend -n 20 -in gene_exons.bed | BedMerge -out gene_exons_pad20.bed

housekeeping_genes_hg38:
	# based on https://academic.oup.com/nar/article/49/D1/D947/5871367
	wget -O - https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv | egrep -v "^Ensembl" | cut -d ";" -f2 | GenesToApproved | cut -f1 | sort | uniq > housekeeping_genes.txt
	GenesToBed -source ensembl -mode exon -in housekeeping_genes.txt | BedSort | BedMerge -out housekeeping_genes_hg38.bed