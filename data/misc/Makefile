help:
	cat Makefile
	
cytoBand:
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip > cytoBand.txt

cn_polymorphisms_300genomes_imgag:
	cp /mnt/SRV018/users/ahsturm1/Sandbox/ClinCNV/cnps_genomes/cn_polymorphisms_300genomes_imgag.igv .
	cut -f1,2,3 cn_polymorphisms_300genomes_imgag.igv | grep -v "Chromosome" | BedMerge > cn_polymorphisms_300genomes_imgag.bed

cn_polymorphisms_700genomes_pcawg:
	cp /mnt/SRV018/users/ahsturm1/Sandbox/ClinCNV/cnps_genomes/cn_polymorphisms_700genomes_pcawg.igv .
	cut -f1,2,3 cn_polymorphisms_700genomes_pcawg.igv | grep -v "Chromosome" | BedMerge > cn_polymorphisms_700genomes_pcawg.bed
	
cn_pathogenic:
	cd /mnt/share/data/dbs/CNVs_pathogenic/ && make convert deploy
	
nobase_regions_GRCh37:
	php extract_nobase.php GRCh37 > nobase_regions_GRCh37.bed
	