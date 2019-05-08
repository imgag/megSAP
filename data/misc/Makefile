help:
	cat Makefile
	
cytoBand:
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip > cytoBand.txt

cn_polymorphisms_zarrei:
	wget -O - https://raw.githubusercontent.com/imgag/megSAP/master/data/dbs/CNPs/chrm_imb_map_str.gff3 | cut -f1,4,5 | grep -v "#" | BedSort > cn_polymorphisms_zarrei.bed

cn_polymorphisms_demidov:
	TODO
	
cn_pathogenic:
	cd /mnt/share/data/dbs/CNVs_pathogenic/ && make convert deploy