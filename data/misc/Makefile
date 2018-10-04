help:
	cat Makefile
	
cytoBand:
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip > cytoBand.txt

copy_number_map_strict:
	wget -O - https://raw.githubusercontent.com/imgag/megSAP/master/data/dbs/CNPs/chrm_imb_map_str.gff3 | cut -f1,4,5 | grep -v "#" | BedSort > copy_number_map_strict.bed
