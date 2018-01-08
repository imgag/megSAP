help:
	@cat Makefile
	
copy_number_map_strict:
	wget -O - https://raw.githubusercontent.com/imgag/megSAP/master/data/dbs/CNPs/chrm_imb_map_str.gff3 | cut -f1,4,5 | grep -v "#" | BedSort > copy_number_map_strict.bed
