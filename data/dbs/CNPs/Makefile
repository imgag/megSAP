help:
	@cat Makefile
	
copy_number_map_strict:
	cut -f1,4,5 chrm_imb_map_str.gff3 | grep -v "#" | BedSort > copy_number_map_strict.bed
