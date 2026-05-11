help:
	@cat Makefile

download:
	wget -O cancerhotspots.v2.maf.gz https://cbioportal-download.s3.amazonaws.com/cancerhotspots.v2.maf.gz

convert:
	php ../../../src/Install/db_converter_cancerhotspots.php -in cancerhotspots.v2.maf.gz -out cancerhotspots_snv_v2.tsv --log cancerhotspots_snv.log


v3_addon:
	#download the new hotspots from the website www.cancerhotspots.org while searching for the new publication "Bandlamudi et al., 2026"
	#download saved in hotspots_v3_website.tsv
	cut -f1 hotspots_v3_website.tsv | uniq > genes.txt
	#GenesToTranscripts -mode best -in genes.txt -out transcript_info.txt
	cut -f2 transcript_info.txt | TranscriptToProtein -out protein_sequence.tsv
	echo "#gene\ttranscript\taa_sequence" > aa_sequence.tsv
	tail -n +2 genes.txt | paste - protein_sequence.tsv >> aa_sequence.tsv
	
	#as our best transcript wasn't always the one used in paper the transcripts for PTPRT and TGFBR2 where changed by hand:
	#PTPRT: 	ENST00000373187 changed to ENST00000373198
	#TGFBR2:	ENST00000295754 changed to ENST00000359013
	
	php ../../../src/Install/db_cancerhotspots_check_aminoacid.php -hotspot_table hotspots_v3_website.tsv -sequences aa_sequence.tsv -out db_v3_addon.tsv
	cat cancerhotspots_snv_v2.tsv db_v3_addon.tsv > cancerhotspots_snv_unsorted.tsv
	#sort but keep header at the top
	head -n 1 cancerhotspots_snv_unsorted.tsv > cancerhotspots_snv.tsv
	tail -n +2 cancerhotspots_snv_unsorted.tsv | sort >> cancerhotspots_snv.tsv
	
	#remove temporary files
	rm genes.txt transcript_info.txt protein_sequence.tsv aa_sequence.tsv db_v3_addon.tsv cancerhotspots_snv_unsorted.tsv 

gzip:
	gzip -9 cancerhotspots_snv.tsv