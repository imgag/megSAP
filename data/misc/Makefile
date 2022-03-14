help:
	cat Makefile
	
cytoBand:
	wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz | gunzip > cytoBand.txt

cnps_genomes_imgag:
	cp /mnt/users/ahsturm1/Sandbox/ClinCNV/cnps_genomes/cnps_genomes_imgag.igv ./cnps_genomes_imgag.igv
	cut -f1,2,3,5 cnps_genomes_imgag.igv | grep -v "Chromosome" > cnps_genomes_imgag.bed

af_genomes_imgag:
	/mnt/storage2/GRCh38/share/opt/ngs-bits-current/NGSDExportCnvTrack -system TruSeqPCRfree -min_af 0.02 -min_dp 30 -max_cnvs 2000 -caller_version 'v1.17.0' -out af_genomes_imgag.igv -stats af_genomes_imgag.log
	php merge_cnps.php af_genomes_imgag.igv cnps_genomes_imgag.igv > tmp.igv
	mv tmp.igv af_genomes_imgag.igv
	cut -f1,2,3,5 af_genomes_imgag.igv | grep -v "Chromosome" > af_genomes_imgag.bed

af_ssHAEv7_imgag:
	/mnt/storage2/GRCh38/share/opt/ngs-bits-current/NGSDExportCnvTrack -system ssHAEv7 -min_af 0.02 -min_dp 80 -max_cnvs 200 -caller_version 'v1.17.0' -out af_ssHAEv7_imgag.igv -stats af_ssHAEv7_imgag.log
	cut -f1,2,3,5 af_ssHAEv7_imgag.igv | grep -v "Chromosome" > af_ssHAEv7_imgag.bed
	
nobase_regions:
	php extract_nobase.php GRCh38 > nobase_regions.bed

pfam_description:
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz | zcat | grep -e '^#=GF AC' -e '^#=GF DE' | tr -d '\n' | sed 's/#=GF AC   /\n/g' | sed 's/\.[0-9]\+#=GF DE   /\t/g' | sed '/^$$/d' | sed '1 i\#PFAM_ID\tDESCRIPTION' | tr '\n\n' '\n'> pfam_description.tsv
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.dead.gz | zcat | grep -e '^#=GF AC' -e '^#=GF FW' | tr -d '\n' | sed 's/#=GF AC   /\n/g' | sed 's/#=GF FW   /\t/g' | sed '/^$$/d' | sed '1 i\#OLD_PFAM_ID\tNEW_PFAM_ID' | tr '\n\n' '\n'> pfam_replacements.tsv
	
telomere_centromere:
	# download Gap (filter: gap.type = 'telomere') and Centromeres table from https://genome.ucsc.edu/cgi-bin/hgTables
	egrep -v "^#" telomere.tsv | cut -f2-4,8 > telomere.bed
	egrep -v "^#" centromere.tsv | cut -f2-5 | BedExtend -n 500000 | BedMerge | BedShrink -n 500000 | awk 'BEGIN {OFS="\t"}; {print $$0, "centromere"}' > centromere.bed 
	cat telomere.bed centromere.bed | BedSort -out centromer_telomer.bed
	rm telomere.bed centromere.bed

circos_karyotype_files:
	echo "These files where created manually by merging the chromosome regions and the centromer_telomer.bed (extended by 3 columns and replaced names by color)"
	
create_custom_somatic_panel:
	BedIntersect -in /mnt/storage2/GRCh38/share/data/enrichment/ssSC_v5.bed -in2 /mnt/storage2/GRCh38/share/data/enrichment/twistCustomExomeV2_2021_12_14.bed -mode in -out somatic_custom_panel.bed
	BedAnnotateGenes -in somatic_custom_panel.bed | cut -f4 | sort | uniq > somatic_custom_panel_genes.txt
