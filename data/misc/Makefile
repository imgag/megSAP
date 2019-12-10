help:
	cat Makefile
	
cytoBand:
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip > cytoBand.txt

cnps_genomes_imgag:
	cp /mnt/SRV018/users/ahsturm1/Sandbox/ClinCNV/cnps_genomes/cnps_genomes_imgag.igv .
	cut -f1,2,3,5 cnps_genomes_imgag.igv | grep -v "Chromosome" > cnps_genomes_imgag.bed

af_genomes_imgag:
	NGSDExportCnvTrack -system TruSeqPCRfree -min_af 0.02 -min_dp 30 -max_cnvs 2000 -caller_version 'v1.16.4' -out af_genomes_imgag.igv -stats af_genomes_imgag.log
	php merge_cnps.php  af_genomes_imgag.igv cnps_genomes_imgag.igv > tmp.igv
	mv tmp.igv af_genomes_imgag.igv
	cut -f1,2,3,5 af_genomes_imgag.igv | grep -v "Chromosome" > af_genomes_imgag.bed
	
cn_pathogenic:
	cd /mnt/share/data/dbs/CNVs_pathogenic/ && make convert deploy
	
nobase_regions_GRCh37:
	php extract_nobase.php GRCh37 > nobase_regions_GRCh37.bed

pfam_description:
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz | zcat | grep -e '^#=GF AC' -e '^#=GF DE' | tr -d '\n' | sed 's/#=GF AC   /\n/g' | sed 's/\.[0-9]\+#=GF DE   /\t/g' | sed '/^$$/d' | sed '1 i\#PFAM_ID\tDESCRIPTION' | tr '\n\n' '\n'> pfam_description.tsv
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.dead.gz | zcat | grep -e '^#=GF AC' -e '^#=GF FW' | tr -d '\n' | sed 's/#=GF AC   /\n/g' | sed 's/#=GF FW   /\t/g' | sed '/^$$/d' | sed '1 i\#OLD_PFAM_ID\tNEW_PFAM_ID' | tr '\n\n' '\n'> pfam_replacements.tsv
	