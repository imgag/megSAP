help:
	cat Makefile
	
cytoBand:
	wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz | gunzip > cytoBand.txt

cnps_genomes_imgag:
	cp /mnt/SRV018/users/ahsturm1/Sandbox/ClinCNV/cnps_genomes/cnps_genomes_imgag.igv .
	cut -f1,2,3,5 cnps_genomes_imgag.igv | grep -v "Chromosome" > cnps_genomes_imgag.bed

af_genomes_imgag:
	/mnt/share/opt/ngs-bits-current/NGSDExportCnvTrack -system TruSeqPCRfree -min_af 0.02 -min_dp 30 -max_cnvs 2000 -out af_genomes_imgag.igv -stats af_genomes_imgag.log
	#TODO take PAR region from this export as long as ClinCNV does not handle PAR correctly
	#/mnt/share/opt/ngs-bits-current/NGSDExportCnvTrack -system TruSeqPCRfree -min_af 0.02 -min_dp 30 -max_cnvs 2000 -out af_genomes_imgag_no_males.igv -stats af_genomes_imgag_no_males.log -skip_males
	#grep -v "Chromosome" af_genomes_imgag_no_males.igv | BedIntersect -in2 par_regions_hg19.bed > par.bed
	cut -f1,2,3,5 af_genomes_imgag.igv | grep -v "Chromosome" > af_genomes_imgag.bed
	
cn_pathogenic:
	cd /mnt/share/data/dbs/CNVs_pathogenic/ && make convert deploy
	
nobase_regions_GRCh37:
	php extract_nobase.php GRCh37 > nobase_regions_GRCh37.bed
	