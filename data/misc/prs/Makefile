help:
	cat Makefile


################################### OCAC36 (ovarian cancer) ###################################

OCAC36_liftover:
	php ../../../src/Tools/vcf_annotate_coords.php -in OCAC_36_GRCh37.vcf -out OCAC_36_GRCh37_ann.vcf
	CrossMap.py vcf /mnt/storage2/megSAP/tools/liftOver/hg19ToHg38.over.chain.gz OCAC_36_GRCh37_ann.vcf /mnt/storage2/megSAP/data/genomes/GRCh38.fa OCAC_36_GRCh38_lifted.vcf
	VcfSort -in OCAC_36_GRCh38_lifted.vcf -out OCAC_36_GRCh38_lifted.vcf
	@echo "WARNING: Check lifting manually and lift remaining by hand!"

OCAC36_calculate_percentiles:
	php ../../../src/Tools/prs2vcf.php -vcf OCAC_36_v1_fixed.vcf -out OCAC_36_GRCh38.vcf -build GRCh38 -exclude_disease_group Neoplasms -processing_system TruSeqPCRfree

################################### BRIDGES_306 (breast cancer) ###################################
BRIDGES_306_download:
	# create header information:
	echo "#pgs_id=BRIDGES_306" > BRIDGES_306_GRCh37.txt
	echo "#trait_reported=Breast cancer" >> BRIDGES_306_GRCh37.txt
	echo "#variants_number=306" >> BRIDGES_306_GRCh37.txt
	echo "#pgp_id=undefined" >> BRIDGES_306_GRCh37.txt
	echo "#citation=CanRisk (https://canrisk.atlassian.net/wiki/spaces/FAQS/pages/35979266/What+variants+are+used+in+the+PRS)" >> BRIDGES_306_GRCh37.txt
	# download score
	wget -O - https://canrisk.atlassian.net/wiki/download/attachments/35979266/BRIDGES_306_PRS.prs | sed "s/,/	/g" >> BRIDGES_306_GRCh37.txt
	# adapt/fix header
	sed -i "s/Chromosome/chr_name/g" BRIDGES_306_GRCh37.txt
	sed -i "s/Position/chr_position/g" BRIDGES_306_GRCh37.txt
	sed -i "s/Reference_Allele/reference_allele/g" BRIDGES_306_GRCh37.txt
	sed -i "s/Effect_Allele_Frequency/allelefrequency_effect/g" BRIDGES_306_GRCh37.txt
	sed -i "s/Effect_Allele/effect_allele/g" BRIDGES_306_GRCh37.txt	
	sed -i "s/Log_Odds_Ratio/effect_weight/g" BRIDGES_306_GRCh37.txt
	sed -i "s/alpha = 0.498/#genome_build=GRCh37/g" BRIDGES_306_GRCh37.txt

BRIDGES_306_to_vcf:
	php ../../../src/Tools/prs2vcf.php -pgs BRIDGES_306_GRCh37.txt -out BRIDGES_306_GRCh37.vcf -build GRCh37 -skip_percentiles
	@echo "Notice: Adding header lines:"
	sed -i "s/Chromosome/chr_name/g" BRIDGES_306_GRCh37.vcf

BRIDGES_306_liftover:
	php ../../../src/Tools/vcf_annotate_coords.php -in BRIDGES_306_GRCh37.vcf -out BRIDGES_306_GRCh37_ann.vcf
	CrossMap.py vcf --no-comp-alleles /mnt/storage2/megSAP/tools/liftOver/hg19ToHg38.over.chain.gz BRIDGES_306_GRCh37_ann.vcf /mnt/storage2/megSAP/data/genomes/GRCh38.fa BRIDGES_306_GRCh38_lifted_tmp.vcf
	# transfer header info from GRCh37 and remove contig header lines:
	egrep "^##" BRIDGES_306_GRCh37_ann.vcf | egrep -v "^##INFO=" > BRIDGES_306_GRCh38_lifted.vcf
	egrep -v "^##contig" BRIDGES_306_GRCh38_lifted_tmp.vcf | egrep -v "^##fileformat" >> BRIDGES_306_GRCh38_lifted.vcf
	rm BRIDGES_306_GRCh38_lifted_tmp.vcf
	VcfLeftNormalize -in BRIDGES_306_GRCh38_lifted.vcf -out BRIDGES_306_GRCh38_lifted.vcf
	VcfSort -in BRIDGES_306_GRCh38_lifted.vcf -out BRIDGES_306_GRCh38_lifted.vcf
	VcfCheck -in BRIDGES_306_GRCh38_lifted.vcf -lines 0
	@echo "WARNING: Check lifting manually and lift remaining by hand!"

BRIDGES_306_calculate_percentiles_GRCh38:
	php ../../../src/Tools/prs2vcf.php -vcf BRIDGES_306_v2_fixed.vcf -out BRIDGES_306_GRCh38.vcf -build GRCh38 -exclude_disease_group Neoplasms -processing_system TruSeqPCRfree
