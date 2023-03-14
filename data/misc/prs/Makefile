help:
	cat Makefile


################################### OCAC36 (ovarian cancer) ###################################

OCAC36_liftover:
	php ../../../src/Tools/vcf_annotate_coords.php -in OCAC_36_GRCh37.vcf -out OCAC_36_GRCh37_ann.vcf
	CrossMap.py vcf /mnt/storage2/megSAP/tools/liftOver/hg19ToHg38.over.chain.gz OCAC_36_GRCh37_ann.vcf /mnt/storage2/megSAP/data/genomes/GRCh38.fa OCAC_36_GRCh38_lifted.vcf
	VcfSort -in OCAC_36_GRCh38_lifted.vcf -out OCAC_36_GRCh38_lifted.vcf
	@echo "WARNING: Check lifting manually and lift remaining by hand!"

OCAC36_calculate_percentiles:
	php ../../../src/Tools/prs2vcf.php -vcf OCAC_36_GRCh38_lifted.vcf -out OCAC_36_GRCh38.vcf -build GRCh38

################################### PGS000004 (breast cancer) ###################################

PGS000004_download:
	wget -O - http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/PGS000004.txt.gz | zcat > PGS000004.txt

PGS000004_to_vcf:
	php ../../../src/Tools/prs2vcf.php -pgs PGS000004.txt -out PGS000004_GRCh37.vcf -build GRCh37 -skip_percentiles

PGS000004_liftover:
	php ../../../src/Tools/vcf_annotate_coords.php -in PGS000004_GRCh37.vcf -out PGS000004_GRCh37_ann.vcf
	CrossMap.py vcf /mnt/storage2/megSAP/tools/liftOver/hg19ToHg38.over.chain.gz PGS000004_GRCh37_ann.vcf /mnt/storage2/megSAP/data/genomes/GRCh38.fa PGS000004_GRCh38_lifted.vcf
	VcfSort -in PGS000004_GRCh38_lifted.vcf -out PGS000004_GRCh38_lifted.vcf
	@echo "WARNING: Check lifting manually and lift remaining by hand!"

calculate_percentiles_GRCh38:
	php ../../../src/Tools/prs2vcf.php -vcf PGS000004_GRCh38_lifted.vcf -out PGS000004_GRCh38.vcf -build GRCh38
