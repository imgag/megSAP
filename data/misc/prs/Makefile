help:
	cat Makefile

# download example PRS from PGS Catalog:
download_prs:
	wget -O - http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/PGS000004.txt.gz | zcat > PGS000004.txt

# convert PRS file to VCF and calculate score distribution over all NGSD WGS samples:
setup_prs_GRCh37:
	php ../../../src/Tools/prs2vcf.php -pgs PGS000004.txt -out PGS000004_GRCh37.vcf -build GRCh37

# calculate score distribution over all NGSD WGS samples for GRCh38:
# (requires manually lifted VCF PRS file since there is no GRCh38 version of PGS000004 available yet)
setup_prs_GRCh38:
	php ../../../src/Tools/prs2vcf.php -pgs PGS000004.txt -out PGS000004_GRCh37_nodist.vcf -build GRCh37 -skip_percentiles
	php ../../../src/Tools/prs2vcf.php -vcf PGS000004_GRCh38_manual_lifted.vcf -out PGS000004_GRCh38.vcf -build GRCh38