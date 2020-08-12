help:
	cat Makefile

# download example PRS from PGS Catalog:
download_prs:
	wget -O - http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/PGS000004.txt.gz | zcat > PGS000004.txt

# convert PRS file to VCF and calculate score distribution over all NGSD WGS samples:
setup_prs:
	php ../../../src/Tools/prs2vcf.php -in PGS000004.txt -out PGS000004.vcf