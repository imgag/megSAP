help:
	cat Makefile

# download example PRS from PGS Catalog:
download_prs:
	wget -O - http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/PGS000004.txt.gz | zcat > PGS000004.txt


# GRCh38
create_pgs_vcf_GRCh38:
	php ../../../src/Tools/prs2vcf.php -pgs PGS000004.txt -out PGS000004_GRCh37_nodist.vcf -build GRCh37 -skip_percentiles

lift_vcf_files:
	echo "Lift the GRCh37 VCFs. (Often requires manual post-processing!)"

# calculate score distribution over all NGSD WGS samples for GRCh38:
# (requires manually lifted VCF PRS file since there is no GRCh38 version of PGS000004 available yet)
calculate_percentiles_GRCh38:
	php ../../../src/Tools/prs2vcf.php -vcf PGS000004_GRCh38_lifted.vcf -out PGS000004_GRCh38.vcf -build GRCh38

# OCAC_36
# Requires installed CrossMap and chain files
# File paths may have to be adjusted 
OCAC_36:
	php ../../../src/Tools/vcf_annotate_coords.php -in OCAC_36_GRCh37.vcf -out OCAC_36_GRCh37_ann.vcf
	CrossMap.py vcf /mnt/share/opt/liftOver/hg19ToHg38.over.chain.gz OCAC_36_GRCh37_ann.vcf /mnt/share/data/genomes/GRCh38.fa OCAC_36_GRCh38_lifted.vcf
	echo "WARNING: Check lifting manually and lift remaining by hand!"

OCAC_36_calculate_percentiles:
	php ../../../src/Tools/prs2vcf.php -vcf OCAC_36_GRCh38_lifted.vcf -out OCAC_36_GRCh38.vcf -build GRCh38
