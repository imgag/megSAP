help:
	@cat Makefile

download:
	wget -O cancerhotspots.v2.maf.gz https://cbioportal-download.s3.amazonaws.com/cancerhotspots.v2.maf.gz

convert:
	php ../../../src/Install/db_converter_cancerhotspots.php -in cancerhotspots.v2.maf.gz -out cancerhotspots_snv.tsv --log cancerhotspots_snv.log

gzip:
	gzip -9 cancerhotspots_snv.tsv