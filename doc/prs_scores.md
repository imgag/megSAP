The megSAP pipeline calculates PRS scores in the WGS single sample analysis. 

## Requirements
To calculate the PRS score of a WGS sample the pipeline needs specialized VCFs in the folder `data/misc/prs/`. For each PRS score a separate VCF is required. The VCF contains all variants which contributes to the PRS score and the corresponding weights in the INFO part of the VCF. Additionally some information about the PRS (like trait, paper,...) is stored in the header part of the VCF. An example of the VCF containing a PRS score for breast cancer is already included in the git repository at `data/misc/prs/PGS000004.vcf`:
```
##fileformat=VCFv4.2
##fileDate=20200723
##pgs_id=PGS000004
##trait=Breast Cancer
##build=GRCh38
##n_var=313
##pgp_id=PGP000002
##citation=Mavaddat N et al. Am J Hum Genet (2018). doi:10.1016/j.ajhg.2018.11.002
##sample_count=1344
##percentiles=-1.828376,-1.612056,-1.488528,-1.406704,-1.33842...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100880328	.	A	T	.	.	WEIGHT=0.0373
chr1	10566215	.	A	G	.	.	WEIGHT=-0.0586
chr1	110198129	.	CAAA	C	.	.	WEIGHT=0.0458
...
```
To create such a VCF file the php script `prs2vcf.php` in `src/Tools`can be used. It takes a PRS text file from [https://www.pgscatalog.org/](https://www.pgscatalog.org/) as input and generates a VCF which can be used in the megSAP pipeline. If the script has also access to the NGSD it automatically computes the percentiles for the given PRS score based on the matching samples in the NGSD.

## PRS calculation
The megSAP pipeline will automatically generate a TSV file during each WGS single sample analysis 