# Correction of false duplications in GRCh38

About 1.3 MB of the GRCh38 reference was identified as falsely duplicated, including some medically relevant genes.  
For details, see this paper: https://www.nature.com/articles/s41587-021-01158-1

To fix this the regions are masked in the reference genome by replacing the sequence with 'N'.  
The regions are specified [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed).

## Changes in small variant calling

Several regions are uniquely mappable after the masking (MAPQ>0), but were not before (MAPQ=0).  
In these region, there will be new variant calls.

Mainly the regions also affected in CNV calling (see below) will contain changes in small variants.

## Changes in CNV calling

Only when reference samples for CNV calling already exists, which were analyzed with the uncorrected reference genome, there is a problem with CNVs.  
In that case you will see these artefacts:

* Masked regions will be shown with copy-number 0 after the change.
* New uniquely mappable regions are shown with copy-number 4.

This affects the following regions:

|chr  |start    |end      |genes                                                       |OMIM genes and phenotypes                                                         |
|-----|---------|---------|------------------------------------------------------------|----------------------------------------------------------------------------------|
|chr8 |30393636 |30407636 |RBPMS                                                       |                                                                                  |
|chr9 |70701717 |70737717 |TRPM3                                                       |                                                                                  |
|chr11|1945622  |2041622  |H19, LINC01219, MIR675, MRPL23, MRPL23-AS1                  |                                                                                  |
|chr13|86253328 |86269328 |                                                            |                                                                                  |
|chr13|111669328|111703328|                                                            |                                                                                  |
|chr13|111754328|111793328|                                                            |                                                                                  |
|chr16|19786345 |19801345 |IQCK                                                        |                                                                                  |
|chr16|34339345 |34500345 |                                                            |                                                                                  |
|chr16|34867345 |34896345 |                                                            |                                                                                  |
|chr16|34960345 |35072345 |CCNYL3, CLUHP11, LINC02184                                  |                                                                                  |
|chr21|5009983  |5165983  |                                                            |                                                                                  |
|chr21|5966983  |6160983  |                                                            |                                                                                  |
|chr21|6426983  |6579983  |                                                            |                                                                                  |
|chr21|6789983  |6933983  |                                                            |                                                                                  |
|chr21|7743983  |7865983  |                                                            |                                                                                  |
|chr21|13713983 |13718983 |                                                            |                                                                                  |
|chr21|13752983 |13798983 |FAM207CP, FEM1AP1, TERF1P1                                  |                                                                                  |
|chr21|34374983 |34494983 |C21ORF140, KCNE1, SMIM11, SMIM34                            |KCNE1: Jervell and Lange-Nielsen syndrome 2, Long QT syndrome                     |
|chr21|43035983 |43186983 |CBS, CRYAA, FRGCA, MRPL51P2, U2AF1                          |CBS: Thrombosis, Homocystinuria; CRYAA:Cataract                                   |
|chr21|43376983 |43570983 |H2BC12L, HSF2BP, LINC00313, LINC00319, RPL31P1, SIK1        |SIK1: Developmental and epileptic encephalopathy; HSF2BP:Premature ovarian failure|
|chr21|44095983 |44252983 |DNMT3L, DNMT3L-AS1, GATD3, ICOSLG, LINC01678, PWP2, TRAPPC10|                                                                                  |
|chr22|18885468 |18939468 |DGCR6, FAM230F, PRODH                                       |PRODH: Schizophrenia, Hyperprolinemia                                             |
|chrX |37084895 |37098895 |                                                            |                                                                                  |

## Checking the genome version

To check if your genome is correctly masked, use this command:  

	> samtools faidx [genome] chr21:5010000-5166246

If the output consists of 'N' only, the genome is correctly masked.

To check if a genome/exome BAM/CRAM file was created with masked reference genome, use this command:  

	> samtools view [bam/cram] chr21:5966593-6161371

If the output is empty, a masked reference genome was used.

## Internal documentation

GitHub issue: https://github.com/imgag/megSAP/issues/147  
Benchmarks were performed in the folder: /mnt/storage2/users/ahsturm1/scripts/2022\_09\_12\_grch38\_false\_duplications
