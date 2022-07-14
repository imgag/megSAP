# Tool and database overview

## Databases

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ClinGen                |pipeline                |12.2021 (latest - not versioned)|                                                                                                    |https://ftp.clinicalgenome.org/                                         |
|RepeatMasker           |pipeline                |12.2021 (140131)                |                                                                                                    |http://www.repeatmasker.org/species/hg.html                             |
|ClinVar (SNVs and CNVs)|pipeline, NGSDImportHPO |07.2022 (20220702)              |update GSvar IGV file                                                                               |https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2021/   |
|HGNC                   |pipeline, NGSDImportHGNC|12.2021 (latest - not versioned)|                                                                                                    |https://ftp.ebi.ac.uk/pub/databases/genenames/                          |
|gnomAD (genome)        |pipeline                |05.2021 (3.1.1)                 |3.1.2 available!!!                                                                                  |http://gnomad.broadinstitute.org/downloads                              |
|gnomAD (constraints)   |NGSDImportGeneInfo      |01.2021 (2.1.1)                 |                                                                                                    |http://gnomad.broadinstitute.org/downloads                              |
|phyloP - for VEP       |pipeline                |12.2021 (05.2015)               |                                                                                                    |http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/            |
|CADD - for VEP         |pipeline                |12.2021 (1.6)                   |                                                                                                    |http://cadd.gs.washington.edu/download                                  |
|REVEL - for VEP        |pipeline                |12.2021 (1.3)                   |                                                                                                    |https://sites.google.com/site/revelgenomics/downloads                   |
|OMIM                   |pipeline, NGSDImportHPO |07.2022 (latest - not versioned)|                                                                                                    |https://omim.org/downloads/                                             |
|HGMD (SNVs and CNVs)   |pipeline                |07.2021 (2022.2)                |update GSvar IGV file                                                                               |https://apps.ingenuity.com/ingsso/login                                 |
|Ensembl                |NGSDImportEnsembl       |12.2021 (GRCh38.105)            |                                                                                                    |https://ftp.ensembl.org/pub/release-105/gff3/homo_sapiens/              |
|HPO                    |NGSDImportHPO           |12.2021 (latest - not versioned)|Bug: https://github.com/obophenotype/human-phenotype-ontology/issues/4916                           |https://hpo.jax.org/app/                                                |
|GenCC                  |NGSDImportHPO           |12.2021 (latest - not versioned)|                                                                                                    |https://search.thegencc.org/download                                    |
|DECIPER                |NGSDImportHPO           |12.2021 (22_12_2021)            |                                                                                                    |https://www.deciphergenomics.org/about/downloads/data                   |
|ORPHA                  |NGSDImportORPHA         |12.2021 (latest - not versioned)|                                                                                                    |https://github.com/Orphanet/Orphadata.org/                              |


## Databases for somatic pipelines

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|CancerHotspots         |pipeline (somatic)      |01.2021                         |version is final and does not change                                                                |https://www.cancerhotspots.org                                          |
|NCG6.0                 |pipeline (somatic)      |01.2021                         |                                                                                                    |http://ncg.kcl.ac.uk/                                                   |
|COSMIC CMC             |pipeline (somatic)      |01.2021                         |                                                                                                    |https://cancer.sanger.ac.uk/cmc                                         |


## Tools

|tool                   |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ngs-bits               |misc                                                 |-                 |                                                                                     |                                                                        |
|freebayes              |variant calling                                      |01.2021 (1.3.3)   |test GATK on twin de-novo data: /mnt/share/evaluations/2020_07_29_twin_denovo/       |https://github.com/ekg/freebayes                                        |
|samtools               |BAM sorting                                          |01.2021 (1.11)    |                                                                                     |http://www.htslib.org/                                                  |
|vcflib                 |VCF normalization                                    |01.2021 (1.0.2)   |                                                                                     |https://github.com/vcflib/vcflib                                        |
|abra2                  |indel realignment                                    |01.2021 (2.23)    |                                                                                     |https://github.com/mozack/abra2                                         |
|samblaster             |duplicate removal                                    |01.2021 (0.1.26)  |                                                                                     |https://github.com/GregoryFaust/samblaster                              |
|bwa                    |mapping                                              |01.2021 (0.7.17)  |                                                                                     |https://github.com/lh3/bwa/                                             |
|bwa2                   |mapping                                              |01.2021 (2.1)     |                                                                                     |https://github.com/bwa-mem2/bwa-mem2                                    |
|ClinCNV                |CNV calling                                          |01.2021 (1.17.0)  |                                                                                     |https://github.com/imgag/ClinCNV                                        |
|VEP                    |annotation                                           |05.2021 (104.3)   | -> 105.0 available (PrimateAI)                                                      |https://github.com/Ensembl/ensembl-vep/releases                         |
|manta                  |structural variant calling (tumor/normal)            |01.2021 (1.6.0)   |                                                                                     |https://github.com/Illumina/manta                                       |
|InterOp                |reading InterOp metric files (Illumina NextSeq 1k/2k)|01.2021 (1.0.25)  | -> 1.1.16 available                                                                 |                                                                        |
|Circos                 |circos plot with CNVs,ROHS,etc                       |01.2021 (0.69-9)  |                                                                                     |http://circos.ca/software/download/                                     |
|ExpansionHunter        |Repeat expansion calling                             |01.2021 (4.0.2)   |                                                                                     |https://github.com/Illumina/ExpansionHunter                             |
|SpliceAI               |Predict splicing variant effect                      |01.2021 (1.3.1)   |                                                                                     |https://github.com/Illumina/SpliceAI                                    |



## Tools for somatic pipelines

|tool somatic           |used for                                                 |last update/check   |notes  |url                                                                     |
|-----------------------|---------------------------------------------------------|--------------------|-------|------------------------------------------------------------------------|
|strelka2               |variant calling (tumor/normal)                           |01.2021 (2.9.10)    |       |https://github.com/Illumina/strelka                                     |
|mantis                 |microsatelite instability (tumor/normal)                 |01.2021 (v1.0.5)    |       |https://github.com/OSU-SRLab/MANTIS/releases                            |
|varscan2               |variant calling                                          |01.2021 (2.4.4)     |       |https://github.com/dkoboldt/varscan                                     |


## Tools for RNA pipeline

|tool RNA               |used for                   |last update/check               |notes                                                                                               |url                                                                            |
|-----------------------|---------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
|STAR                   |mapping                    |08.2020 (2.7.3a)                |2.7.5a available                                                                                    |https://github.com/alexdobin/STAR                                              |
|subread                |read counting              |08.2020 (2.0.0)                 |2.0.1 available                                                                                     |http://subread.sourceforge.net/                                                |
|skewer                 |single-end adapter trimming|07.2020 (0.2.2)                 |                                                                                                    |https://github.com/relipmoc/skewer                                             |
|STAR-Fusion            |fusion detection           |07.2020 (1.9.0)                 |                                                                                                    |https://github.com/STAR-Fusion/STAR-Fusion                                     |
|samtools               |fusion detection           |10.2020 (1.7)                   |STAR-Fusion needs specific samtools version for certain RNA samples                                 |https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2|
