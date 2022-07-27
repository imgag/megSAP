# Tool and database overview

## Databases general

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


## Tools general

|tool                   |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ngs-bits               |annotation, quality control, ...                     |-                 |                                                                                     |                                                                        |
|samtools               |BAM sorting                                          |07.2022 (1.15.1)  |                                                                                     |http://www.htslib.org/                                                  |
|bwa2                   |mapping (default)                                    |07.2022 (2.2.1)   |                                                                                     |https://github.com/bwa-mem2/bwa-mem2                                    |
|bwa                    |mapping (if use_bwa1 is true in settings)            |07.2022 (0.7.17)  |                                                                                     |https://github.com/lh3/bwa/                                             |
|samblaster             |duplicate removal                                    |07.2022 (0.1.26)  |                                                                                     |https://github.com/GregoryFaust/samblaster                              |
|abra2                  |indel realignment                                    |07.2022 (2.23)    | update to 2.24 not possible because of https://github.com/mozack/abra2/issues/46    |https://github.com/mozack/abra2                                         |
|freebayes              |variant calling                                      |07.2022 (1.3.6)   |                                                                                     |https://github.com/ekg/freebayes                                        |
|vcflib                 |VCF normalization                                    |07.2022 (1.0.3)   |                                                                                     |https://github.com/vcflib/vcflib                                        |
|VEP                    |annotation                                           |01.2021 (104.4)   | update to 107 when VcfAnnotateconsequences is used for consequence annotation       |https://github.com/Ensembl/ensembl-vep/releases                         |
|SpliceAI               |Predict splicing variant effect                      |07.2022 (1.3.1)   |                                                                                     |https://github.com/Illumina/SpliceAI                                    |
|ClinCNV                |CNV calling                                          |07.2022 (1.17.2)  |                                                                                     |https://github.com/imgag/ClinCNV                                        |
|manta                  |structural variant calling                           |07.2022 (1.6.0)   |                                                                                     |https://github.com/Illumina/manta                                       |
|ExpansionHunter        |Repeat expansion calling                             |07.2022 (5.0.0)   |                                                                                     |https://github.com/Illumina/ExpansionHunter                             |
|REViewer               |Repeat expansion visualization                       |07.2022 (0.2.7)   |                                                                                     |https://github.com/Illumina/REViewer                                    |
|Circos                 |circos plot with CNVs,ROHS,etc                       |07.2022 (0.69-9)  |                                                                                     |http://circos.ca/software/download/                                     |
|InterOp                |reading InterOp metric files (Illumina NextSeq 1k/2k)|07.2022 (1.0.25)  | v1.1.25 available, but update is not necessary as the QC import works               |                                                                        |


## Tools for somatic pipelines

|tool somatic           |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|strelka2               |variant calling (tumor/normal)                       |01.2021 (2.9.10)  |                                                                                     |https://github.com/Illumina/strelka                                     |
|mantis                 |microsatelite instability (tumor/normal)             |01.2021 (v1.0.5)  |                                                                                     |https://github.com/OSU-SRLab/MANTIS/releases                            |
|varscan2               |variant calling                                      |01.2021 (2.4.4)   |                                                                                     |https://github.com/dkoboldt/varscan                                     |


## Tools for RNA pipeline

|tool RNA               |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|STAR                   |mapping                                              |08.2020 (2.7.3a)  |2.7.5a available                                                                     |https://github.com/alexdobin/STAR                                       |
|subread                |read counting                                        |08.2020 (2.0.0)   |2.0.1 available                                                                      |http://subread.sourceforge.net/                                         |
|Arriba                 |fusion detection                                     |01.2022 (2.2.1)   |2.3.0 available                                                                      |https://github.com/suhrig/arriba                                        |
|Kraken2                |???                                                  |07.2022 (2.1.2)   |                                                                                     |https://github.com/DerrickWood/kraken2                                  |
