# Tool and database overview

## Databases general

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ClinGen                |pipeline                |07.2022 (latest - not versioned)|                                                                                                    |https://ftp.clinicalgenome.org/                                         |
|RepeatMasker           |pipeline                |07.2022 (140131)                |                                                                                                    |http://www.repeatmasker.org/species/hg.html                             |
|ClinVar (SNVs and CNVs)|pipeline, NGSDImportHPO |10.2022 (20220924)              |update GSvar IGV file                                                                               |https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2021/   |
|HGNC                   |pipeline, NGSDImportHGNC|07.2022 (latest - not versioned)|                                                                                                    |https://ftp.ebi.ac.uk/pub/databases/genenames/                          |
|gnomAD (genome)        |pipeline                |07.2022 (3.1.2)                 |                                                                                                    |http://gnomad.broadinstitute.org/downloads                              |
|gnomAD (constraints)   |NGSDImportGeneInfo      |07.2022 (2.1.1)                 |                                                                                                    |http://gnomad.broadinstitute.org/downloads                              |
|phyloP - for VEP       |pipeline                |07.2022 (05.2015)               |                                                                                                    |http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/            |
|CADD - for VEP         |pipeline                |07.2022 (1.6)                   |                                                                                                    |http://cadd.gs.washington.edu/download                                  |
|REVEL - for VEP        |pipeline                |07.2022 (1.3)                   |                                                                                                    |https://sites.google.com/site/revelgenomics/downloads                   |
|OMIM                   |pipeline, NGSDImportHPO |10.2022 (latest - not versioned)|                                                                                                    |https://omim.org/downloads/                                             |
|HGMD (SNVs and CNVs)   |pipeline                |10.2022 (2022.3)                |update GSvar IGV file                                                                               |https://apps.ingenuity.com/ingsso/login                                 |
|Ensembl                |NGSDImportEnsembl       |07.2022 (GRCh38.107)            |                                                                                                    |https://ftp.ensembl.org/pub/release-107/gff3/homo_sapiens/              |
|HPO                    |NGSDImportHPO           |07.2022 (latest - not versioned)|send updated HPO list to Anne (/mnt/users/ahsturm1/Sandbox/2021_10_21_hpo_update/)                  |https://hpo.jax.org/app/                                                |
|GenCC                  |NGSDImportHPO           |07.2022 (latest - not versioned)|                                                                                                    |https://search.thegencc.org/download                                    |
|DECIPER                |NGSDImportHPO           |07.2022 (22_12_2021)            |                                                                                                    |https://www.deciphergenomics.org/about/downloads/data                   |
|ORPHA                  |NGSDImportORPHA         |07.2022 (latest - not versioned)|                                                                                                    |https://github.com/Orphanet/Orphadata.org/                              |


## Databases for somatic pipelines

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|CancerHotspots         |pipeline (somatic)      |01.2021                         |version is final and does not change                                                                |https://www.cancerhotspots.org                                          |
|NCG6.0                 |pipeline (somatic)      |01.2021                         |version 7 available                                                                                 |http://ncg.kcl.ac.uk/                                                   |
|COSMIC CMC             |pipeline (somatic)      |01.2021                         |                                                                                                    |https://cancer.sanger.ac.uk/cmc                                         |
|Human Protein Atlas    |pipeline (somatic)      |08.2022 (v21.1)                 |                                                                                                    |https://www.proteinatlas.org/about/download                             |

## Databases for RNA pipelines

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|Ensembl GTF file       |pipeline (RNA)  	     |08.2022 (release 107)           |                                                                                                    |http://ftp.ensembl.org/pub/                                             |

## Tools general

|tool                   |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ngs-bits               |annotation, quality control, ...                     |10.2022 (2022_10) |                                                                                     |                                                                        |
|samtools               |BAM sorting                                          |07.2022 (1.15.1)  |                                                                                     |http://www.htslib.org/                                                  |
|bwa2                   |mapping (default)                                    |07.2022 (2.2.1)   |                                                                                     |https://github.com/bwa-mem2/bwa-mem2                                    |
|bwa                    |mapping (if use_bwa1 is true in settings)            |07.2022 (0.7.17)  |                                                                                     |https://github.com/lh3/bwa/                                             |
|samblaster             |duplicate removal                                    |07.2022 (0.1.26)  |                                                                                     |https://github.com/GregoryFaust/samblaster                              |
|abra2                  |indel realignment                                    |07.2022 (2.23)    | update to 2.24 not possible because of https://github.com/mozack/abra2/issues/46    |https://github.com/mozack/abra2                                         |
|freebayes              |variant calling                                      |07.2022 (1.3.6)   |                                                                                     |https://github.com/ekg/freebayes                                        |
|vcflib                 |VCF normalization                                    |07.2022 (1.0.3)   |                                                                                     |https://github.com/vcflib/vcflib                                        |
|VEP                    |annotation                                           |08.2022 (107.0)   |                                                                                     |https://github.com/Ensembl/ensembl-vep/releases                         |
|SpliceAI               |Predict splicing variant effect                      |07.2022 (1.3.1)   |                                                                                     |https://github.com/Illumina/SpliceAI                                    |
|ClinCNV                |CNV calling                                          |10.2022 (1.18.0)  |                                                                                     |https://github.com/imgag/ClinCNV                                        |
|manta                  |structural variant calling                           |07.2022 (1.6.0)   |                                                                                     |https://github.com/Illumina/manta                                       |
|ExpansionHunter        |Repeat expansion calling                             |07.2022 (5.0.0)   |                                                                                     |https://github.com/Illumina/ExpansionHunter                             |
|REViewer               |Repeat expansion visualization                       |07.2022 (0.2.7)   |                                                                                     |https://github.com/Illumina/REViewer                                    |
|Circos                 |circos plot with CNVs,ROHS,etc                       |07.2022 (0.69-9)  |                                                                                     |http://circos.ca/software/download/                                     |
|InterOp                |reading InterOp metric files (Illumina NextSeq 1k/2k)|07.2022 (1.0.25)  | v1.1.25 available, but update is not necessary as the QC import works               |                                                                        |


## Tools for somatic pipelines

|tool somatic           |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|strelka2               |variant calling (tumor/normal)                       |07.2022 (2.9.10)  |                                                                                     |https://github.com/Illumina/strelka                                     |
|mantis                 |microsatelite instability (tumor/normal)             |07.2022 (v1.0.5)  |                                                                                     |https://github.com/OSU-SRLab/MANTIS/releases                            |
|varscan2               |variant calling                                      |07.2022 (2.4.4)   |                                                                                     |https://github.com/dkoboldt/varscan                                     |
|umiVar2                |variant calling cfDNA                                |07.2022 (latest)  |                                                                                     |https://github.com/dkoboldt/varscan                                     |


## Tools for RNA pipeline

|tool RNA               |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|STAR                   |mapping                                              |08.2022 (2.7.10a) |                                                                                     |https://github.com/alexdobin/STAR                                       |
|subread                |read counting                                        |08.2022 (2.0.3)   |                                                                                     |http://subread.sourceforge.net/                                         |
|Arriba                 |fusion detection                                     |08.2022 (2.2.1)   |2.3.0 available                                                                      |https://github.com/suhrig/arriba                                        |
|Kraken2                |fastq filtering (?)                                  |07.2022 (2.1.2)   |                                                                                     |https://github.com/DerrickWood/kraken2                                  |
|hla-genotyper          |HLA genotyping                                       |07.2022 (2022_05) |                                                                                     |https://github.com/axelgschwind/hla-genotyper                           |
|umi_tools              |UMI extraction                                       |08.2022 (1.1.2)   |                                                                                     |https://github.com/CGATOxford/UMI-tools                                 |
