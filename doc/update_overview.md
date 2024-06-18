# Tool and database overview

## Databases general

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ClinGen                |pipeline                |12.2023 (latest - not versioned)|                                                                                                    |https://ftp.clinicalgenome.org/                                         |
|RepeatMasker           |pipeline                |12.2023 (140131)                |                                                                                                    |http://www.repeatmasker.org/species/hg.html                             |
|ClinVar (SNVs and CNVs)|pipeline, NGSDImportHPO |06.2024 (20240603)              |update IGV custom tacks in GSvar                                                                    |https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/   |
|HGNC                   |pipeline, NGSDImportHGNC|06.2024 (latest - not versioned)|                                                                                                    |https://ftp.ebi.ac.uk/pub/databases/genenames/                          |
|gnomAD (genome)        |pipeline                |12.2023 (3.1.2)                 |version v4.1 available > test and update                                                            |http://gnomad.broadinstitute.org/downloads                              |
|gnomAD (constraints)   |NGSDImportGeneInfo      |06.2024 (2.1.1)                 |v4 still beta https://gnomad.broadinstitute.org/news/2024-03-gnomad-v4-0-gene-constraint/           |http://gnomad.broadinstitute.org/downloads                              |
|phyloP                 |pipeline                |12.2023 (05.2015)               |                                                                                                    |http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/            |
|CADD                   |pipeline                |12.2023 (1.6)                   |                                                                                                    |http://cadd.gs.washington.edu/download                                  |
|REVEL                  |pipeline                |12.2023 (1.3)                   |                                                                                                    |https://sites.google.com/site/revelgenomics/downloads                   |
|AlphaMissense          |pipeline                |12.2023 (03.08.2023)            |                                                                                                    |https://console.cloud.google.com/storage/browser/dm_alphamissense       |
|OMIM                   |pipeline, NGSDImportHPO |07.2023 (latest - not versioned)|                                                                                                    |https://omim.org/downloads/                                             |
|HGMD (SNVs and CNVs)   |pipeline                |06.2024 (2024.1)                |update IGV custom tacks in GSvar                                                                    |https://apps.ingenuity.com/ingsso/login                                 |
|Ensembl                |NGSDImportEnsembl       |06.2024 (GRCh38.112)            |update IGV genome                                                                                   |https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/              |
|HPO                    |NGSDImportHPO           |06.2024 (2024-04-26)            |send updated HPO list to Anne (scripts/2021_10_21_hpo_update/)                                      |https://github.com/obophenotype/human-phenotype-ontology                |
|GenCC                  |NGSDImportHPO           |06.2024 (latest - not versioned)|                                                                                                    |https://search.thegencc.org/download                                    |
|DECIPER                |NGSDImportHPO           |12.2023 (14_6_2024)             |                                                                                                    |https://www.deciphergenomics.org/about/downloads/data                   |
|ORPHA                  |NGSDImportORPHA         |06.2024 (latest - not versioned)|Products 1 and 6 are updated tweice a year only (July and December)                                 |https://github.com/Orphanet                                             |

## Databases for somatic pipelines

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|CancerHotspots         |pipeline (somatic)      |01.2024                         |version is final and does not change                                                                |https://www.cancerhotspots.org                                          |
|NCG7.1                 |pipeline (somatic)      |01.2024 (v7.1)                  |					                                                                                   |http://ncg.kcl.ac.uk/                                                   |
|COSMIC CMC             |pipeline (somatic)      |01.2024 (v98)                   |                                                                                                    |https://cancer.sanger.ac.uk/cmc                                         |
|Human Protein Atlas    |pipeline (somatic)      |01.2024 (v23)                   |                                                                                                    |https://www.proteinatlas.org/about/download                             |

## Databases for RNA pipelines

|database               |used for                |last update/check               |notes                                                                                               |url                                                                     |
|-----------------------|------------------------|--------------------------------|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|Ensembl GTF file       |pipeline (RNA)          |03.2023 (GRCh38.109)            |                                                                                                    |https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/               |

## Tools general

|tool                   |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|ngs-bits               |annotation, quality control, ...                     |12.2023 (2023_11) |                                                                                     |                                                                        |
|samtools               |BAM sorting                                          |12.2023 (1.19)    |                                                                                     |http://www.htslib.org/                                                  |
|bwa2                   |mapping (default)                                    |12.2023 (2.2.1)   |                                                                                     |https://github.com/bwa-mem2/bwa-mem2                                    |
|bwa                    |mapping (if use_bwa1 is true in settings)            |12.2023 (0.7.17)  |                                                                                     |https://github.com/lh3/bwa/                                             |
|samblaster             |duplicate removal                                    |12.2023 (0.1.26)  |                                                                                     |https://github.com/GregoryFaust/samblaster                              |
|abra2                  |indel realignment                                    |12.2023 (2.23)    | update to 2.24 not possible: https://github.com/mozack/abra2/issues/46              |https://github.com/mozack/abra2                                         |
|freebayes              |variant calling                                      |12.2023 (1.3.6)   | update to 1.3.7 not possible: https://github.com/freebayes/freebayes/issues/765     |https://github.com/ekg/freebayes                                        |
|vcflib                 |VCF normalization                                    |12.2023 (1.0.3)   | update to 1.0.9 not possible: https://github.com/vcflib/vcflib/issues/398           |https://github.com/vcflib/vcflib                                        |
|ClinCNV                |CNV calling                                          |12.2023 (1.18.3)  |                                                                                     |https://github.com/imgag/ClinCNV                                        |
|VEP                    |annotation                                           |12.2023 (110.1)   |                                                                                     |https://github.com/Ensembl/ensembl-vep/releases                         |
|manta                  |structural variant calling                           |12.2023 (1.6.0)   |                                                                                     |https://github.com/Illumina/manta                                       |
|InterOp                |reading InterOp metric files (Illumina NextSeq 1k/2k)|12.2023 (1.2.4)   | version 1.3.0 available but update not necessary                                    |https://github.com/Illumina/interop                                     |
|Circos                 |circos plot with CNVs,ROHS,etc                       |12.2023 (0.69.9)  |                                                                                     |http://circos.ca/software/download/                                     |
|ExpansionHunter        |Repeat expansion calling                             |12.2023 (5.0.0)   |                                                                                     |https://github.com/Illumina/ExpansionHunter                             |
|SpliceAI               |Predict splicing variant effect                      |12.2023 (1.3.1)   |                                                                                     |https://github.com/Illumina/SpliceAI                                    |
|REViewer               |Repeat expansion visualization                       |12.2023 (0.2.7)   |                                                                                     |https://github.com/Illumina/REViewer                                    |
|bedtools               |Masking false duplications in genome                 |12.2023 (2.31.1)  |                                                                                     |https://github.com/arq5x/bedtools2/releases/                            |
|ORAD                   |Illumina ORA file decompression                      |12.2023 (2.6.1)   |                                                                                     |                                                                        |

## Tools for somatic pipelines

|tool somatic           |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|strelka2               |variant calling (tumor/normal)                       |01.2024 (2.9.10)  |                                                                                     |https://github.com/Illumina/strelka                                     |
|mantis                 |microsatelite instability (tumor/normal)             |01.2024 (v1.0.5)  |                                                                                     |https://github.com/OSU-SRLab/MANTIS/releases                            |
|varscan2               |variant calling                                      |01.2024 (2.4.6)   |                                                                                     |https://github.com/dkoboldt/varscan                                     |
|umiVar2                |variant calling cfDNA                                |01.2024 (2023_11) |                                                                                     |https://github.com/imgag/umiVar2                                        |
|hla-genotyper          |HLA genotyping                                       |01.2024 (2022_05) |                                                                                     |https://github.com/axelgschwind/hla-genotyper                           |
|SigProfilerExtractor   |mutational signatures                                |01.2024 (1.1.23)  |                                                                                     |https://github.com/AlexandrovLab/SigProfilerExtractor                   |


## Tools for RNA pipeline

|tool RNA               |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|STAR                   |mapping                                              |03.2023 (2.7.10b) |                                                                                     |https://github.com/alexdobin/STAR                                       |
|subread                |read counting                                        |03.2023 (2.0.4)   |                                                                                     |http://subread.sourceforge.net/                                         |
|Arriba                 |fusion detection                                     |03.2023 (2.4.0)   |                                                                                     |https://github.com/suhrig/arriba                                        |
|Kraken2                |fastq filtering                                      |03.2023 (2.1.2)   |                                                                                     |https://github.com/DerrickWood/kraken2                                  |
|umi_tools              |UMI extraction                                       |03.2023 (1.1.4)   |                                                                                     |https://github.com/CGATOxford/UMI-tools                                 |

## Tools for longread pipeline

|tool longread          |used for                                             |last update/check |notes                                                                                |url                                                                     |
|-----------------------|-----------------------------------------------------|------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------|
|minimap2               |mapping                                              |01.2024 (2.26)    |                                                                                     |https://github.com/lh3/minimap2                                         |
|clair3                 |small variant calling                                |01.2024 (1.0.5)   |                                                                                     |https://github.com/HKU-BAL/Clair3                                       |
|whatshap               |phasing (helper tool for clair3)                     |01.2024 (2.1)     | installed via pip in megSAP python environment                                      |https://github.com/whatshap/whatshap                                    |
|pypy3                  |alt. python implementation (helper tool for clair3)  |01.2024 (v7.3.14) |                                                                                     |https://www.pypy.org/                                                   |
|parallel				|tool to execute in parallel (helper tool for clair3) |06.2023 (20230522)| website down/ certificate expired                                                   |https://www.gnu.org/software/parallel/                                  |
|longphase              |phasing                                              |01.2024 (v1.6)    |                                                                                     |https://github.com/twolinin/longphase                                   |
|sniffles               |structural variant calling                           |01.2024 (v2.2)    | installed via pip in megSAP python environment                                      |https://github.com/fritzsedlazeck/Sniffles                              |
