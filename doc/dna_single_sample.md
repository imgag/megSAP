# megSAP - DNA analysis (single sample)

### Basics

Single sample DNA analysis is performed using the `analyze.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/analyze.php --help

The main parameters that you have to provide are:

* `folder` - The sample folder, which contains the the FASTQ files as produced by bcl2fastq2.
* `name` - The sample name, which must be a prefix of the FASTQ files.
* `steps` -  Analysis steps to perform. Please use `ma,vc` to perform mapping and variant calling (with annotation).
* `system` - The [processing system INI file](processing_system_ini_file.md).

### Running an analysis

The analysis pipeline assumes that that all data to analyze resides in a sample folder as produced by Illumina's [bcl2fastq](http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) tool. If that is the case, the whole analysis is performed with one command, for example like this:

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system twist_exome.ini -steps ma,vc

In the example above, the configuration of the pipeline is done using the `twist_exome.ini` file, which contains all necessary information (see [processing system INI file](processing_system_ini_file.md)).

### Running an analysis with DRAGEN

A short instruction how to setup the DRAGEN can be found [here](setup_dragen.md).
To run an analysis with DRAGEN mapping you simply has to pass the parameter `-use_dragen` to the `analysis.php`: 

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system twist_exome.ini -steps ma,vc -use_dragen

### Tools used in this analysis pipeline

The following tools are used for mapping and calling of small variants and annotation of small variants:

| step                                           | tool                 | comments                                         |
|------------------------------------------------|----------------------|--------------------------------------------------|
| mapping - adapter and quality trimming         | SeqPurge             |                                                  |
| mapping - mapping and alignment                | bwa-mem2             | Performed by Dragen if '-use_dragen' is used.    |
| mapping - duplicate marking                    | samblaster           | Performed by Dragen if '-use_dragen' is used.    |
| mapping - indel realignment                    | ABRA2                |                                                  |
| variant calling - calling of SNVs and InDels   | freebayes            |                                                  |
| variant calling - decompose complex variants   | vcfallelicprimitives | Performed by Dragen if '-use_dragen' is used.    |
| variant calling - break multi-allelic variants | vcfbreakmulti        |                                                  |
| variant calling - left-normalization of InDels | VcfLeftNormalize     |                                                  |
| annotation                                     | VEP                  |                                                  |

CNV calling and annotation is performed using these tools:

| step                                               | tool                 | comments                                            |
|----------------------------------------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              |                                                     |
| annotation - general                               | BedAnnotateFromBed   | Several data sources are annotated using this tool. |
| annotation - gene information                      | CnvGeneAnnotation    |                                                     |
| annotation - overlapping pathogenic CNVs from NGSD | NGSDAnnotateCNV      |                                                     |

SV calling and annotation is performed using these tools:

| step                                      | tool                            | comments                                                                                                  |
|-------------------------------------------|---------------------------------|-----------------------------------------------------------------------------------------------------------|
| SV calling                                | Manta                           | Performed by Dragen if '-use_dragen' is used and 'dragen_sv_calling' is enabled in the 'settings.ini'.    |
| annotation - gene information             | BedpeGeneAnnotation             |                                                                                                           |
| annotation - matching SVs from NGSD       | BedpeAnnotateCounts             |                                                                                                           |
| annotation - breakpoint density from NGSD | BedpeAnnotateBreakpointDensity  |                                                                                                           |

RE calling using these tools:

| step                                      | tool                            | comments                                            |
|-------------------------------------------|---------------------------------|-----------------------------------------------------|
| RE calling                                | ExpansionHunter                 |                                                     |


A complete list of all tools and databases used in megSAP including version and when they were last updated can be found [here](update_overview.md).

### Performance

Performance benchmarks of the the megSAP pipeline can be found [here](performance.md).

### Output

After the analysis, these files are created in the output folder:

1. mapped reads in BAM format  
2. a variant list in VCF format
3. a variant list in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md)
4. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which can be opened with a web browser

### Test data

Example data which can be analyzed using the command above can be downloaded from [here](https://download.imgag.de/NA12878_01.zip).

[back to the start page](../README.md)
