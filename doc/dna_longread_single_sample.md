# megSAP - DNA analysis (single sample)

### Basics

Single sample long-read DNA analysis is performed using the `analyze_longread.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/analyze_longread.php --help

The main parameters that you have to provide are:

* `folder` - The sample folder, which contains the the FASTQ files as produced by bcl2fastq2.
* `name` - The sample name, which must be a prefix of the FASTQ files.
* `steps` -  Analysis steps to perform. Please use `ma,vc,sv,re,an,db` to perform mapping, variant calling (SNV/InDels, CNVs, SVs and REs), annotation and database import.
* `system` - The [processing system INI file](processing_system_ini_file.md).

### Running an analysis

The analysis pipeline assumes that that all data to analyze (FastQ files) resides in a sample folder. If that is the case, the whole analysis is performed with one command, for example like this:

	php megSAP/src/Pipelines/analyze_longread.php -folder Sample_NA12878_01 -name NA12878_01 -system SQK-114.ini -steps ma,vc,cn,sv,re,an

In the example above, the configuration of the pipeline is done using the `SQK-114.ini` file, which contains all necessary information (see [processing system INI file](processing_system_ini_file.md)).



### Tools used in this analysis pipeline

The following tools are used for mapping and calling of small variants and annotation of small variants:

| step                                           | tool                 | comments  |
|------------------------------------------------|----------------------|-----------|
| mapping                                        | minimap2             |           |
| mapping - extract methylation info             | modkit               |           |
| variant calling - calling of SNVs and InDels   | clair3               |           |
| variant calling - decompose complex variants   | vcfallelicprimitives |           |
| variant calling - break multi-allelic variants | vcfbreakmulti        |           |
| variant calling - left-normalization of InDels | VcfLeftNormalize     |           |
| annotation                                     | VEP                  |           |

CNV calling and annotation is performed using these tools:

| step                                               | tool                 | comments                                            |
|----------------------------------------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              |                                                     |
| annotation - general                               | BedAnnotateFromBed   | Several data sources are annotated using this tool. |
| annotation - gene information                      | CnvGeneAnnotation    |                                                     |
| annotation - overlapping pathogenic CNVs from NGSD | NGSDAnnotateCNV      |                                                     |

SV calling and annotation is performed using these tools:

| step                                      | tool                            | comments   |
|-------------------------------------------|---------------------------------|------------|
| SV calling                                | Sniffles                        |            |
| annotation - gene information             | BedpeGeneAnnotation             |            |
| annotation - matching SVs from NGSD       | BedpeAnnotateCounts             |            |
| annotation - breakpoint density from NGSD | BedpeAnnotateBreakpointDensity  |            |

RE calling using these tools:

| step                                      | tool                            | comments                                            |
|-------------------------------------------|---------------------------------|-----------------------------------------------------|
| RE calling                                | Straglr                         |                                                     |

Phasing is performed using these tools:

| step                                      | tool                            | comments   |
|-------------------------------------------|---------------------------------|------------|
| phasing                                   | longphase                       |            |



A complete list of all tools and databases used in megSAP including version and when they were last updated can be found [here](update_overview.md).

### Performance

Performance benchmarks of the the megSAP pipeline can be found [here](performance.md).

### Output

After the analysis, these files are created in the output folder:

1. mapped reads in BAM format  
2. a variant list in VCF format
3. a variant list in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md)
4. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which can be opened with a web browser


