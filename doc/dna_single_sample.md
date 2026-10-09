# megSAP - DNA analysis (single sample)

### Basics

Single sample DNA analysis is performed using the `analyze.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/analyze.php --help

The main parameters that you have to provide are:

* `folder` - The sample folder, which contains the raw data FASTQ files or ORA file.
* `name` - The sample name, which must be a prefix of the FASTQ files.
* `steps` -  Analysis steps to perform. For example `ma,vc,sv` to perform mapping, small variant calling and structural variant calling (with annotation).
* `system` - The [processing system INI file](processing_system_ini_file.md).

### Running an analysis

The analysis pipeline assumes that that all raw data to analyze (FASTQ, ORA or BAM) is in the sample folder. If that is the case, the whole analysis is performed with one command, for example like this:

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system twist_exome.ini -steps ma,vc,sv,re

In the example above, the configuration of the pipeline is done using the `twist_exome.ini` file, which contains all necessary information (see [processing system INI file](processing_system_ini_file.md)).

### Running an analysis based on Illumina DRAGEN output

A short instruction how to setup the DRAGEN can be found [here](setup_dragen.md).
To run an analysis with DRAGEN you first run the DRAGEN analysis using `analyze_dragen.php`:

	php megSAP/src/Pipelines/analyze_dragen.php -folder Sample_NA12878_01 -name NA12878_01 -system twist_exome.ini -steps ma -no_queuing

Next, you run the main megSAP analysis using `analysis.php`: 

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system twist_exome.ini -steps vc,sv

Note: If you have a queing system configured for megSAP, you can also skip the parameter `-no_queuing` to automaticall queue the main megSAP analysis when the DRAGEN analysis is done.

### Tools used in this analysis pipeline

The following tools are used for mapping (step `ma`) and small variants calling/annotation (step `vc`):

| task                                           | tool                 | comments                                         |
|------------------------------------------------|----------------------|--------------------------------------------------|
| mapping - adapter and quality trimming         | SeqPurge             | Skipped in in DRAGEN mode.                       |
| mapping - mapping and alignment                | bwa-mem2             | Performed by DRAGEN in DRAGEN mode.              |
| mapping - duplicate marking                    | samblaster           | Performed by DRAGEN in DRAGEN mode.              |
| variant calling - calling of SNVs and InDels   | DeepVariant          | Performed by DRAGEN in DRAGEN mode.              |
| variant calling - decompose complex variants   | vcfallelicprimitives | Performed by DRAGEN in DRAGEN mode.              |
| variant calling - break multi-allelic variants | VcfBreakMulti        |                                                  |
| variant calling - left-normalization of InDels | VcfLeftNormalize     |                                                  |
| annotation - general                           | VcfAnnotateFrom*     |                                                  |

CNV calling/annotation (step `cn`) is performed using these tools:

| task                                               | tool                 | comments                                            |
|----------------------------------------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              |                                                     |
| annotation - general                               | BedAnnotateFromBed   | Several data sources are annotated using this tool. |
| annotation - gene information                      | CnvGeneAnnotation    |                                                     |
| annotation - overlapping pathogenic CNVs from NGSD | NGSDAnnotateCNV      |                                                     |

SV calling/annotation (step `sv`) is performed using these tools:

| task                                      | tool                            | comments                                         |
|-------------------------------------------|---------------------------------|--------------------------------------------------|
| SV calling                                | Manta                           | Performed by DRAGEN in DRAGEN mode.              |
| annotation - gene information             | BedpeGeneAnnotation             |                                                  |
| annotation - matching SVs from NGSD       | BedpeAnnotateCounts             |                                                  |
| annotation - breakpoint density from NGSD | BedpeAnnotateBreakpointDensity  |                                                  |

RE calling (step `re`) is performed using these tools:

| task                                      | tool                            | comments                                            |
|-------------------------------------------|---------------------------------|-----------------------------------------------------|
| RE calling                                | ExpansionHunter                 |                                                     |


A complete list of all tools and databases used in megSAP including version and when they were last updated can be found [here](development/update_overview.md).

### Performance

Performance benchmarks of the the megSAP pipeline can be found [here](performance.md).

### Output

After the analysis with the steps `ma,vc,sv`, these files are created in the output folder:

1. mapped reads in CRAM format  
2. small variant list in VCF format
3. small variant list in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md)
4. structural variant list in VCF format
5. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which can be opened with with GSvar or Firefox

### Test data

Example data which can be analyzed using the command above can be downloaded from [here](https://megsap.de/download/examples/NA12878_01.zip).

[back to the start page](../README.md)



