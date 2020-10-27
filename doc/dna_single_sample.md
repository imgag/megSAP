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

### Poster

A poster about megSAP which describes the all steps of the single-sample analysis pipeline can be found 
[here](Poster_April_2017.pdf).

### Running an analysis

The analysis pipeline assumes that that all data to analyze resides in a sample folder as produced by Illumina's [bcl2fastq](http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) tool. If that is the case, the whole analysis is performed with one command, for example like this:

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system hpHBOCv5.ini -steps ma,vc

In the example above, the configuration of the pipeline is done using the `hpHBOCv5.ini` file, which contains all necessary information (see [processing system INI file](processing_system_ini_file.md)).

### Running an analysis with DRAGEN

A short instruction how to setup the DRAGEN can be found [here](setup_dragen.md).
To run an analysis with DRAGEN mapping you simply has to pass the parameter `-use_dragen` to the `analysis.php`: 

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system hpHBOCv5.ini -steps ma,vc -use_dragen

### Performance

A performance comparison between the megSAP pipeline using DRAGEN vs. bwa can be found [here](performance.md)

### Test data

Example data which can be analyzed using the command above can be downloaded from [here](https://download.imgag.de/NA12878_01.zip).

### Output

After the analysis, these files are created in the output folder:

1. mapped reads in BAM format  
2. a variant list in VCF format
3. a variant list in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md)
4. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which can be opened with a web browser

[back to the start page](../README.md)








