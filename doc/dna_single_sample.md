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

### Configure mapping using DRAGEN

The megSAP pipeline also supports mapping using the illumina DRAGEN server. 
For that you have to install and configure `php` and `Sun GridEngine` on the DRAGEN server. Additionally some settings in the `settings.ini` of megSAP needs to be configured:

* `dragen_user` - User which is used to run the analysis on the DRAGEN server. It has to be the same user who started the complete analysis and has to have read and write access to the folders defined below.

* `dragen_in`/`dragen_out` - Transfer folders which have to be accessible from both the server which performs the analysis and the DRAGEN server. These two foloders are used to transfer data to the DRAGEN server (e. g. FastQ files) and transfer data from the DRAGEN server to the analysis server (e. g. BAM files).

* `dragen_data` - Temporary folder on the DRAGEN server in which the mapping is performed. This folder should be located on the fast SSD storage of the DRAGEN server (usually: `/staging/...`) and is created for each mapping and deleted after the mapped data has been moved to the tranfer folder.

* `dragen_genomes` - Path to the genome reference hash tables. Should also be stored on the DRAGEN SSD storage.

* `queues_dragen` - Queue(s) with 1 slot where the DRAGEN mapping jobs are submitted to and are run on the DRAGEN server(s).

### Running an analysis with DRAGEN

After megSAP is configured correctly using the DRAGEN server you can perform a analysis using the DRAGEN mapping by passing the parameter `-use_dragen` to the `analysis.php`: 

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system hpHBOCv5.ini -steps ma,vc -use_dragen

### Performance

A performance comparison between the megSAP pipeline using DRAGEN vs. bwa can be found [here](performance.md)

### Test data

Example data which can be analyzed using the command above can be downloaded from [here](https://download.imgag.de/NA12878_01.zip).

### Output

After the analysis, these files are created in the output folder:

1. mapped reads in BAM format  
2. a variant list in VCF format
3. a variant list in [GSvar format](gsvar_format.md)
4. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which can be opened with a web browser

[back to the start page](../README.md)








