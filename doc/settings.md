# Settings

## Section 'tools-ngs'

Paths of general NGS tools.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'tools-longread'

Paths of tools used in long-read pipeline.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'tools-ngs-somatic'

Paths of tools used in somatic pipeline.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'tools-ngs-rna'

Paths of tools used in RNA pipeline.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'tools-primer'

Paths of tools used for primer design.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'singularity-container'

Paths of tools encapsulated in singularity containers.  
Paths are preset and should not be changed unless you want to exchange a container.

## Section 'folders'

This section contains path to important folders:

* `data_folder` - general data folder with downloaded database flat files.
* `local_data` - local data folder. This folder is used to copy database files for annotation from `data_folder` to a local temporary folder. The local temporary folder is usually on a SSD/NVME, so this will speed up the annotation.
* `test_data_folder` - folder for test data. This folder is needed for unit tests and tool tests during pipeline development.
* `project_folder` - folder that contains the NGS data. This folder has to be set for the NGSD integration (see below). Instead of providing a single path, it is also possible to provide project-type-specific paths - in this case the dictionary syntax of PHP is used, e.g. `project_folder['diagnostic'] = /mnt/projects/diagnostic/`

## Section 'general'

This section contains general settings:

* `copy_dbs_to_local_data` - Flag (true/false) that indicates if large database files (gnomAD, etc) are copied to the local folder.
* `delete_fastq_files` - Flag (true/false) that indicates if FASTQ files are to be deleted after the BAM/CRAM file is created after mapping.
* `cnv_bin_size_wgs` - Bin size used for CNV analysis of WGS samples.
* `cnv_bin_size_shallow_wgs` - Bin size used for CNV analysis of shallow WGS samples.
* `cnv_bin_size_longread_wgs` - Bin size used for CNV analysis of long-read WGS samples.
* `use_bwa1` - Flag (true/false) that indicates if BWA mem should be used for mapping instead of BWA mem 2 (BWA mem 2 is faster, but needs more RAM).
* `annotate_refseq_consequences` - Flag (true/false) that indicates if variant consequences based on RefSeq transcripts should be annotated in addition to variant consequences based on Ensembl transcripts.
* `custom_colums` - Used to add custom annotations to the output VCF/GSvar file. Each entry consists of a colon-speparated list of VCF file, INFO field name in the source VCF (prefixed with `CUSTOM_` in the annotated VCF) and column description. Provide the annotation using the dictionary syntax of PHP, e.g. `custom_colums['NGSD_counts'] = "/mnt/data/dbs/NGSD/NGSD_germline.vcf.gz;COUNTS;NGSD counts"`.

## Section 'dragen'

This section contains settings to run the germline/somatic data analysis on a on-site Illumina Dragen server:

* `dragen_in` - Folder into which FASTQ data is copied for Dragen data analysis.
* `dragen_out` - Folder into which Dragen output files are copied after the analysis.
* `dragen_data` - Folder used as working directory on the Dragen server.
* `dragen_genomes` - Folder containing genome files for Dragen on the Dragen server.
* `dragen_log` - Folder used to write SGE logs of Dragen analyses.
* `queues_dragen` - Comma-separated list of SGE queues for Dragen servers (one per server).
* `use_dragen_sv_calling` - Flag (true/false) that indicates if SV calls produced by Dragen are used. If set to `false`, SV calling based on Manta is used. 

## Section 'mysql-databases'

This section contains NGSD data (host, database name, credentials).

## Section 'grid_engine'

This section contains settings for automated data analysis via SGE:

* `queues_default` - Default SGE queues.
* `queues_research` - Research SGE queues.
* `queues_high_priority` - High-priority queues.
* `queues_high_mem` - High-memory consumption queues (RNA).
