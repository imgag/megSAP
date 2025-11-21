# Settings

## Section 'tools'

Paths of locally installed tools.  
Paths are preset and should not be changed unless you want to exchange a tool.

## Section 'tools-data'

Paths of downloaded data used by specific tools.
Paths are preset and should not be changed unless you want to change where the data is stored.
* `ngs-bits_local` - the path to a local ngs-bits installation can be provided here. If no path is provided the ngs-bits apptainer container listed in the 'apptainer-container' section is used.

## Section 'apptainer-container'

Paths of tools encapsulated in apptainer containers.  
Paths are preset and should not be changed unless you want to exchange a container.

## Section 'folders'

This section contains path to important folders:

* `data_folder` - general data folder with downloaded database flat files.
* `local_data` - local data folder. This folder is used to copy database files for annotation from `data_folder` to a local temporary folder. The local temporary folder is usually on a SSD/NVME, so this will speed up the annotation.
* `container_folder` - general container folder with downloaded apptainer container images. Images can be copied to a subdirectory of `locla_data` to improve access speed.
* `test_data_folder` - folder for test data. This folder is needed for unit tests and tool tests during pipeline development.
* `project_folder` - folder that contains the NGS data. This folder has to be set for the NGSD integration (see below). Instead of providing a single path, it is also possible to provide project-type-specific paths - in this case the dictionary syntax of PHP is used, e.g. `project_folder['diagnostic'] = /mnt/projects/diagnostic/`

## Section 'general'

This section contains general settings:

* `copy_dbs_to_local_data` - Flag (true/false) that indicates if large database files (gnomAD, etc) and apptainer container are copied to the `local_data` folder.
* `delete_fastq_files` - Flag (true/false) that indicates if FASTQ files are to be deleted after the BAM/CRAM file is created after mapping.
* `cnv_bin_size_wgs` - Bin size used for CNV analysis of WGS samples.
* `cnv_bin_size_shallow_wgs` - Bin size used for CNV analysis of shallow WGS samples.
* `cnv_bin_size_longread_wgs` - Bin size used for CNV analysis of long-read WGS samples.
* `use_freebayes` - Flag (true/false) that indicates if freebayes should be used instead of DeepVariant for short-read variant calling.
* `use_deepsomatic` Flag (true/false) that indicates if deepsomatic should be used for tumor-normal and tumor-only variant calling.
* `use_bwa1` - Flag (true/false) that indicates if BWA mem should be used for mapping instead of BWA mem 2 (BWA mem 2 is faster, but needs more RAM).
* `bwa_mem2_suffix` - Suffix appended to the BWA-mem2 executable, e.g. `.avx2` for `AMD EPYC9654`.
* `annotate_refseq_consequences` - Flag (true/false) that indicates if variant consequences based on RefSeq transcripts should be annotated in addition to variant consequences based on Ensembl transcripts.
* `custom_columns` - Used to add custom annotations to the output VCF/GSvar file. Each entry consists of a colon-speparated list of VCF file, INFO field name in the source VCF (prefixed with `CUSTOM_` in the annotated VCF) and column description. The VCF file has to be gzipped and indexed. Provide the annotation using the dictionary syntax of PHP, e.g. `custom_columns['NGSD_counts'] = "/mnt/data/dbs/NGSD/NGSD_germline.vcf.gz;COUNTS;NGSD counts"`.
* `locaton` - If set enables site-specific functionality and tests. Set only if you are a collaborator with a specific site name.
* `rna_allowed_systems` - Used to allow multiple processing systems in the RNA analysis for the cohort. If set the samples of all given processing systems are used to build the cohort and a batch correction is run to correct for differences. E.G. `rna_allowed_systems['processing_system_short_name1']="processing_system_short_name2,processing_system_short_name3"`
## Section 'dragen'

This section contains settings to run the germline/somatic data analysis on a on-site Illumina Dragen server:

* `dragen_version` - Version of Dragen to use, e.g. `4.3.17`.
* `dragen_in` - Folder into which input data is copied for Dragen data analysis. Only used for somatic tumor-normal analysis right now.
* `dragen_out` - Folder into which Dragen output files are copied after the analysis. Only used for somatic tumor-normal analysis right now.
* `dragen_data` - Folder used as working directory on the Dragen server.
* `dragen_genome` - Path to the Dragen genome reference hash tables.
* `dragen_log` - Folder used to write SGE logs of Dragen analyses.
* `queues_dragen` - Comma-separated list of SGE queues for Dragen servers (one per server).

## Section 'mysql-databases'

This section contains NGSD data (host, database name, credentials).

## Section 'grid_engine'

This section contains settings for automated data analysis via SGE:

* `queues_default` - Default SGE queues.
* `queues_research` - Research SGE queues.
* `queues_high_priority` - High-priority queues.
* `queues_high_mem` - High-memory consumption queues (RNA).
