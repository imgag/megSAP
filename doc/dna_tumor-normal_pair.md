#megSAP - DNA analysis (tumor-normal pair)

##Basics

Tumor-normal pairs can be analyzed by using the somatic_dna.php script. This script is also used by several project specific data analysis pipelines. For example, the script `somatic_emed.php` can be used for analysis of combined DNA/RNA tumor-normal data. 

Please also look at the help using:

	> php megSAP/src/Pipelines/somatic_dna.php --help

The main parameters that you have to provide are:

* `p_folder` - Project folder that contains all sample folders that will be used for the current analysis. The sample folders should be named like `Sample_ID` while ID is the corresponding sample ID. Sample IDs are given by the parameters `t_id` and `n_id`. Sample folders contain the FASTQ files produced by bcl2fastq2.
* `t_id` - The tumor sample name, which must be a prefix of the FASTQ files.
* `n_id` -  The normal sample name, which must be a prefix of the FASTQ files.
* `o_folder` - The output folder.
* `steps` - Analysis steps to perform. Please use `ma,vc,an` to perform mapping, variant calling and variant annotation.
* `filter_set` - Filter set to use for post-call variant filtering. Applies only if annotation step is selected. Multiple filters can be comma separated. Same options like in `filter_vcf.php`
* `min_af` - Minimum variant allele frequency to detect.
* `t_sys` - The processing system configuration INI file for the tumor sample (see next section).
* `n_sys` - The processing system configuration INI file for the normal sample (see next section).
* `nsc` - Skip sample correlation check. This is useful of correlation is low and the pipeline will give an error otherwise.

##Processing system configuration INI file

The system INI file described the wet-lab processing of the sample and defines the main parameters for the data analysis:

* `name_short` - Processing system short name (must be a valid file name).
* `name_manufacturer` - Processing system full name.
* `target_file` - Target region BED file path.
* `adapter1_p5` - Read 1 adapter sequence (Illumina standard is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`).
* `adapter2_p7` - Read 1 adapter sequence (Illumina standard is `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`).
* `type` - Processing system type: `WGS`, `WES`, `Panel` or `Panel Haloplex`.
* `shotgun` - `true` for randomly-fragmented reads,  `false` for amplicon-based reads.
* `build` - Currently only 'GRCh37' is supported.

[back to the start page](../README.md)
