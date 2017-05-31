#megSAP - DNA analysis (tumor-normal pair)

###Basics

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
* `t_sys` - The [processing system INI file](processing_system_ini_file.md) for the tumor sample.
* `n_sys` - The [processing system INI file](processing_system_ini_file.md) for the normal sample.
* `nsc` - Skip sample correlation check. This is useful if sample correlation is low and the pipeline will give an error otherwise.

###Output

coming soon

[back to the start page](../README.md)


