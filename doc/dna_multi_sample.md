# megSAP - DNA analysis (multi-sample and trio)

## Multi-sample pipeline

### Basics

Multi-sample DNA analysis is performed using the `multisample.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/multisample.php --help

The main parameters that you have to provide are:

* `bams` - A list of input BAM files.
* `status` - A comma-separated list of affected status corresponding to the input BAM files.
* `out_folder` -  Output folder.
* `system` - The processing system configuration INI file (see next section).

*Note:* The processing system of the *first* sample is used to determine the target region for the multi-sample and trio analyses.

### Tools used in this analysis pipeline

Mapping is not part of this pipline - it has to be perfomed beforehand using the single-sample pipeline.
The tools used for variant calling and annotation are the same as for the [single-sample pipeline](dna_single_sample.md).

### Output

After the analysis, these files are created in the output folder:

1. a multi-sample variant list `all.vcf.gz` in VCF format.
2. a multi-sample variant list `multi.GSvar` in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md).

## Trio pipeline

### Basics

Trio DNA analysis is a special case of the multi-sample analysis for an **affected child and healthy** parents.  
It is performed using the `trio.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/trio.php --help

The main parameters that you have to provide are:

* `f` - BAM file of father.
* `m` - BAM file of mother.
* `c` - BAM file of child (index).
* `out_folder` -  Output folder.
* `system` - The processing system configuration INI file (see next section).

### Tools used in this analysis pipeline

Mapping is not part of this pipline - it has to be perfomed beforehand using the single-sample pipeline.
The tools used for variant calling and annotation are the same as for the [single-sample pipeline](dna_single_sample.md).

### Output

After the analysis, these files are created in the output folder:

1. a multi-sample variant list `all.vcf.gz` in VCF format.
2. a multi-sample variant list `trio.GSvar` in [GSvar format](https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/gsvar_format.md).

[back to the start page](../README.md)
