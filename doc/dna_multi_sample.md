# megSAP - DNA analysis (multi-sample and trio)

## Prerequisite: Single sample analysis

The multi-sample and trio pipelines require files form the sinlge sample analysis (BAM, VCF, ...).  
Thus, each samples has to be analyzed using the [single sample analysis](dna_single_sample.md) first.  
After the single-sample analyses, the multi-sample or trio analysis is perfomed.

## Multi-sample pipeline

### Basics

Multi-sample DNA analysis is performed using the `multisample.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/multisample.php --help

The main parameters that you have to provide are:

* `bams` - A list of input BAM files.
* `status` - A comma-separated list of affected status corresponding to the input BAM files.
* `out_folder` -  Output folder.
* `system` - The [processing system INI file](processing_system_ini_file.md).

*Note:* The processing system of the *first* sample is used to determine the target region for the multi-sample and trio analyses.

### Output

After the analysis, these files are created in the output folder:

1. small variants VCF/GSvar file
1. CNV TSV file
1. structural variants VCF/BEDPE file

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
* `system` - The [processing system INI file](processing_system_ini_file.md).

### Output

The same files as in the multi-sample analysis are produced.

[back to the start page](../README.md)

