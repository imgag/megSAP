# megSAP - DNA analysis (multi-sample and trio)

## multi-sample pipeline

### Basics

Multi-sample DNA analysis is performed using the `multisample.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/multisample.php --help

The main parameters that you have to provide are:

* `bams` - A list of input BAM files.
* `status` - A comma-separated list of affected status corresponding to the input BAM files.
* `out_folder` -  Output folder.
* `system` - The processing system configuration INI file (see next section).

### Output

After the analysis, these files are created in the output folder:

1. a multi-sample variant list `all.vcf.gz` in VCF format.
2. a multi-sample variant list `multi.GSvar` in [GSvar format](gsvar_format.md).

## trio pipeline

### Basics

Trio DNA analysis is performed using the `trio.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/trio.php --help

The main parameters that you have to provide are:

* `f` - BAM file of father.
* `m` - BAM file of mother.
* `c` - BAM file of child (index).
* `out_folder` -  Output folder.
* `system` - The processing system configuration INI file (see next section).

### Output

After the analysis, these files are created in the output folder:

1. a multi-sample variant list `all.vcf.gz` in VCF format.
2. a multi-sample variant list `trio.GSvar` in [GSvar format](gsvar_format.md).

[back to the start page](../README.md)





