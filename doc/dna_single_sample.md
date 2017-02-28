#megSAP - DNA analysis (single sample)

##Basics

Single sample DNA analysis is performed using the `analyze.php` script.  
Please have a look at the help using:

	> php megSAP/src/Pipelines/analyze.php --help

The main parameters that you have to provide are:

* `folder` - The sample folder, which contains the the FASTQ files as produced by bcl2fastq2.
* `name` - The sample name, which must be a prefix of the FASTQ files.
* `steps` -  Analysis steps to perform. Please use `ma,vc,an` to perform mapping, variant calling and variant annotation.
* `system` - The processing system configuration INI file (see next section).


##Processing system configuration INI file

The system INI file described the wet-lab processing of the sample and defines the main parameters for the data analysis:

* `name_short` - Processing system short name (must be a valid file name).
* `name_manufacturer` - Processing system full name.
* `target_file` - Target region BED file path.
* `adapter1_p5` - Read 1 adapter sequence (Illumina standard is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`).
* `adapter2_p7` - Read 1 adapter sequence (Illumina standard is `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`).
* `type` - Processing system type: `WGS`, `WES`, `Panel` or `Panel Haloplex`.
* `shotgun` - `true` for randomly-fragmented reads,  `false` for amplicon-based reads.
* `build` - Currently only 'hg19' is supported.

[back to the start page](../Readme.md)

