# megSAP - Processing system INI file

The processing system INI file described the wet-lab processing of the sample (adapter sequences, target region, data type, etc) and defines the main parameters for the data analysis:

* `name_short` - Processing system short name (must be a valid file name).
* `name_manufacturer` - Processing system full name (can be any name, including characters that are invalid in file names).
* `target_file` - Target region BED file path. The target region is used to determine where indel realignment and variant calling are done.
* `adapter1_p5` - Read 1 adapter sequence (Illumina standard is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`).
* `adapter2_p7` - Read 2 adapter sequence (Illumina standard is `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`).
* `type` - Processing system type: `WGS`, `WES`, `Panel`, `Panel Haloplex`, `Panel MIPs` or `RNA`.
* `shotgun` - `true` for randomly-fragmented reads,  `false` for amplicon-based reads.
* `umi_type` - Unique molecular identifier type: `n/a`, `HaloPlex HS`, `SureSelect HS`, `ThruPLEX`, `Safe-SeqS` or `MIPs`.
* `build` - Only `GRCh38` is supported.

Notes for the [RNA analysis pipeline](rna_single_sample.md):

* `target_file` - Target region is only used for mapping quality control.
* `build` - Basename of the genome build. FASTA, STAR index and GTF annotation have to be present.

## Examples

Example processing system INI files are provided here:

* WES: [Agilent SureSelectXT Human All Exon V7](ini_examples/ssHAEv7.ini)
* WGS: [Illumina TruSeq DNA PCR-Free](ini_examples/TruSeqPCRfree.ini) (for WGS you can use this [generic WGS target region](ini_examples/WGS_hg38.bed) independent of the kit)
* lrGS: [ONT Amplicon by Ligation ONT SQK-LSK114](ini_examples/PCR-SQK-LSK114.ini)

**Note:** For WES and WGS you need to adapt the paths to the target region before using these examples!
	
	
	

[back to the start page](../README.md)
