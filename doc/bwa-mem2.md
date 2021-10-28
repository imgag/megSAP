# bwa-mem2

To map FastQ files megSAP can also use `bwa-mem2`. It is faster than the original `bwa-mem` implementation and can be found here: 
[https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)


## installation and usage

`bwa-mem2` is automatically installed besides the standard `bwa` during the [initial setup](install_unix.md#initial-setup) (it is part of the `download_tools.sh` script in `data/`) and is activated by default. To use the standard `bwa` you simply have to change the `use_bwa1` parameter in the megSAP `settings.ini` to `true`.

## update from older megSAP version
`bwa-mem2` requires additional index files for reference genome. If you update from an older megSAP version these index files have to be created either by (re-)run the reference download script:
	
	> data/download_GRCh38.sh

Or run the genome index script:

    > php src/Tools/index_genome.php -in reference.fa


**Note:** Since the index generation of `bwa-mem2` is automatically adapted to the available CPU features the index should be created by the same machine which is used for mapping.



## performance comparison

`bwa-mem2` is faster than the original `bwa` algorithm, but requires more memory. A performance comparison of `bwa` and `bwa-mem2` in the `megSAP` pipeline (Standard analysis with mapping step, but without InDel realignment) is shown below.     
The NA12878 WGS sample was sequenced with the illumina TruSeq DNA PCR-Free kit and the NA12878 WES sample was sequences using the Agilent SureSelectXT Human All Exon V7 kit.

| sample     | type | algorithm  | threads | runtime (hh:mm) | memory usage |
|------------|------|------------|--------:|----------------:|-------------:|
| NA12878_58 | WES  |bwa 0.7.14  |     5   |             1:34|       6.28 GB|
| NA12878_58 | WES  |bwa-mem2 2.1|     5   |             1:02|      18.04 GB|
| NA12878_45 | WGS  |bwa 0.7.14  |     5   |            62:25|      20.04 GB|
| NA12878_45 | WGS  |bwa-mem2 2.1|     5   |            37:06|      22.04 GB|

