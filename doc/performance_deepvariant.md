# megSAP benchmarks DeepVariant

## single-sample short-read pipline

All performance benchmarks are performed on the GIAB reference sample NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1).

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Auxilary/validate_NA12878.php).

### Whole genome sequencing

The WGS sample was processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq6000 using 151PE at 40x average depth.  
The benchmarks were performed on GIAB high-confidence regions **with at least 15x coverage**.

<table>
	<tr>
		<th rowspan=2>Test</th>
		<th colspan=3>SNVs</th>
		<th colspan=3>InDels</th>
	</tr>
	<tr>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>genotyping</th>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>genotyping</th>
	</tr>
	<tr>
		<td>megSAP freebayes v1.3.6</td> <!--- dataset: NA12878x2_93 @ 40x --->
		<td>99.62%</td>
		<td>99.45%</td>
		<td>99.94%</td>
		<td>97.48%</td>
		<td>99.42%</td>
		<td>98.46%</td>
	</tr>
	<tr>
		<td>megSAP DeepVariant v1.8.0</td> <!--- dataset: NA12878x2_93 @ 40x --->
		<td>99.65%</td>
		<td>99.93%</td>
		<td>99.98%</td>
		<td>98.77%</td>
		<td>99.72%</td>
		<td>99.81%</td>
	</tr>
	<tr>
		<td>megSAP DRAGEN v4.2.4</td> <!--- dataset: NA12878x2_93 @ 40x --->
		<td>99.74%</td>
		<td>99.74%</td>
		<td>99.97%</td>
		<td>99.68%</td>
		<td>99.66%</td>
		<td>99.91%</td>
	</tr>
</table>

### Exome sequencing (Twist)

The WES sample was processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content) and sequenced on NovaSeqXPlus.   
The benchmarks were performed on GIAB high-confidence regions **with at least 15x coverage**.

<table>
	<tr>
		<th rowspan=2>Test</th>
		<th colspan=3>SNV</th>
		<th colspan=3>InDel</th>
	</tr>
	<tr>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>genotyping</th>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>genotyping</th>
	</tr>
	<tr>
		<td>megSAP freebayes v1.3.6</td> <!--- dataset: NA12878x3_23 @ 188x --->
		<td>99.00%</td>
		<td>98.48%</td>
		<td>99.78%</td>
		<td>96.79%</td>
		<td>96.57%</td>
		<td>96.07%</td>
	</tr>
	<tr>
		<td>megSAP DeepVariant v1.8.0</td> <!--- dataset: NA12878x3_23 @ 188x --->
		<td>99.02%</td>
		<td>99.93%</td>
		<td>99.91%</td>
		<td>97.96%</td>
		<td>99.45%</td>
		<td>99.86%</td>
	<tr>
	<tr>
		<td>megSAP DRAGEN v4.2.4</td> <!--- dataset: NA12878x3_23 @ 188x --->
		<td>99.31%</td>
		<td>98.57%</td>
		<td>99.78%</td>
		<td>99.09%</td>
		<td>97.20%</td>
		<td>99.45%</td>
	</tr>
</table>

### Runtime and memory usage


| sample     | type | caller  | threads | runtime (hh:mm) | memory usage |
|------------|------|------------|--------:|----------------:|-------------:|
| NA12878x2_93 | WGS  |freebayes v1.3.6 	|     1   | 08:86            | 3.31 GB      |
| NA12878x2_93 | WGS  |freebayes v1.3.6		|     5   | 01:50            | 4.78 GB       |
| NA12878x2_93 | WGS  |freebayes v1.3.6  	|     10  | 01:35            | 5.81 GB      |
| NA12878x2_93 | WGS  |DeepVariant v1.8.0	|     1   | 70:47            | 14.96 GB      |
| NA12878x2_93 | WGS  |DeepVariant v1.8.0  	|     5   | 15:55            | 16.28 GB      |
| NA12878x2_93 | WGS  |DeepVariant v1.8.0	|     10  | 08:36       	 	| 17.51 GB		|
| NA12878x3_23 | WES  |freebayes v1.3.6 	|     1   | 01:17            | 860 MB      |
| NA12878x3_23 | WES  |freebayes v1.3.6		|     5   | 00:17            | 2.2 GB      |
| NA12878x3_23 | WES  |freebayes v1.3.6  	|     10  | 00:10            | 3.33 GB      |
| NA12878x3_23 | WES  |DeepVariant v1.8.0	|     1   | 03:53            | 9.11 GB      |
| NA12878x3_23 | WES  |DeepVariant v1.8.0  	|     5   | 01:11            | 9.32 GB      |
| NA12878x3_23 | WES  |DeepVariant v1.8.0	|     10  | 00:41       		| 9.46 GB		|
