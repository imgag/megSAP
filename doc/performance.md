# megSAP benchmarks

## single-sample short-read pipline

All performance benchmarks are performed on the GIAB reference sample NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1).

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Tools/validate_NA12878.php).

### Whole genome sequencing

The WGS samples were processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq6000 using 151PE.  
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
    <td>megSAP default</td> <!--- dataset: NA12878_45 @ 40x --->
 	<td>99.63%</td>
 	<td>99.43%</td>
 	<td>99.93%</td>
 	<td>96.75%</td>
 	<td>99.38%</td>
 	<td>97.95%</td>
  </tr>
  <tr>
    <td>megSAP DRAGEN v4.0.3 (-use_dragen)</td> <!--- dataset: NA12878_45 @ 40x --->
 	<td>99.75%</td>
 	<td>99.75%</td>
 	<td>99.98%</td>
 	<td>99.69%</td>
 	<td>99.66%</td>
 	<td>99.94%</td>
  </tr>
</table>

### Whole genome sequencing - non-default parameters

The WGS samples were processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq6000 using 151PE.  
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
		<td>default parameters (indel realignment, AF=0.1, MQ=20, BQ=20)</td> <!--- dataset: NA12878_45 @ 40x --->
		<td>99.70%</td>
		<td>99.44%</td>
		<td>99.94%</td>
		<td>96.77%</td>
		<td>99.38%</td>
		<td>97.95%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_noabra</td>
		<td>99.70%</td>
		<td>99.43%</td>
		<td>99.94%</td>
		<td>96.86%</td>
		<td>99.37%</td>
		<td>97.98%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_af005</td>
		<td>99.70%</td>
		<td>99.44%</td>
		<td>99.95%</td>
		<td>96.83%</td>
		<td>99.37%</td>
		<td>97.97%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_af015</td>
		<td>99.70%</td>
		<td>99.43%</td>
		<td>99.94%</td>
		<td>96.18%</td>
		<td>99.41%</td>
		<td>97.96%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_af020</td>
		<td>99.69%</td>
		<td>99.49%</td>
		<td>99.95%</td>
		<td>94.20%</td>
		<td>99.49%</td>
		<td>98.01%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_mq01</td>
		<td>99.71%</td>
		<td>99.42%</td>
		<td>99.94%</td>
		<td>96.74%</td>
		<td>99.39%</td>
		<td>97.95%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_mq40</td>
		<td>99.57%</td>
		<td>99.57%</td>
		<td>99.95%</td>
		<td>96.73%</td>
		<td>99.42%</td>
		<td>97.96%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_mq50</td>
		<td>98.04%</td>
		<td>99.68%</td>
		<td>99.95%</td>
		<td>96.18%</td>
		<td>99.45%</td>
		<td>97.96%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_bq10</td>
		<td>99.70%</td>
		<td>99.48%</td>
		<td>99.95%</td>
		<td>97.80%</td>
		<td>99.40%</td>
		<td>98.46%</td>
	</tr>
	<tr>
		<td>NA12878_45_40x_bq30</td>
		<td>99.69%</td>
		<td>99.47%</td>
		<td>99.94%</td>
		<td>96.05%</td>
		<td>99.41%</td>
		<td>97.74%</td>
	</tr>
</table>


### Exome sequencing

The WES NA12878 sample was processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content) and sequenced on NovaSeq6000 using 109PE.  
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
    <td>megSAP default</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.23%</td>
 	<td>98.30%</td>
 	<td>99.74%</td>
 	<td>96.85%</td>
 	<td>94.11%</td>
 	<td>96.23%</td>
  </tr>
  <tr>
    <td>megSAP DRAGEN (v4.0.3) (-use_dragen)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.21%</td>
 	<td>98.82%</td>
 	<td>99.83%</td>
 	<td>98.81%</td>
 	<td>97.26%</td>
 	<td>99.45%</td>
  </tr>
</table>


### Exome sequencing - non-default parameters

The WES samples were processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content).  
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
    <td>default parameters (indel realignment, AF=0.1, MQ=20, BQ=20)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.23%</td>
 	<td>98.30%</td>
 	<td>99.74%</td>
 	<td>96.85%</td>
 	<td>94.11%</td>
 	<td>96.23%</td>
  </tr>
  <tr>
    <td>No indel-realignment (-no_abra)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.23%</td>
 	<td>98.31%</td>
 	<td>99.74%</td>
 	<td>96.81%</td>
 	<td>93.81%</td>
 	<td>96.04%</td>
  </tr>
  <tr>
    <td>5% AF cutoff for variant calling (-min_af 0.05)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.22%</td>
 	<td>98.33%</td>
 	<td>99.74%</td>
 	<td>97.13%</td>
 	<td>94.04%</td>
 	<td>96.38%</td>
  </tr>
  <tr>
    <td>15% AF cutoff for variant calling (-min_af 0.15)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.23%</td>
 	<td>98.31%</td>
 	<td>99.72%</td>
 	<td>96.12%</td>
 	<td>95.60%</td>
 	<td>96.20%</td>
  </tr>
  <tr>
    <td>20% AF cutoff for variant calling (-min_af 0.2)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.17%</td>
 	<td>98.68%</td>
 	<td>99.73%</td>
 	<td>94.62%</td>
 	<td>97.46%</td>
 	<td>96.24%</td>
  </tr>
  <tr>
    <td>Minimum MQ of 1 for variant calling (-min_mq 1)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.24%</td>
 	<td>98.23%</td>
 	<td>99.73%</td>
 	<td>96.81%</td>
 	<td>94.19%</td>
 	<td>96.23%</td>
  </tr
  <tr>
    <td>Minimum MQ of 40 for variant calling (-min_mq 40)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>98.85%</td>
 	<td>98.93%</td>
 	<td>99.83%</td>
 	<td>96.58%</td>
 	<td>94.38%</td>
 	<td>96.22%</td>
  </tr>
  <tr>
    <td>Minimum MQ of 50 for variant calling (-min_mq 50)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>96.45%</td>
 	<td>99.36%</td>
 	<td>99.85%</td>
 	<td>95.21%</td>
 	<td>94.65%</td>
 	<td>96.17%</td>
  </tr>
  <tr>
    <td>Minimum BQ of 10 for variant calling (-min_bq 10)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.21%</td>
 	<td>98.35%</td>
 	<td>99.74%</td>
 	<td>97.49%</td>
 	<td>93.77%</td>
 	<td>96.82%</td>
  </tr>
  <tr>
    <td>Minimum BQ of 30 for variant calling (-min_bq 30)</td> <!--- dataset: NA12878x2_82 @ 113x --->
 	<td>99.21%</td>
 	<td>98.38%</td>
 	<td>99.73%</td>
 	<td>94.94%</td>
 	<td>95.37%</td>
 	<td>95.15%</td>
  </tr>
</table>

## single-sample long-read pipline
The analyses were performed with the [long-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze_longread.php) on the GRCh38 reference genome.

The lrGS samples were processed with the Oxford Nanopore Tech. Ligation Sequencing Kit V14 (SQK-LSK114).  
The benchmarks were performed on GIAB high-confidence regions **with at least 3x coverage**.
 <!--- dataset: 23014LRa023L2_01, 60x average coverge --->

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
    <td>default parameters</td>
    <td>99.92%</td>
    <td>99.78%</td>
    <td>99.94%</td>
    <td>78.38%</td>
    <td>95.67%</td>
    <td>99.11%</td>
  </tr>
  <tr>
    <td>validation only on coding region</td>
    <td>99.96%</td>
    <td>99.81%</td>
    <td>99.97%</td>
    <td>93.22%</td>
    <td>97.56%</td>
    <td>99.77%</td>
  </tr>
</table>
