# megSAP benchmarks

## single-sample short-read pipline

All performance benchmarks are performed on the GIAB reference sample NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1).

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Tools/validate_NA12878.php).

### Whole genome sequencing

The WGS samples were processed with the Illumina TruSeq DNA PCR-Free kit.  
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
    <td>NovaSeq6000, 151PE, megSAP default</td> <!--- dataset: NA12878_45 @ 40x --->
    <td>99.62%</td>
    <td>99.42%</td>
    <td>99.93%</td>
    <td>96.74%</td>
    <td>99.38%</td>
    <td>97.94%</td>
  </tr>
  <tr>
    <td>NovaSeq6000, 151PE, megSAP DRAGEN v4.0.3 (-use_dragen)</td> <!--- dataset: NA12878_45 @ 40x --->
	<td>99.74%</td>
    <td>99.75%</td>
    <td>99.97%</td>
    <td>99.68%</td>
    <td>99.66%</td>
    <td>98.93%</td>
  </tr>
</table>


### Exome sequencing

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
    <td>NovaSeq6000, 109PE, megSAP default</td> <!--- dataset: NA12878x2_82 @ 113x --->
    <td>99.22%</td>
    <td>98.30%</td>
    <td>99.73%</td>
    <td>96.82%</td>
    <td>94.10%</td>
    <td>96.23%</td>
  </tr>
  <tr>
    <td>NovaSeq6000, 109PE, megSAP DRAGEN (v4.0.3) (-use_dragen)</td> <!--- dataset: NA12878x2_82 @ 113x --->
    <td>99.21%</td>
    <td>98.81%</td>
    <td>99.82%</td>
    <td>98.81%</td>
    <td>97.26%</td>
    <td>99.44%</td>
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
    <td>default parameters</td>
    <td>99.22%</td>
    <td>98.30%</td>
    <td>99.73%</td>
    <td>96.82%</td>
    <td>94.10%</td>
    <td>96.23%</td>
  </tr>
  <tr>
    <td>No indel-realignment (-no_abra)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>5% AF cutoff for variant calling (-min_af 0.05)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>15% AF cutoff for variant calling (-min_af 0.15)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>20% AF cutoff for variant calling (-min_af 0.2)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>Minimum mapping quality of 1 for variant calling (-min_mq 1)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr
  <tr>
    <td>Minimum mapping quality of 40 for variant calling (-min_mq 40)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>Minimum mapping quality of 50 for variant calling (-min_mq 50)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>Minimum base quality of 1 for variant calling (-min_bq 1)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
  </tr>
  <tr>
    <td>Minimum base quality of 30 for variant calling (-min_bq 30)</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
    <td>TODO</td>
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
