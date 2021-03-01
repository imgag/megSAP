# megSAP benchmarks

## megSAP single-sample pipline

All performance benchmarks are performed on the GIAB reference sample NA12878.  
The analyses were performed with the [single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh37 reference genome.

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Tools/validate_NA12878.php).

### Whole genome sequencing

The WGS samples were processed with the Illumina TruSeq DNA PCR-Free kit.  
The benchmarks were performed on GIAB high-confidence regions - **no depth cutoff was used**.
 <!--- dataset: NA12878_45, 95x average coverge --->

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
    <td>99.96%</td>
    <td>99.72%</td>
    <td>99.98%</td>
    <td>97.73%</td>
    <td>99.52%</td>
    <td>98.34%</td>
  </tr>
  <tr>
    <td>mapping with DRAGEN (-use_dragen)</td>
    <td>99.96%</td>
    <td>99.74%</td>
    <td>99.98%</td>
    <td>97.71%</td>
    <td>99.52%</td>
    <td>98.36%</td>
  </tr>
</table>


### Exome sequencing

The WES samples were processed with the Agilent SureSelectXT Human All Exon V7 kit.  
The benchmarks were performed on GIAB high-confidence regions - **no depth cutoff was used**.
 <!--- dataset: NA12878_58, 150x average coverge --->

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
    <td>99.62%</td>
    <td>99.53%</td>
    <td>99.97%</td>
    <td>93.13%</td>
    <td>95.05%</td>
    <td>96.50%</td>
  </tr>
  <tr>
    <td>mapping with DRAGEN (-use_dragen)</td>
    <td>99.61%</td>
    <td>99.59%</td>
    <td>99.97%</td>
    <td>93.23%</td>
    <td>94.96%</td>
    <td>96.45%</td>
  </tr>
</table>


### Exome sequencing - non-default parameters

The WES samples were processed with the Agilent SureSelectXT Human All Exon V7 kit.  
The benchmarks were performed on GIAB high-confidence regions **with at least 20x coverage**.
 <!--- dataset: NA12878_58, 150x average coverge --->

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
    <td>99.96%</td>
    <td>99.60%</td>
    <td>99.98%</td>
    <td>96.56%</td>
    <td>94.10%</td>
    <td>96.22%</td>
  </tr>
  <tr>
    <td>In indel-realignment (-no_abra)</td>
    <td>99.94%</td>
    <td>99.57%</td>
    <td>99.98%</td>
    <td>94.37%</td>
    <td>95.01%</td>
    <td>94.14%</td>
  </tr>
  <tr>
    <td>No adapter/quality trimming (-no_trim)</td>
    <td>99.95%</td>
    <td>99.63%</td>
    <td>99.98%</td>
    <td>96.56%</td>
    <td>94.15%</td>
    <td>96.22%</td>
  </tr>
  <tr>
    <td>5% AF cutoff for variant calling (-min_af 0.05)</td>
    <td>99.95%</td>
    <td>99.65%</td>
    <td>99.98%</td>
    <td>96.77%</td>
    <td>93.64%</td>
    <td>96.07%</td>
  </tr>
  <tr>
    <td>15% AF cutoff for variant calling (-min_af 0.15)</td>
    <td>99.95%</td>
    <td>99.56%</td>
    <td>99.97%</td>
    <td>95.46%</td>
    <td>95.51%</td>
    <td>96.28%</td>
  </tr>
  <tr>
    <td>20% AF cutoff for variant calling (-min_af 0.2)</td>
    <td>99.94%</td>
    <td>99.58%</td>
    <td>99.97%</td>
    <td>93.95%</td>
    <td>97.56%</td>
    <td>96.67%</td>
  </tr>
</table>
