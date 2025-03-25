# megSAP benchmarks

## Small variant calling benchmarks

All small variant benchmarks are done on the GIAB reference sample NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the short-read and long-read single sample pipelines on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1).

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Auxilary/validate_NA12878.php).

### short-read WES
<!--- NA12878x2_80 --->

The sample was processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content) and sequenced on NovaSeq6000 using 109PE at 113x average depth.  
All benchmarks were performed on GIAB high-confidence regions **with at least 15x coverage**.

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
		<td>megSAP 2025_03 - bwa-mem2, ABRA2, freebayes</td>
		<td>99.21%</td>
		<td>98.33%</td>
		<td>99.74%</td>
		<td>97.49%</td>
		<td>93.98%</td>
		<td>96.49%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - bwa-mem2, DeepVariant</td>
		<td>99.37%</td>
		<td>99.86%</td>
		<td>99.95%</td>
		<td>97.85%</td>
		<td>99.63%</td>
		<td>99.81%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - DRAGEN 4.2.4 no ML model</td>
		<td>99.22%</td>
		<td>98.82%</td>
		<td>99.83%</td>
		<td>98.81%</td>
		<td>97.18%</td>
		<td>99.45%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - DRAGEN 4.2.4 with ML model</td>
		<td>98.63%</td>
		<td>99.85%</td>
		<td>99.90%</td>
		<td>98.18%</td>
		<td>98.90%</td>
		<td>99.54%</td>
	</tr>
</table>

### short-read WGS
<!--- NA12878x3_28--->

The sample was processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq X Plus using 159PE at 35x average depth.  
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
		<td>megSAP 2025_03 - bwa-mem2, ABRA2, freebayes</td>
		<td>99.66%</td>
		<td>99.29%</td>
		<td>99.94%</td>
		<td>98.05%</td>
		<td>99.30%</td>
		<td>98.29%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - bwa-mem2, DeepVariant</td>
		<td>99.73%</td>
		<td>99.92%</td>
		<td>99.98%</td>
		<td>99.23%</td>
		<td>99.70%</td>
		<td>99.84%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - DRAGEN 4.2.4 no ML model</td>
		<td>99.78%</td>
		<td>99.70%</td>
		<td>99.98%</td>
		<td>99.64%</td>
		<td>99.60%</td>
		<td>99.85%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - DRAGEN 4.2.4 with ML model</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
	</tr>
</table>


### long-read WGS
<!--- 24067LRa008_01 --->

The sample was processed with the Oxford Nanopore Tech. Ligation Sequencing Kit V14 (SQK-LSK114) and sequenced at 42x average depth.  
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
		<td>megSAP 2025_03 (high accuracy)</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 (high accuracy) - coding region only</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 (super accuracy)</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
	</tr>
</table>

### small variants CMRG benchmarks
<!--- NA24385_32 / 24067LRa002_01 --->

For genome sequening, we also performed the [CMRG benchmark](https://www.nature.com/articles/s41587-021-01158-1).  

The short-read WGS sample was processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq X Plus using 159PE at 39x average depth.  
The long-read WGS sample  with the Oxford Nanopore Tech. Ligation Sequencing Kit V14 (SQK-LSK114) and sequenced at 40x average depth.  

All benchmarks were performed on GIAB high-confidence regions **with at least 15x coverage**.

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
		<td>megSAP 2025_03 - bwa-mem2, ABRA2, freebayes</td>
		<td>98.66%</td>
		<td>96.55%</td>
		<td>99.36%</td>
		<td>90.42%</td>
		<td>92.77%</td>
		<td>95.78%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - bwa-mem2, DeepVariant</td>
		<td>98.13%</td>
		<td>99.45%</td>
		<td>99.71%</td>
		<td>93.22%</td>
		<td>93.80%</td>
		<td>99.11%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - DRAGEN 4.2.4 no ML model</td>
		<td>98.10%</td>
		<td>98.48%</td>
		<td>99.73%</td>
		<td>94.59%</td>
		<td>92.06%</td>
		<td>99.27%</td>
	</tr>
	<tr>
		<td>megSAP 2025_03 - long-read WGS</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
		<td>TODO%</td>
	</tr>
</table>

## Structural variant calling benchmarks

source sample: `HG002` (pre-release truth set)  
caller: `Manta` (short-read), `Sniffles 2` (long-read)  
validation tool: `hap_eval`

<table>
	<tr>
		<th>Sample</th>
		<th>processing system</th>
		<th>coverage</th>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>F1</th>
	</tr>
	<tr>
		<td>short-read WGS</td>
		<td>TruSeq DNA PCR-Free</td>
		<td>43x</td>
		<td>32.84%</td>
		<td>96.96%</td>
		<td>47.28%</td>
	</tr>
	<tr>
		<td>long-read WGS (R10)</td>
		<td>ONT Ligation Sequencing Kit V14</td>
		<td>51x</td>
		<td>91.81%</td>
		<td>96.27%</td>
		<td>92.03%</td>
	</tr>
</table>
