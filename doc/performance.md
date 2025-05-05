# Benchmarks

All benchmarks are perfomed with the megSAP release [2025_03](https://github.com/imgag/megSAP/releases/tag/2025_03).  
As reference genome GRCh38 with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1) was used.

## Small variant calling benchmarks

All small variant benchmarks are done on the GIAB reference sample NA12878/HG001 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the short-read and long-read single sample pipelines.

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
		<td>short-read WES - bwa-mem2, ABRA2, freebayes</td>
		<td>99.21%</td>
		<td>98.33%</td>
		<td>99.74%</td>
		<td>97.49%</td>
		<td>93.98%</td>
		<td>96.49%</td>
	</tr>
	<tr>
		<td>short-read WES - bwa-mem2, DeepVariant</td>
		<td>99.37%</td>
		<td>99.86%</td>
		<td>99.95%</td>
		<td>97.85%</td>
		<td>99.63%</td>
		<td>99.81%</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.3.13 no ML model</td>
		<td>99.23%</td>
		<td>98.77%</td>
		<td>99.83%</td>
		<td>98.77%</td>
		<td>97.47%</td>
		<td>99.49%</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.3.13 with ML model</td>
		<td>98.87%</td>
		<td>99.72%</td>
		<td>99.95%</td>
		<td>98.49%</td>
		<td>98.58%</td>
		<td>99.44%</td>
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
		<td>short-read WGS - bwa-mem2, ABRA2, freebayes</td>
		<td>99.66%</td>
		<td>99.29%</td>
		<td>99.94%</td>
		<td>98.05%</td>
		<td>99.30%</td>
		<td>98.29%</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>99.73%</td>
		<td>99.92%</td>
		<td>99.98%</td>
		<td>99.23%</td>
		<td>99.70%</td>
		<td>99.84%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13 no ML model</td>
		<td>99.78%</td>
		<td>99.70%</td>
		<td>99.98%</td>
		<td>99.65%</td>
		<td>99.61%</td>
		<td>99.85%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13 with ML model</td>
		<td>99.80%</td>
		<td>99.81%</td>
		<td>99.99%</td>
		<td>99.68%</td>
		<td>99.55%</td>
		<td>99.89%</td>
	</tr>
</table>

#### Note on DeepVariant GPU version
<!--- NA12878x2_93 --->
Calling variants with DeepVariant consists of three steps: `make_examples`, `call_variants`, `postprocess`.  
GPU acceleration only affects the `call_variants` step. Since `make_examples` takes up most of the overall runtime we don't provide a DeepVariant GPU option for now.

Step-wise runtime comparison: GPU vs. CPU (WGS):

| Step           | GPU Time | CPU Time |
|----------------|------------|------------|
| `make_examples` | 377min19s | 466min 18s |
| `call_variants` | 10min 44s   | 43min 33s    |
| `postprocess`   | 2min 48s | 3min 29s        |

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
		<td>long-read WGS (high accuracy)</td>
		<td>99.84%</td>
		<td>99.83%</td>
		<td>99.98%</td>
		<td>78.61%</td>
		<td>97.40%</td>
		<td>99.30%</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>99.82%</td>
		<td>99.87%</td>
		<td>99.99%</td>
		<td>89.58%</td>
		<td>96.19%</td>
		<td>98.58%</td>
	</tr>
</table>

	
### coding region benchmark

All small variant benchmarks above are done on regions with at least with 15x coverage.  
To allow a comparison of WES, WGS and lrGS independent of the coverage, we also perfomed a benchmark without depth cutoff on the coding region of all protein-coding genes padded by two bases to include the consensus splice sites.

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
		<td>short-read WES - bwa-mem2, ABRA2, freebayes</td>
		<td>97.77%</td>
		<td>98.39%</td>
		<td>99.75%</td>
		<td>89.15%</td>
		<td>95.44%</td>
		<td>97.37%</td>
	</tr>
	<tr>
		<td>short-read WES - bwa-mem2, DeepVariant</td>
		<td>97.98%</td>
		<td>99.95%</td>
		<td>99.98%</td>
		<td>91.28%</td>
		<td>99.54%</td>
		<td>99.53%</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.3.13 no ML model</td>
		<td>97.74%</td>
		<td>98.80%</td>
		<td>99.85%</td>
		<td>92.13%</td>
		<td>95.16%</td>
		<td>99.54%</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.3.13 with ML model</td>
		<td>97.57%</td>
		<td>99.73%</td>
		<td>99.94%</td>
		<td>91.28%</td>
		<td>97.72%</td>
		<td>99.77%</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, ABRA2, freebayes</td>
		<td>99.20%</td>
		<td>98.68%</td>
		<td>99.85%</td>
		<td>95.96%</td>
		<td>97.62%</td>
		<td>98.67%</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>99.20%</td>
		<td>99.91%</td>
		<td>99.96%</td>
		<td>97.02%</td>
		<td>99.56%</td>
		<td>99.78%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13 no ML model</td>
		<td>99.26%</td>
		<td>99.11%</td>
		<td>99.91%</td>
		<td>98.51%</td>
		<td>97.47%</td>
		<td>99.78%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13 with ML model</td>
		<td>99.31%</td>
		<td>99.68%</td>
		<td>99.95%</td>
		<td>98.30%</td>
		<td>97.88%</td>
		<td>99.78%</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>99.91%</td>
		<td>99.80%</td>
		<td>99.99%</td>
		<td>92.77%</td>
		<td>98.20%</td>
		<td>99.54%</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>99.83%</td>
		<td>99.75%</td>
		<td>99.99%</td>
		<td>96.60%</td>
		<td>98.27%</td>
		<td>99.78%</td>
	</tr>
</table>

### CMRG benchmark
<!--- NA24385_32 / 24067LRa002_01 --->

For genome sequening, we also performed the [CMRG benchmark](https://www.nature.com/articles/s41587-021-01158-1) based on the NA24385/HG002 sample.

The short-read WGS sample was processed with the Illumina TruSeq DNA PCR-Free kit and sequenced on NovaSeq X Plus using 159PE at 39.5x average depth.  
The long-read WGS sample with the Oxford Nanopore Tech. Ligation Sequencing Kit V14 (SQK-LSK114) and sequenced at 40.5x average depth.  

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
		<td>short-read WGS - bwa-mem2, ABRA2, freebayes</td>
		<td>98.66%</td>
		<td>96.55%</td>
		<td>99.36%</td>
		<td>90.42%</td>
		<td>92.77%</td>
		<td>95.78%</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>98.13%</td>
		<td>99.45%</td>
		<td>99.71%</td>
		<td>93.22%</td>
		<td>93.80%</td>
		<td>99.11%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13 no ML model</td>
		<td>98.12%</td>
		<td>98.50%</td>
		<td>99.69%</td>
		<td>94.58%</td>
		<td>92.00%</td>
		<td>99.27%</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>98.83%</td>
		<td>95.78%</td>
		<td>99.77%</td>
		<td>65.92%</td>
		<td>86.58%</td>
		<td>99.21%</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>98.83%</td>
		<td>95.78%</td>
		<td>99.77%</td>
		<td>65.92%</td>
		<td>86.58%</td>
		<td>99.21%</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>98.76%</td>
		<td>96.39%</td>
		<td>99.87%</td>
		<td>78.20%</td>
		<td>85.24%</td>
		<td>97.56%</td>
	</tr>
</table>

## Structural variant calling benchmarks

All structural variant benchmarks are done on the GIAB reference sample NA24385/HG002 using the [draft SV benchmark v1.1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/).  
The analyses were performed with the short-read and long-read single sample pipelines.

Sensitivity and positive predictive value (PPV) were measured using [Hap-Eval](https://github.com/Sentieon/hap-eval).

<table>
	<tr>
		<th>Test</th>
		<th>coverage</th>
		<th>sensitivity</th>
		<th>PPV</th>
	</tr>
	<tr>
		<td>short-read WGS - Manta 1.6.0</td>
		<td>39.5x</td>
		<td>36.30%</td>
		<td>96.53%</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.3.13</td>
		<td>39.5x</td>
		<td>49.90%</td>
		<td>97.68%</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy) - Sniffles 2.4</td>
		<td>40.5x</td>
		<td>90.59%</td>
		<td>98.03%</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy) - Sniffles 2.4</td>
		<td>40.5x</td>
		<td>91.02%</td>
		<td>98.05%</td>
	</tr>
</table>
