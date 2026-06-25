# megSAP benchmarks for germline pipelines

All benchmarks are perfomed with the megSAP release [2025_03](https://github.com/imgag/megSAP/releases/tag/2025_03).  
As reference genome GRCh38 with decoy chromosomes, without ALT chromosomes and with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1) was used.

## Small variants benchmark

All small variant benchmarks are done on the GIAB reference sample NA12878/HG001 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/).  
The analyses were performed with the short-read and long-read single sample pipelines.

Sensitivity, positive predictive value (PPV) and genotyping accuracy were measured using our [validation tool](https://github.com/imgag/megSAP/blob/master/src/Auxilary/validate_small_variants.php).

<!--- TODO: sample names --->

The following data was used for the benchmark:

<table>
	<tr>
		<th>Type</th>
		<th>DNA Fragmentation</th>
		<th>Kit</th>
		<th>Sequencer</th>
		<th>Mean depth</th>
		<th>Mean insert size</th>
	</tr>
	<tr>
		<td>short-read WES</td>
		<td>Covaris</td>
		<td>Twist custom exome kit (Core, RefSeq, Mito and custom content)</td>
		<td>TODO  - TODO PE</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS</td>
		<td>Covaris</td>
		<td>Illumina TruSeq DNA PCR-Free</td>
		<td>TODO  - TODO PE</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS</td>
		<td>-</td>
		<td>Oxford Nanopore Tech. Ligation Sequencing Kit V14e (SQK-LSK114)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>-</td>
	</tr>
</table>

The benchmarks were performed on the GIAB high-confidence region **with at least 15x coverage**:

<table>
	<tr>
		<th rowspan=2>Test</th>
		<th rowspan=2>%roi covered 15x</th>
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
		<td>short-read WES - bwa-mem2, DeepVariant</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4/td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
</table>

	
## Small variants benchmark - coding region

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
		<td>short-read WES - bwa-mem2, DeepVariant</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
</table>

### CMRG benchmark

For genome sequening, we also performed the [CMRG benchmark](https://www.nature.com/articles/s41587-021-01158-1) based on the NA24385/HG002 sample.

TODO: sample table

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
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
</table>

## Structural variant calling benchmarks

All structural variant benchmarks are done on the GIAB reference sample NA24385/HG002 using the [draft SV benchmark v1.1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/).  
The analyses were performed with the short-read and long-read single sample pipelines.

TODO: sample table or reference to table in CMRG

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
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy) - Sniffles 2.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy) - Sniffles 2.4</td>
		<td>TODO</td>
		<td>TODO</td>
		<td>TODO</td>
	</tr>
</table>
