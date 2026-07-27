# megSAP benchmarks for germline pipelines

All benchmarks are perfomed with the megSAP release [2026_06](https://github.com/imgag/megSAP/releases/tag/2026_06).  
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
		<td>NovaSeq 6000 - 2x105bp PE</td>
		<td>101.1</td>
		<td>304.4</td>
	</tr>
	<tr>
		<td>short-read WGS</td>
		<td>Covaris</td>
		<td>Illumina TruSeq DNA PCR-Free</td>
		<td>NovaSeq 6000 - 2x159bp PE</td>
		<td>42.6</td>
		<td>377.3</td>
	</tr>
	<tr>
		<td>long-read WGS</td>
		<td>-</td>
		<td>Oxford Nanopore Tech. Ligation Sequencing Kit V14e (SQK-LSK114)</td>
		<td>PromethION P24</td>
		<td>44.3</td>
		<td>-</td>
	</tr>
</table>

The benchmarks were performed on the GIAB high-confidence region **with at least 15x coverage**:

<table>
	<tr>
		<th rowspan=2>Test</th>
		<th rowspan=2>%genom in high conf region and covered 15x</th>
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
		<td>88.29%</td>
		<td>0.9949</td>
        <td>0.9848</td>
        <td>0.9994</td>
        <td>0.9727</td>
        <td>0.9714</td>
        <td>0.9986</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.4</td>
		<td>88.29%</td>
		<td>0.9906</td>
        <td>0.9968</td>
        <td>0.9991</td>
        <td>0.9835</td>
        <td>0.9849</td>
        <td>0.9944</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>81.09%</td>
		<td>0.9963</td>
        <td>0.9992</td>
        <td>0.9998</td>
        <td>0.9914</td>
        <td>0.9971</td>
        <td>0.9991</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>81.09%</td>
		<td>0.9979</td>
		<td>0.9983</td>
		<td>0.9998</td>
		<td>0.9972</td>
		<td>0.9966</td>
		<td>0.9994</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>81.34%</td>
		<td>0.9989</td>
		<td>0.9987</td>
		<td>0.9999</td>
		<td>0.8764</td>
		<td>0.9456</td>
		<td>0.9859</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>81.34%</td>
		<td>0.9997</td>
		<td>0.9959</td>
		<td>0.9998</td>
		<td>0.9042</td>
		<td>0.9420</td>
		<td>0.9843</td>
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
		<td>0.9897</td>
        <td>0.9868</td>
        <td>0.9999</td>
        <td>0.9739</td>
        <td>0.9760</td>
        <td>1.0000</td>
	</tr>
	<tr>
		<td>short-read WES - DRAGEN 4.4</td>
		<td>0.9878</td>
        <td>0.9974</td>
        <td>0.9993</td>
        <td>0.9799</td>
        <td>0.9799</td>
        <td>0.9977</td>
	</tr>
	<tr>
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>0.9957</td>
        <td>0.9915</td>
        <td>0.9997</td>
        <td>0.9845</td>
        <td>0.9903</td>
        <td>0.9980</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>0.9964</td>
        <td>0.9973</td>
        <td>0.9996</td>
        <td>0.9884</td>
        <td>0.9771</td>
        <td>0.9980</td>
	</tr>
	<tr>
		<td>long-read WGS (high accuracy)</td>
		<td>0.9989</td>
        <td>0.9984</td>
        <td>0.9999</td>
        <td>0.9674</td>
        <td>0.9825</td>
        <td>0.9980</td>
	</tr>
	<tr>
		<td>long-read WGS (super accuracy)</td>
		<td>0.9998</td>
        <td>0.9928</td>
        <td>0.9999</td>
        <td>0.9655</td>
        <td>0.9618</td>
        <td>0.9980</td>
	</tr>
</table>

### CMRG benchmark

For genome sequening, we also performed the [CMRG benchmark](https://www.nature.com/articles/s41587-021-01158-1) based on the NA24385/HG002 sample.

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
		<td>short-read WGS</td>
		<td>Covaris</td>
		<td>Illumina TruSeq DNA PCR-Free</td>
		<td>NovaSeq 6000 - 2x159bp PE</td>
		<td>42.4</td>
		<td>348.3</td>
	</tr>
	<tr>
		<td>long-read WGS (ONT)</td>
		<td>-</td>
		<td>Oxford Nanopore Tech. Ligation Sequencing Kit V14e (SQK-LSK114)</td>
		<td>PromethION P24</td>
		<td>44.3</td>
		<td>-</td>
	</tr>
	<tr>
		<td>long-read WGS (PacBio)</td>
		<td>-</td>
		<td>LongRead Sequencing with PacBio (LR-PB-highcov)</td>
		<td>Revio</td>
		<td>38.0</td>
		<td>-</td>
	</tr>
</table>

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
		<td>short-read WGS - bwa-mem2, DeepVariant</td>
		<td>0.9803</td>
        <td>0.9969</td>
        <td>0.9981</td>
        <td>0.9333</td>
        <td>0.9526</td>
        <td>0.9957</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>0.9791</td>
        <td>0.9962</td>
        <td>0.9979</td>
        <td>0.9467</td>
        <td>0.9301</td>
        <td>0.9955</td>
	</tr>
	<tr>
		<td>ONT long-read WGS (high accuracy)</td>
		<td>0.9877</td>
        <td>0.9603</td>
        <td>0.9989</td>
        <td>0.7582</td>
        <td>0.8195</td>
        <td>0.9735</td>
	</tr>
	<tr>
		<td> ONT long-read WGS (super accuracy)</td>
		<td>0.9905</td>
        <td>0.9417</td>
        <td>0.9980</td>
        <td>0.7993</td>
        <td>0.8058</td>
        <td>0.9712</td>
	</tr>
	<tr>
		<td>PacBio long-read WGS </td>
		<td>0.9972</td>
        <td>0.9941</td>
        <td>0.9992</td>
        <td>0.9732</td>
        <td>0.9769</td>
        <td>0.9958</td>
	</tr>
</table>

## Structural variant calling benchmarks

All structural variant benchmarks are done on the GIAB reference sample NA24385/HG002 using the [draft SV benchmark v1.1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/).  
The analyses were performed with the short-read and long-read single sample pipelines.

samples: see [CMRG benchmark](#cmrg-benchmark)

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
		<td>44.55</td>
		<td>0.3649</td>
		<td>0.9656</td>
	</tr>
	<tr>
		<td>short-read WGS - DRAGEN 4.4</td>
		<td>44.22</td>
		<td>0.6453</td>
		<td>0.9474</td>
	</tr>
	<tr>
		<td>ONT long-read WGS (high accuracy) - Sniffles 2.4</td>
		<td>42.32</td>
		<td>0.8986</td>
		<td>0.9636</td>
	</tr>
	<tr>
		<td>ONT long-read WGS (super accuracy) - Sniffles 2.4</td>
		<td>42.36</td>
		<td>0.9027</td>
		<td>0.9648</td>
	</tr>
	<tr>
		<td>PacBio long-read WGS - Sniffles 2.4</td>
		<td>38.0</td>
		<td>0.9104</td>
		<td>0.9729</td>
	</tr>
</table>
