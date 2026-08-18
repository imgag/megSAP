# megSAP benchmarks for somatic pipelines

## tumor-normal short-read pipline

All performance benchmarks are performed on the GIAB reference samples NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/) and the NA12877 using the [PlatimunGenomes variant list](https://github.com/Illumina/PlatinumGenomes). The two samples were mixed to prepare a in-silico tumor sample with a specified tumor content of 5, 10, 20 or 40 %.
The samples were mapped with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1) and the calling was done using [short-read tumor normal pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/somatic_tumor_normal.php).

Sensitivity and positive predictive value (PPV) were measured using our [somatic validation tool](https://github.com/imgag/megSAP/blob/master/src/Auxilary/validate_somatic.php).

### Whole exome sequencing

The WES samples were processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content) and sequenced on NovaSeq X Plus using 151PE.
The normal sample had a depth of 110x while the mixed 'tumor' samples had a depth of 180-215x.
The benchmarks were performed on GIAB / PlatinumGenomes high-confidence regions with at least 60x coverage.

#### Strelka2 calling
<table>
	<tr>
		<th rowspan=2>Test - BWA-MEM2 + Strelka2 calling</th>
		<th colspan=2>SNVs</th>
		<th colspan=2>InDels</th>
    	<th colspan=2>SNVs+InDels</th>
	</tr>
	<tr>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>sensitivity</th>
		<th>PPV</th>
    	<th>sensitivity</th>
		<th>PPV</th>
	</tr>
	<tr>
		<td>Variants >= 5% allele freq</td> 
		<td>92.83%</td>
		<td>98.89%</td>
		<td>52.16%</td>
		<td>81.45%</td>
    	<td>90.72%</td>
		<td>98.27%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td>
		<td>97.53%</td>
		<td>99.70%</td>
		<td>67.28%</td>
		<td>97.32%</td>
    	<td>95.97%</td>
		<td>99.61%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td>
		<td>98.52%</td>
		<td>99.98%</td>
		<td>74.54%</td>
		<td>99.79%</td>
    	<td>97.28%</td>
		<td>99.98%</td>
	</tr>
</table>


#### Dragen 4.4 calling

<table>
	<tr>
		<th rowspan=2>Test - Dragen mapping + Dragen calling</th>
		<th colspan=2>SNVs</th>
		<th colspan=2>InDels</th>
    	<th colspan=2>SNVs+InDels</th>
	</tr>
	<tr>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>sensitivity</th>
		<th>PPV</th>
    	<th>sensitivity</th>
		<th>PPV</th>
	</tr>
	<tr>
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>94.70%</td>
		<td>98.89%</td>
		<td>58.80%</td>
		<td>92.48%</td>
    	<td>92.83%</td>
		<td>98.66%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>98.62%</td>
		<td>99.69%</td>
		<td>73.30%</td>
		<td>96.35%</td>
    	<td>97.30%</td>
		<td>99.56%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>99.20%</td>
		<td>99.97%</td>
		<td>85.96%</td>
		<td>98.93%</td>
    	<td>98.51%</td>
		<td>99.93%</td>
	</tr>
</table> 

### Conclusion

The comparision between the open-source calling with bwa-mem2 + strelka2 and the illumina dragen shows the dragen calling provides better results.
The SNV sensitivity increases by about 1% but the main difference is the INDEL calling where the dragen calling hast 5-10% more sensitivity while keeping a high prescision especially for low allele frequency variants.


### Whole genome sequencing

The WGS samples were processed with the "NEBNext UltraShear FFPE DNA Library Prep Kit" from "New England Biolabs" and sequenced on NovaSeq X Plus using 151PE.
The normal sample had a depth of 30x while the mixed 'tumor' samples had a depth of 100x.
The benchmarks were performed on GIAB / PlatinumGenomes high-confidence regions with at least 60x coverage.

Note: bwa-mem2 + Strelka2 was not validated as the runtime for the Strelka2 calling is to long to be practical given our environement.

#### Dragen 4.4 calling
<table>
	<tr>
		<th rowspan=2>Test - Dragen mapping + Dragen calling</th>
		<th colspan=2>SNVs</th>
		<th colspan=2>InDels</th>
    	<th colspan=2>SNVs+InDels</th>
	</tr>
	<tr>
		<th>sensitivity</th>
		<th>PPV</th>
		<th>sensitivity</th>
		<th>PPV</th>
    	<th>sensitivity</th>
		<th>PPV</th>
	</tr>
	<tr>
		<td>Variants >= 5% allele freq</td> 
		<td>69.05%</td>
		<td>98.83%</td>
		<td>32.90%</td>
		<td>97.11%</td>
    	<td>64.81%</td>
		<td>98.72%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td>
		<td>96.00%</td>
		<td>99.65%</td>
		<td>64.62%</td>
		<td>99.57%</td>
    	<td>92.32%</td>
		<td>99.64%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td>
		<td>99.34%</td>
		<td>99.96%</td>
		<td>88.89%</td>
		<td>99.95%</td>
    	<td>98.12%</td>
		<td>99.96%</td>
	</tr>
</table>

