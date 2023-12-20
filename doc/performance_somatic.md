# megSAP benchmarks

## tumor-normal short-read pipline

All performance benchmarks are performed on the GIAB reference samples NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/) and the NA12877 using the [PlatimunGenomes variant list](https://github.com/Illumina/PlatinumGenomes). The two samples were mixed to prepare a in-silico tumor sample with a specified tumor content of 5, 10, 20 or 40 %.
The samples were mapped with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1) and the calling was done using [short-read tumor normal pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/somatic_dna.php).

Sensitivity and positive predictive value (PPV) were measured using our [somatic validation tool](https://github.com/imgag/megSAP/blob/master/src/Tools/validate_somatic.php).

### Whole exome sequencing

The WES samples were processed with a custom exome kit based on a Twist enrichment (Core, RefSeq, Mito and custom content) and sequenced on NovaSeq X Plus using 151PE.
The normal samples had a depth of 110x while the mixed 'tumor' samples had a depth of 220x.
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>90.37%</td>
		<td>99.21%</td>
		<td>45.58%</td>
		<td>89.10%</td>
    		<td>88.32%</td>
		<td>98.95%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>96.82%</td>
		<td>99.79%</td>
		<td>69.42%</td>
		<td>97.83%</td>
    		<td>95.57%</td>
		<td>99.73%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>98.55%</td>
		<td>99.98%</td>
		<td>77.12%</td>
		<td>99.50%</td>
    		<td>97.57%</td>
		<td>99.96%</td>
	</tr>
</table>

<table>
	<tr>
		<th rowspan=2>Test - Dragen mapping + Strelka2 calling</th>
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>90.12%</td>
		<td>99.49%</td>
		<td>45.77%</td>
		<td>90.15%</td>
    <td>88.09%</td>
		<td>99.25%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>96.61%</td>
		<td>99.84%</td>
		<td>69.42%</td>
		<td>97.57%</td>
    <td>95.37%</td>
		<td>99.76%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>98.34%</td>
		<td>99.98%</td>
		<td>76.92%</td>
		<td>99.50%</td>
    <td>97.36%</td>
		<td>99.96%</td>
	</tr>
</table>


#### Dragen calling
<table>
	<tr>
		<th rowspan=2>Test - BWA-MEM2 + Dragen calling</th>
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>93.78%</td>
		<td>98.61%</td>
		<td>60.39%</td>
		<td>95.44%</td>
    		<td>92.26%</td>
		<td>98.51%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>98.29%</td>
		<td>99.66%</td>
		<td>75.19%</td>
		<td>98.24%</td>
    		<td>97.23%</td>
		<td>99.61%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>99.07%</td>
		<td>99.99%</td>
		<td>88.27%</td>
		<td>99.78%</td>
    		<td>98.58%</td>
		<td>99.98%</td>
	</tr>
</table>

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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>94.14%</td>
		<td>98.84%</td>
		<td>63.01%</td>
		<td>92.14%</td>
    <td>92.72%</td>
		<td>98.62%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>98.32%</td>
		<td>99.74%</td>
		<td>79.04%</td>
		<td>97.39%</td>
    <td>97.44%</td>
		<td>99.65%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_22_NA12877_23 --->
		<td>99.05%</td>
		<td>99.99%</td>
		<td>91.73%</td>
		<td>99.38%</td>
    <td>98.72%</td>
		<td>99.96%</td>
	</tr>
</table> 

### Conclusion

While the mapping has only a small influence on the sensitivity and precision. The dragen calling improves the sensitivity especially for variants with low allel frequency.
The sensitivity for InDel variants with Dragen calling improves over Strelka2 by over 10% and for SNVs by ~2%.   
