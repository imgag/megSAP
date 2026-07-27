# megSAP benchmarks for somatic pipelines

## tumor-normal short-read pipline

All performance benchmarks are performed on the GIAB reference samples NA12878 using the [gold-standard variant list v4.2.1](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/) and the NA12877 using the [PlatimunGenomes variant list](https://github.com/Illumina/PlatinumGenomes). The two samples were mixed to prepare a in-silico tumor sample with a specified tumor content of 5, 10, 20 or 40 %.
The samples were mapped with the [short-read single sample pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/analyze.php) on the GRCh38 reference genome with [masked false duplications](https://www.nature.com/articles/s41587-021-01158-1) and the calling was done using [short-read tumor normal pipeline](https://github.com/imgag/megSAP/blob/master/src/Pipelines/somatic_tumor_normal.php).

Sensitivity and positive predictive value (PPV) were measured using our [somatic validation tool](https://github.com/imgag/megSAP/blob/master/src/Auxilary/validate_somatic.php).

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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>93.37%</td>
		<td>98.86%</td>
		<td>51.08%</td>
		<td>76.50%</td>
    		<td>90.91%</td>
		<td>98.02%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>97.56%</td>
		<td>99.75%</td>
		<td>68.00%</td>
		<td>99.75%</td>
    		<td>96.02%</td>
		<td>99.55%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>98.58%</td>
		<td>99.98%</td>
		<td>74.46%</td>
		<td>99.79%</td>
    		<td>97.33%</td>
		<td>99.98%</td>
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>92.97%</td>
		<td>99.19%</td>
		<td>47.31%</td>
		<td>82.35%</td>
    <td>90.59%</td>
		<td>98.64%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>97.43%</td>
		<td>99.84%</td>
		<td>67.28%</td>
		<td>95.42%</td>
    <td>95.86%</td>
		<td>99.67%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>98.41%</td>
		<td>99.99%</td>
		<td>74.65%</td>
		<td>100.00%</td>
    <td>97.18%</td>
		<td>99.99%</td>
	</tr>
</table>


#### Dragen 4.4 calling
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>94.75%</td>
		<td>98.38%</td>
		<td>55.08%</td>
		<td>93.47%</td>
    		<td>92.75%</td>
		<td>98.22%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>98.67%</td>
		<td>99.67%</td>
		<td>69.38%</td>
		<td>97.83%</td>
    		<td>97.15%</td>
		<td>99.60%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>99.23%</td>
		<td>99.97%</td>
		<td>80.62%</td>
		<td>99.62%</td>
    		<td>98.26%</td>
		<td>99.96%</td>
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
		<td>Variants >= 5% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>94.72%</td>
		<td>98.79%</td>
		<td>57.14%</td>
		<td>90.73%</td>
    <td>92.77%</td>
		<td>98.51%</td>
	</tr>
	<tr>
		<td>Variants >= 10% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>98.59%</td>
		<td>99.73%</td>
		<td>73.27%</td>
		<td>96.36%</td>
    <td>97.27%</td>
		<td>99.59%</td>
	</tr>
	<tr>
		<td>Variants >= 20% allele freq</td> <!--- dataset: NA12878x3_68_NA12877_46 --->
		<td>99.14%</td>
		<td>99.97%</td>
		<td>84.95%</td>
		<td>99.10%</td>
    <td>98.40%</td>
		<td>99.93%</td>
	</tr>
</table> 

### Conclusion

While the mapping has only a small influence on the sensitivity and precision. The dragen calling improves the sensitivity especially for variants with low allel frequency.
The sensitivity for InDel variants with Dragen calling improves over Strelka2 by over 5% and for SNVs by ~2%.   
