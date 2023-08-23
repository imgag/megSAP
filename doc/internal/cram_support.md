# CRAM support

To save storage space, megSAP now uses CRAM instead of BAM files as result of the mapping pipeline for DNA.  
This page show some benchmarking results between BAM and CRAM.


## size comparison

CRAM files are only 30-40% the size of BAMs:

<table>
	<tr>
		<th>sample</th>
		<th>BAM [GB]</th>
		<th>CRAM [GB]</th>
	</tr>
	<tr>
		<th>NA12878x2_93 (WGS @ 40x, NovaSeq6000)</th>
		<th>51.95</th>
		<th>19.99</th>
	</tr>
	<tr>
		<th>NA12878_58 (Agilent V7 @ 142x, NovaSeq6000)</th>
		<th>6.13</th>
		<th>2.07</th>
	</tr>
	<tr>
		<th>NA12878_58 (Twist Custom V2 @ 113, NovaSeq6000)</th>
		<th>5.95</th>
		<th>2.00</th>
	</tr>
</table>

## conversion benchmarks

### BAM to CRAM

BAM to CRAM  conversion for a genome takes about an hour, but scales linear with the number of threads:

<table>
	<tr>
		<th>sample</th>
		<th>runtime [min]</th>
		<th>memory [GB]</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 1thread</th>
		<th>54.71</th>
		<th>277.00</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 2thread</th>
		<th>26.43</th>
		<th>770.23</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 4thread</th>
		<th>13.77</th>
		<th>837.60</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 8thread</th>
		<th>7.79</th>
		<th>927.748</th>
	</tr>
</table>

### CRAM to BAM

CRAM to BAM conversion for a genome takes more than an hour, but scales linear with the number of threads:

<table>
	<tr>
		<th>sample</th>
		<th>runtime [min]</th>
		<th>memory [GB]</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 1thread</th>
		<th>86.35</th>
		<th>277.00</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 2thread</th>
		<th>43.82</th>
		<th>770.23</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 4thread</th>
		<th>22.84</th>
		<th>837.60</th>
	</tr>
	<tr>
		<th>NA12878x2_93 @ 8thread</th>
		<th>12.49</th>
		<th>927.748</th>
	</tr>
</table>


## variant calling benchmarks

To give an idea about the runtime difference between BAM and CRAM we performed typical variant calling taks on the same WES dataset.

### freebayes (small variant calling)

The runtime of freebayes increases by about 100%:

<table>
	<tr>
		<th>format</th>
		<th>runtime [min]</th>
		<th>memory [GB]</th>
	</tr>
	<tr>
		<th>BAM</th>
		<th>34.83</th>
		<th>788.39</th>
	</tr>
	<tr>
		<th>CRAM</th>
		<th>64.57</th>
		<th>795.23</th>
	</tr>
</table>

### BedCoverage (CNV calling)

The runtime of BedCoverage increases by about 300%:

<table>
	<tr>
		<th>format</th>
		<th>runtime [min]</th>
		<th>memory [GB]</th>
	</tr>
	<tr>
		<th>BAM</th>
		<th>9.52</th>
		<th>66.13</th>
	</tr>
	<tr>
		<th>CRAM</th>
		<th>41.25</th>
		<th>186.96</th>
	</tr>
</table>

### Manta (SV calling)


The runtime of Manta increases by about 30%:

<table>
	<tr>
		<th>format</th>
		<th>runtime [min]</th>
		<th>memory [GB]</th>
	</tr>
	<tr>
		<th>BAM</th>
		<th>13.11</th>
		<th>58.58</th>
	</tr>
	<tr>
		<th>CRAM</th>
		<th>16.19</th>
		<th>220.11</th>
	</tr>
</table>

Generally, the runtime and memory consumption of the read processing increases when using CRAM.  
Thus, we will perform the processing on BAM files and use CRAM as storage format, i.e. when the data is not analyzed. 
