# hg38 flavours benchmark

# Benchmark

To benchmark the small variant calling on different hg38 flavours, we performed tests with the sample NA24385/H002 processed with the Illumina WGS kit `TruSeq DNA PCR-Free`.  

The tests were performed with the the megSAP pipeline version 2023_11 - either with the default pipelein (BWA-mem2, freebayes 1.3.3) or with Dragen (version 4.2.4).  

These are the results on the default GiaB benchmark dataset: 
<table border=1>
<tr><th></th><th colspan=3 align=left>SNVs</th><th colspan=3 align=left>InDels</th></tr>
<tr><th></th><th align=left>sensitivity</th><th align=left>ppv</th><th align=left>genotyping</th><th align=left>sensitivity</th><th align=left>ppv</th><th align=left>genotyping</th></tr>
<tr><td>default genome of megSAP</td><td>0.9917</td><td>0.9935</td><td>0.9992</td><td>0.9732</td><td>0.9932</td><td>0.9833</td></tr>
<tr><td>no duplicate masking</td><td>0.9916</td><td>0.9935</td><td>0.9992</td><td>0.9733</td><td>0.9932</td><td>0.9833</td></tr>
<tr><td>duplicate masking with Illumina BED</td><td>0.9917</td><td>0.9935</td><td>0.9992</td><td>0.9733</td><td>0.9932</td><td>0.9833</td></tr>
<tr><td>no decoy</td><td>0.9917</td><td>0.9916</td><td>0.9992</td><td>0.9732</td><td>0.9924</td><td>0.9833</td></tr>
<tr><td>with ALT contigs, duplicate masking with Illumina BED</td><td>0.9917</td><td>0.9946</td><td>0.9992</td><td>0.9733</td><td>0.9936</td><td>0.9833</td></tr>
<tr><td>DRAGEN - default</td><td>0.9935</td><td>0.9981</td><td>0.9995</td><td>0.9955</td><td>0.9972</td><td>0.9993</td></tr>
<tr><td>DRAGEN - with ALT contigs, duplicate masking with Illumina BED</td><td>0.9939</td><td>0.9985</td><td>0.9994</td><td>0.9958</td><td>0.9974</td><td>0.9993</td></tr>
<tr><td>DRAGEN - graph</td><td>0.9976</td><td>0.9993</td><td>0.9992</td><td>0.9970</td><td>0.9978</td><td>0.9992</td></tr>
</table>

These are the results on the default GiaB CMRG dataset: 
<table border=1>
<tr><th></th><th colspan=3 align=left>SNVs</th><th colspan=3 align=left>InDels</th></tr>
<tr><th></th><th align=left>sensitivity</th><th align=left>ppv</th><th align=left>genotyping</th><th align=left>sensitivity</th><th align=left>ppv</th><th align=left>genotyping</th></tr>
<tr><td>default genome of megSAP</td><td>0.9646</td><td>0.9672</td><td>0.9935</td><td>0.8940</td><td>0.9344</td><td>0.9618</td></tr>
<tr><td>no duplicate masking</td><td>**0.9565**</td><td>0.967</td><td>0.9928</td><td>0.8898</td><td>0.9341</td><td>0.9619</td></tr>
<tr><td>duplicate masking with Illumina BED</td><td>0.9646</td><td>0.9672</td><td>0.9935</td><td>0.8940</td><td>0.9344</td><td>0.9618</td></tr>
<tr><td>no decoy</td><td>0.9655</td><td>0.9534</td><td>0.9935</td><td>0.8945</td><td>0.9274</td><td>0.9616</td></tr>
<tr><td>with ALT contigs, duplicate masking with Illumina BED</td><td>0.9646</td><td>0.9672</td><td>0.9935</td><td>0.8940</td><td>0.9344</td><td>0.9618</td></tr>
<tr><td>DRAGEN - default</td><td>0.9604</td><td>0.9920</td><td>0.9966</td><td>0.9323</td><td>0.9248</td><td>0.9930</td></tr>
<tr><td>DRAGEN - with ALT contigs, duplicate masking with Illumina BED</td><td>0.9614</td><td>0.9909</td><td>0.9965</td><td>0.9328</td><td>0.9227</td><td>0.9930</td></tr>
<tr><td>DRAGEN - graph</td><td>0.9705</td><td>0.9944</td><td>0.9957</td><td>0.9357</td><td>0.9243</td><td>0.9936</td></tr>
</table>

# Conclusion

For the default pipline, it makes no real difference which genome flavour is used as long as false ducplications are masked.

For DRAGEN, we could improve the variant calling by using the graph genome, instead of the megSAP default genome.

# Notes

Dataset used: NA24385_03  
Folder: /mnt/storage2/users/ahsturm1/scripts/2023_10_18_HG38_flavours/  
