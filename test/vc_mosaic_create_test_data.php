<?php

// "Input"
$ngs_bits = "/mnt/users/ahott1a1/ngs-bits/bin/";

$bam = "/mnt/storage2/projects/test/WES_Onko_test/Sample_NA12877_NA12878x2_10perc/NA12877_NA12878x2_10perc.bam";
$vcf = "/mnt/storage2/projects/test/WES_Onko_test/Sample_NA12877_14/NA12877_14_var.vcf.gz";
$genes = "/mnt/users/ahott1a1/megSAP/test/data/vc_mosaic_in_genes.txt";


//tmp files:
$gene_regions = "data/vc_mosaic_in_genes.bed";
$bed_ext = "data/vc_mosaic_in_genes_ext.bed";
$bed_merged = "data/vc_mosaic_in_genes_ext_merged.bed";
 

//target region:
exec("$ngs_bits/GenesToBed -in $genes -source ccds -mode exon -fallback -out $gene_regions");
exec("$ngs_bits/BedExtend -in $gene_regions -n 20 -out $bed_ext"); 
exec("$ngs_bits/BedMerge -in $bed_ext -out $bed_merged"); 

//limit and downsample bam:
$threads = 5;
$percentage = 0.5;

$limited_bam = "data/vc_mosaic_tmp_limited.bam";
$input_bam = "data/vc_mosaic_in1.bam";
exec("samtools view -b -L $bed_merged -@ $threads $bam > $input_bam");
exec("samtools index -b $input_bam");


//limit base vcf:
$filtered_vcf = "data/vc_mosaic_in1.vcf";
exec("gunzip -c $vcf | $ngs_bits/VcfFilter -reg $bed_ext > $filtered_vcf");

//delete tmp files:
exec("rm $gene_regions $bed_ext")

?>