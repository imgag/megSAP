<?php

require_once("framework.php");

$name = "baf_germline";
start_test($name);

$bam = data_folder()."vc_strelka2_tu_in.bam";
$vcf = data_folder()."/{$name}_snps.vcf.gz";

// test1: default parameters
$out2 = output_folder()."{$name}_out1.igv";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam $bam -out $out2 -vcf $vcf --log ".output_folder()."{$name}_1.log");
check_file($out2, data_folder()."{$name}_out1.igv");

// test2: with downsampling
$out3 = output_folder()."{$name}_out2.igv";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam $bam -out $out3 -vcf $vcf -downsample 10 --log ".output_folder()."{$name}_2.log");
check_file($out3, data_folder()."{$name}_out2.igv");

end_test();

?>