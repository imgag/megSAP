<?php

require_once("framework.php");

$name = "an_filter_dragen_somatic";
start_test($name);

$vcf = data_folder() . "an_filter_dragen_somatic_in1.vcf.gz";
$tumor_name  = "DNA24XXXX1A1_01";
$normal_name = "DNA24XXXX2A1_01";
$out = output_folder() . "an_filter_dragen_somatic_out1.vcf.gz"; 
$target = data_folder() . "an_filter_dragen_somatic_target1.bed"; 
check_exec("php ".src_folder()."/Tools/an_filter_dragen_somatic.php -in $vcf -tumor_name $tumor_name -normal_name $normal_name -out $out -target $target");
check_file($out, data_folder() . "an_filter_dragen_somatic_ref1.vcf.gz", false);


end_test();
?>
