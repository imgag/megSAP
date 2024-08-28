<?php

require_once("framework.php");

//NOTE: test data generated from 21073LRa015 and 23014LRa028_01 using the region defined in BED file vc_clair_in_roi.bed

$name = "vc_clair_singularity";
start_test($name);

########################## SQK-LSK109 ##################
$model2 = "r941_prom_sup_g5014";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in2.bam -folder ".output_folder()." -target ".data_folder().$name."_in_roi.bed -name LSK109_01 -model {$model2} -threads 4 --log ".output_folder().$name."_out2.log");
check_file(output_folder()."LSK109_01_var.vcf.gz", data_folder().$name."_out2.vcf.gz");

end_test();

?>


