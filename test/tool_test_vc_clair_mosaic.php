<?php

require_once("framework.php");

//NOTE: test data generated from 21073LRa015 and 23014LRa028_01 using the region defined in BED file vc_clair_in_roi.bed

$name = "vc_clair_mosaic";
start_test($name);


########################## BAM ##################
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in1.bam -out ".output_folder()."vc_clair_mosaic_out1.vcf.gz -target ".data_folder().$name."_in_roi.bed -name LSK114_01 -model ont_r10_dorado_sup_5khz -threads 4 --log ".output_folder().$name."_out1.log");
check_file(output_folder()."vc_clair_mosaic_out1.vcf.gz", data_folder().$name."_out1.vcf.gz");

########################## SQK-LSK109 ##################
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in1.cram -out ".output_folder()."vc_clair_mosaic_out2.vcf.gz -target ".data_folder().$name."_in_roi.bed -name LSK114_01 -model ont_r10_dorado_sup_5khz -threads 4 --log ".output_folder().$name."_out1.log");
check_file(output_folder()."vc_clair_mosaic_out2.vcf.gz", data_folder().$name."_out1.vcf.gz");

end_test();

?>


