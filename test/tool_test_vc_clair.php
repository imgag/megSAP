<?php

require_once("framework.php");

//NOTE: test data generated from 21073LRa015 and 23014LRa028_01 using the region defined in BED file vc_clair_in_roi.bed

$name = "vc_clair";
start_test($name);


########################## SQK-LSK114 ##################
$model1 = get_path("clair3_models")."/r1041_e82_400bps_hac_g632/";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in1.bam -folder ".output_folder()." -target ".data_folder().$name."_in_roi.bed -name LSK114_01 -model {$model1} -threads 4 --log ".output_folder().$name."_out1.log");
check_file(output_folder()."LSK114_01_var.vcf.gz", data_folder().$name."_out1.vcf.gz");

########################## SQK-LSK109 ##################
$model2 = get_path("clair3_models")."/r941_prom_sup_g5014/";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in2.bam -folder ".output_folder()." -target ".data_folder().$name."_in_roi.bed -name LSK109_01 -model {$model2} -threads 4 --log ".output_folder().$name."_out2.log");
check_file(output_folder()."LSK109_01_var.vcf.gz", data_folder().$name."_out2.vcf.gz");

########################## SQK-LSK114 mito-calling ##################
$model1 = get_path("clair3_models")."/r1041_e82_400bps_hac_g632/";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in3.bam -folder ".output_folder()." -target ".data_folder().$name."_in_roi.bed -name LSK114_02 -model {$model1} -threads 4 --log ".output_folder().$name."_out3.log");
# remove unreliably varian chrMT:7337 G>A
exec("zcat ".output_folder()."LSK114_02_var.vcf.gz | egrep -v '^chrMT\t7337\t' | bgzip -c > ".output_folder()."LSK114_02_var_modified.vcf.gz");
check_file(output_folder()."LSK114_02_var_modified.vcf.gz", data_folder().$name."_out3.vcf.gz");

end_test();

?>


