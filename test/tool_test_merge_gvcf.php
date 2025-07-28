<?php

require_once("framework.php");

//NOTE: test data generated from 21073LRa015_01, 21073LRa013_01 and 21073LRa014_01 in the regions chr[1-5]:1000000-1005000

$name = "merge_gvcf";
start_test($name);

# longread - single-threaded
$inputs_lr = array(data_folder().$name."_in1.gvcf.gz", data_folder().$name."_in2.gvcf.gz", data_folder().$name."_in3.gvcf.gz");
check_exec("php ".src_folder()."/Tools/{$name}.php -gvcfs ".implode(" ", $inputs_lr)." -out ".output_folder().$name."_out1.vcf.gz -status affected control control -threads 1 -mode longread --log ".output_folder().$name."_out1.log");
check_file(output_folder().$name."_out1.vcf.gz", data_folder().$name."_out1.vcf.gz");
check_file(output_folder().$name."_out1.gvcf.gz", data_folder().$name."_out1.gvcf.gz");

# longread - multi-threaded
check_exec("php ".src_folder()."/Tools/{$name}.php -gvcfs ".implode(" ", $inputs_lr)." -out ".output_folder().$name."_out2.vcf.gz -status affected control control -threads 4 -mode longread --log ".output_folder().$name."_out2.log");
check_file(output_folder().$name."_out2.vcf.gz", data_folder().$name."_out1.vcf.gz");
check_file(output_folder().$name."_out2.gvcf.gz", data_folder().$name."_out1.gvcf.gz");

# dragen - single-threaded
$inputs_dragen = array(data_folder().$name."_dragen_in1.gvcf.gz", data_folder().$name."_dragen_in2.gvcf.gz");
check_exec("php ".src_folder()."/Tools/{$name}.php -gvcfs ".implode(" ", $inputs_dragen)." -out ".output_folder().$name."_out3.vcf.gz -status affected control -threads 1 -mode dragen --log ".output_folder().$name."_out3.log");
check_file(output_folder().$name."_out3.vcf.gz", data_folder().$name."_out3.vcf.gz");

end_test();

?>


