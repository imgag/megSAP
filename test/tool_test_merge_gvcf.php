<?php

require_once("framework.php");

//NOTE: test data generated from 21073LRa015_01, 21073LRa013_01 and 21073LRa014_01 in the regions chr[1-5]:1000000-1005000

$name = "merge_gvcf";
start_test($name);

########################## Trio ##################
$inputs = array(data_folder().$name."_in1.gvcf.gz", data_folder().$name."_in2.gvcf.gz", data_folder().$name."_in3.gvcf.gz");
# single-threaded
check_exec("php ".src_folder()."/NGS/{$name}.php -gvcfs ".implode(" ", $inputs)." -out ".output_folder().$name."_out1.vcf.gz -status affected control control -threads 1 --log ".output_folder().$name."_out1.log");
check_file(output_folder().$name."_out1.vcf.gz", data_folder().$name."_out1.vcf.gz");
check_file(output_folder().$name."_out1.gvcf.gz", data_folder().$name."_out1.gvcf.gz");
# multi-threaded
check_exec("php ".src_folder()."/NGS/{$name}.php -gvcfs ".implode(" ", $inputs)." -out ".output_folder().$name."_out2.vcf.gz -status affected control control -threads 4 --log ".output_folder().$name."_out2.log");
check_file(output_folder().$name."_out2.vcf.gz", data_folder().$name."_out1.vcf.gz");
check_file(output_folder().$name."_out2.gvcf.gz", data_folder().$name."_out1.gvcf.gz");

end_test();

?>


