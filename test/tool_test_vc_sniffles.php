<?php

require_once("framework.php");

$name = "vc_sniffles";
start_test($name);

//NOTE: test data generated from 21073LRa015 and 23014LRa028_01 using the region defined in BED file vc_clair_in_roi.bed

########################## germline ##################
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -include_mosaic -bam ".data_folder().$name."_in1.bam -target ".data_folder().$name."_in1.bed -out $out_file1 -threads 2 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

########################## Trio ##################
//NOTE: test data generated from 21073LRa01[3|4|5] in region chr1:1000000-1100000
$out_file2 = output_folder().$name."_out2.vcf.gz";
$input_files = array();
$input_files[] = data_folder().$name."_in2_index.bam";
$input_files[] = data_folder().$name."_in2_father.bam";
$input_files[] = data_folder().$name."_in2_mother.bam";
check_exec("php ".src_folder()."/Tools/{$name}.php -include_mosaic -bam ".implode(" ", $input_files)." -sample_ids 21073LRa015 21073LRa013 21073LRa014 -target ".data_folder().$name."_in1.bed -out $out_file2 -threads 2 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

//TODO
########################## somatic ##################
// $out_file2 = output_folder().$name."_out2.vcf.gz";
// check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in2.bam -out $out_file2 -somatic -threads 2 --log ".output_folder().$name."_out2.log");
// check_file($out_file2, data_folder().$name."_out2.vcf.gz");

end_test();

?>


