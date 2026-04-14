<?php

require_once("framework.php");

$name = "vc_deepsomatic";
start_test($name);

$t_bam = data_folder().$name."_tu_in.bam";
$n_bam = data_folder().$name."_no_in.bam";
$debug_region = data_folder().$name."_debug_region.bed";

########################## tumor-normal with target ##########################

$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -debug_region {$debug_region} -bam_tumor {$t_bam} -out $out_file1 -model_type WES -bam_normal {$n_bam} -threads 4 -target ".data_folder().$name.".bed --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

########################## tumor-normal raw output ##########################

$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -debug_region {$debug_region} -bam_tumor {$t_bam} -out $out_file2 -model_type WES -bam_normal {$n_bam} -threads 4 -target ".data_folder().$name.".bed -raw_output --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

########################## tumor-only ##########################

$out_file3 = output_folder().$name."_out3.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -debug_region {$debug_region} -bam_tumor {$t_bam} -out $out_file3 -model_type WGS_TUMOR_ONLY -tumor_only -threads 4 -target ".data_folder().$name.".bed --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.vcf.gz");

########################## tumor-normal no target ##########################

$out_file4 = output_folder().$name."_out4.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -debug_region {$debug_region} -bam_tumor {$t_bam} -out $out_file4 -model_type WES -bam_normal {$n_bam} -threads 4 --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.vcf.gz");

end_test();

?>


