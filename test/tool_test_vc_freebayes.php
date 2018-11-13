<?php

require_once("framework.php");

$name = "vc_freebayes";
start_test($name);

//test with target region
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file1 -target ".data_folder().$name."_in.bed --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//test without target region - AF>=20%
$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file2 --log ".output_folder().$name."_out2.log -min_af 0.20");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

//test with target region - target_extend
$out_file3 = output_folder().$name."_out3.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file3 -target ".data_folder().$name."_in.bed -target_extend 50 --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.vcf.gz");

end_test();

?>
