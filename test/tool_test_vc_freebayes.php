<?php

require_once("framework.php");

$name = "vc_freebayes";
start_test($name);

//test with target region - AF>=15%
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file1 -target ".data_folder().$name."_in.bed --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//test without target data - AF>=20%
$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file2 --log ".output_folder().$name."_out2.log -min_af 0.20");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

end_test();

?>
