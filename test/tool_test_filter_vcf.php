<?php

require_once("framework.php");

$name = "filter_vcf";

start_test("filter_vcf");

$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in1.vcf -t_id TUMOR -n_id NORMAL -out $out_file1");
check_file($out_file1, data_folder().$name."_out1.vcf");

end_test();
?>