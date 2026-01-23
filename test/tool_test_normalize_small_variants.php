<?php

require_once("framework.php");

$name = "normalize_small_variants";
start_test($name);

//DeepVariant input
$in_file1 = data_folder().$name."_in1.vcf";
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/Tools/{$name}.php -in {$in_file1} -out {$out_file1}");
check_file($out_file1, data_folder().$name."_out1.vcf");

end_test();

?>
