<?php

require_once("framework.php");

$name = "an_somatic_cancerhotspots";

start_test($name);
$out1 = output_folder()."{$name}_out1.vcf";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder()."{$name}_in1.vcf -out {$out1}");
check_file($out1, data_folder().$name."_out1.vcf");
end_test();

?>
