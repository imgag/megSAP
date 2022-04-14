<?php

require_once("framework.php");

$name = "an_spliceai";
start_test($name);

//no variants left to score after filtering > input is equal to output
$out1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out1 --log ".output_folder().$name."_out1.log");
check_file($out1, data_folder().$name."_in1.vcf", true);

//too many variants to score => input is equal to output
$out2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -max_vars 0 -out $out2 --log ".output_folder().$name."_out2.log");
check_file($out2, data_folder().$name."_in2.vcf", true);

//perform scoring (no contig headers)
$out3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -out $out3 --log ".output_folder().$name."_out3.log");
check_file($out3, data_folder().$name."_out3.vcf", true);

end_test();

?>
