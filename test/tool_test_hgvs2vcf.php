<?php

require_once("framework.php");

$name = "hgvs2vcf";
start_test($name);

$out_file1 = output_folder()."/".$name."_out1.vcf";
$output = check_exec("php ".realpath(src_folder())."/Tools/".$name.".php -in ".data_folder()."/".$name."_in1.txt -out $out_file1");
check_file($out_file1, data_folder().$name."_out1.vcf", true);

end_test();

?>