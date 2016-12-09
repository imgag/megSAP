<?php

require_once("framework.php");

$name = "converter_vcf2qci";
start_test($name);

$in_file = data_folder().$name."_in.vcf.gz";
$out_file = output_folder().$name."_out.vcf.gz";

check_exec("php ".src_folder()."/Tools/converter_vcf2qci.php -in '$in_file' -out '$out_file'");
check_file($out_file, data_folder().$name."_out.vcf.gz");

end_test();

?>