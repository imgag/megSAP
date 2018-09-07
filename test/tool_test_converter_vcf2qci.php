<?php

require_once("framework.php");

$name = "converter_vcf2qci";
start_test($name);

$in_file = data_folder().$name."_in.vcf.gz";
$out_file = output_folder().$name."_out.vcf.gz";

//check with tumor-normal
check_exec("php ".src_folder()."/Tools/converter_vcf2qci.php -in $in_file -t_id vc_strelka_tu_in -n_id vc_strelka_no_in -out $out_file");
check_file($out_file, data_folder().$name."_out.vcf.gz");

end_test();

?>