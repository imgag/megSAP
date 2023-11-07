<?php

require_once("framework.php");

$name = "vc_expansionhunter";
start_test($name);

//test
$out_file_vcf = output_folder().$name."_out1.vcf";
$out_file_json = output_folder().$name."_out1.json";
$out_file_svg_prefix = output_folder()."repeat_expansions/{$name}_out1_";
check_exec("php ".src_folder()."/NGS/".$name.".php -in ".data_folder().$name."_in1.bam -out $out_file_vcf");
remove_lines_containing($out_file_vcf, array("##filedate=", "##reference="));
check_file($out_file_vcf, data_folder().$name."_out1.vcf");
check_file($out_file_json, data_folder().$name."_out1.json");
foreach(array("AFF2", "ATXN10", "FXN", "TCF4") as $re)
{
    check_file($out_file_svg_prefix.$re.".svg", data_folder().$name."_out1_".$re.".svg");
}

end_test();

?>