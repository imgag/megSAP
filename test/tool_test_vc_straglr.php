<?php

require_once("framework.php");

$name = "vc_straglr";
start_test($name);

//test
$out_file_vcf = output_folder().$name."_out1.vcf";
$out_file_bed = output_folder().$name."_out1.bed";
$out_file_log = output_folder().$name."_out1.log";
$out_file_svg_prefix = output_folder()."repeat_expansions/{$name}_out1_";
$variant_catalog = data_folder().$name."_variant_catalog_in1.bed";
check_exec("php ".src_folder()."/NGS/".$name.".php -in ".data_folder().$name."_in1.bam -out $out_file_bed -loci {$variant_catalog} --log {$out_file_log}");
remove_lines_containing($out_file_vcf, array("##fileDate=", "##reference="));
check_file($out_file_vcf, data_folder().$name."_out1.vcf");
check_file($out_file_bed, data_folder().$name."_out1.bed");
foreach(array("FXN", "ATXN3", "DMPK") as $re)
{
    check_file($out_file_svg_prefix.$re.".fa", data_folder().$name."_out1_".$re.".fa");
    remove_lines_containing($out_file_svg_prefix.$re."_hist.svg", array("<dc:date>")); 
    check_file($out_file_svg_prefix.$re."_hist.svg", data_folder().$name."_out1_".$re."_hist.svg");
       
}

end_test();

?>