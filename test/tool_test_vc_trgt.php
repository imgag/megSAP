<?php

require_once("framework.php");

$name = "vc_trgt";
start_test($name);

//test
$out_file_vcf = output_folder().$name."_out1.vcf";
$out_file_log = output_folder().$name."_out1.log";
$out_file_svg_prefix = output_folder()."repeat_expansions/out1_HG002_repeats_";
$out_file_bam = output_folder()."repeat_expansions/out1_HG002_spanning.bam";
$variant_catalog = data_folder().$name."_variant_catalog_in1.bed";
check_exec("php ".src_folder()."/Tools/".$name.".php -in ".data_folder().$name."_in1.cram -threads 4 -sample_name out1_HG002 -gender male -out {$out_file_vcf} -loci {$variant_catalog} --log {$out_file_log}");
remove_lines_containing($out_file_vcf, array("##fileDate=", "##reference=", "##trgtCommand"));
check_file($out_file_vcf, data_folder().$name."_out1.vcf");
foreach(array("FXN", "ATXN3", "DMPK") as $re)
{
    remove_lines_containing($out_file_svg_prefix.$re.".svg", array("<dc:date>")); 
    check_file($out_file_svg_prefix.$re.".svg", data_folder().$name."_out1_HG002_repeats_".$re.".svg");
}
check(file_exists($out_file_bam), true);
check(filesize($out_file_bam), 140000, 10000);

end_test();

?>