<?php

require_once("framework.php");

$name = "extract_signatures";
start_test($name);

$seeds = data_folder()."/extract_signatures_seeds.txt";

//snvs
$out = output_folder()."{$name}_out0/";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder()."/extract_signatures_in0.vcf -mode snv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out0.log");

check_file_tsv($out."De_Novo_map_to_COSMIC_SBS96.tsv", data_folder().$name."_out0_SBS96.tsv", false, "", "");
check_file_tsv($out."De_Novo_map_to_COSMIC_ID83.tsv", data_folder().$name."_out0_ID83.tsv", false, "", "");
check_file_tsv($out."De_Novo_map_to_COSMIC_DBS78.tsv", data_folder().$name."_out0_DBS78.tsv", false, "", "");

//snvs - empty file
$out = output_folder()."{$name}_out1/";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder()."/extract_signatures_in3.vcf -mode snv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out1.log");

//cnvs
$out = output_folder()."{$name}_out2/";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder()."/extract_signatures_clincnv_in1.tsv -mode cnv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out2.log");

check_file_tsv($out."De_Novo_map_to_COSMIC_CNV48.tsv", data_folder().$name."_clincnv_out1_CNV48.tsv", false, "", "");

//cnvs empty file (no cnvs):
$out = output_folder()."{$name}_out3/";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder()."/extract_signatures_clincnv_in2.tsv -mode cnv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out3.log");

end_test();

?>