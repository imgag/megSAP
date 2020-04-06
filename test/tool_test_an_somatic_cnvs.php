<?php

require_once("framework.php");

$name = "an_somatic_cnvs";

start_test($name);

//1st case: no old CGI data annotated to *_cnvs.tsv file
$cnv_input_1 = data_folder()."/an_somatic_cnvs_in1_cnvs.tsv"; 
$cgi_input_1 = data_folder()."/an_somatic_cnvs_in1_cgi_cnv_analysis.tsv"; 
$cnv_ref_1 = data_folder()."/an_somatic_cnvs_ref1_cnvs.tsv";
$cnv_output_1 = output_folder()."/an_somatic_cnvs_out1_cnvs.tsv";

check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_1 -cnv_in_cgi $cgi_input_1 -include_ncg -out $cnv_output_1 ");
check_file($cnv_output_1,$cnv_ref_1,true);

//2nd case: old CGI data annotated to *_cnvs.tsv file
$cnv_input_2 = data_folder()."/an_somatic_cnvs_in2_cnvs.tsv"; 
$cgi_input_2 = data_folder()."/an_somatic_cnvs_in2_cgi_cnv_analysis.tsv"; 
$cnv_ref_2 = data_folder()."/an_somatic_cnvs_ref2_cnvs.tsv";
$cnv_output_2 = output_folder()."/an_somatic_cnvs_out2_cnvs.tsv";

check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_2 -cnv_in_cgi $cgi_input_2 -out $cnv_output_2");
check_file($cnv_output_2,$cnv_ref_2,true);

//3rd case: input cgi file was filtered for target region -> count of genes is lower in the CGI analysis file than in the input CNVs file
$cnv_input_3 = data_folder()."/an_somatic_cnvs_in3_cnvs.tsv"; 
$cgi_input_3 = data_folder()."/an_somatic_cnvs_in3_cgi_cnv_analysis.tsv";
$cnv_ref_3 = data_folder()."/an_somatic_cnvs_ref3_cnvs.tsv";
$cnv_output_3 = output_folder()."/an_somatic_cnvs_out3_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_3 -cnv_in_cgi $cgi_input_3 -out $cnv_output_3");
check_file($cnv_output_3,$cnv_ref_3,true);

//4th case: annotate RNA data
$rna_counts = data_folder()."/an_somatic_cnvs_rna_counts.tsv";
$cnv_output_4 = output_folder()."/an_somatic_cnvs_out4_cnvs.tsv";
$cnv_ref_4 = data_folder()."/an_somatic_cnvs_ref4_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_3 -rna_counts $rna_counts -out $cnv_output_4 -rna_id RX01_01 -rna_ref_tissue colon");
check_file($cnv_output_4,$cnv_ref_4,false);

//5th case: annotate cytobands
$cnv_output_5 = output_folder()."/an_somatic_cnvs_out5_cnvs.tsv";
$cnv_ref_5 = data_folder()."/an_somatic_cnvs_ref5_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_3 -out $cnv_output_5 -include_cytoband");
check_file($cnv_output_5,$cnv_ref_5,false);
end_test();
?>