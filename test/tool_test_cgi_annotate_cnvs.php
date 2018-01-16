<?php

require_once("framework.php");

$name = "cgi_annotate_cnvs";

start_test($name);

//1st case: no old CGI data annotated to *_cnvs.tsv file
$cnv_input_1 = data_folder()."/cgi_annotate_cnvs_in1_cnvs.tsv"; 
$cgi_input_1 = data_folder()."/cgi_annotate_cnvs_in1_cgi_cnv_analysis.tsv"; 
$cnv_ref_1 = data_folder()."/cgi_annotate_cnvs_ref1_cnvs.tsv";
$cnv_output_1 = output_folder()."/cgi_annotate_cnvs_out1_cnvs.tsv";

check_exec("php ".src_folder()."/NGS/cgi_annotate_cnvs.php -cnv_in $cnv_input_1 -cnv_in_cgi $cgi_input_1 -out $cnv_output_1");
check_file($cnv_output_1,$cnv_ref_1,true);

//2nd case: old CGI data annotated to *_cnvs.tsv file
$cnv_input_2 = data_folder()."/cgi_annotate_cnvs_in2_cnvs.tsv"; 
$cgi_input_2 = data_folder()."/cgi_annotate_cnvs_in2_cgi_cnv_analysis.tsv"; 
$cnv_ref_2 = data_folder()."/cgi_annotate_cnvs_ref2_cnvs.tsv";
$cnv_output_2 = output_folder()."/cgi_annotate_cnvs_out2_cnvs.tsv";

check_exec("php ".src_folder()."/NGS/cgi_annotate_cnvs.php -cnv_in $cnv_input_2 -cnv_in_cgi $cgi_input_2 -out $cnv_output_2");
check_file($cnv_output_2,$cnv_ref_2,true);

//3rd case: input cgi file was filtered for target region -> count of genes is lower in the CGI analysis file than in the input CNVs file
$cnv_input_3 = data_folder()."/cgi_annotate_cnvs_in3_cnvs.tsv"; 
$cgi_input_3 = data_folder()."/cgi_annotate_cnvs_in3_cgi_cnv_analysis.tsv";
$cnv_ref_3 = data_folder()."/cgi_annotate_cnvs_ref3_cnvs.tsv";
$cnv_output_3 = output_folder()."/cgi_annotate_cnvs_out3_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/cgi_annotate_cnvs.php -cnv_in $cnv_input_3 -cnv_in_cgi $cgi_input_3 -out $cnv_output_3");
check_file($cnv_output_3,$cnv_ref_3,true);

?>