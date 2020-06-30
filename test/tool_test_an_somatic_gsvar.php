<?php

require_once("framework.php");

$name = "an_somatic_gsvar";

start_test($name);

//2. test: input file without known CGI cancer acronym translation
$gsvar_input1 = data_folder() . "an_somatic_gsvar_in1.GSvar";
$gsvar_ref1 = data_folder() . "an_somatic_gsvar_ref1.GSvar";
$cgi_mutation_file_input1 = data_folder() . "an_somatic_gsvar_in1_cgi_mutation_analysis.tsv";
$gsvar_output1 = output_folder() . "an_somatic_gsvar_out1.GSvar"; 
check_exec("php ".src_folder()."/NGS/an_somatic_gsvar.php -gsvar_in $gsvar_input1 -cgi_snv_in $cgi_mutation_file_input1 -out $gsvar_output1");
check_file($gsvar_output1,$gsvar_ref1,true);

//2. test: input file with known CGI cancer acronym translation, include NCG6.0 Oncogene/TSG for annotation
$gsvar_input2 = data_folder() . "an_somatic_gsvar_in2.GSvar";
$gsvar_ref2 = data_folder() . "an_somatic_gsvar_ref2.GSvar";
$cgi_mutation_file_input2 = data_folder() . "an_somatic_gsvar_in2_cgi_mutation_analysis.tsv";
$gsvar_output2 = output_folder() . "an_somatic_gsvar_out2.GSvar"; 
check_exec("php ".src_folder()."/NGS/an_somatic_gsvar.php -gsvar_in $gsvar_input2 -cgi_snv_in $cgi_mutation_file_input2 -out $gsvar_output2 -include_ncg");
check_file($gsvar_output2,$gsvar_ref2,true);

//3.test: annotate RNA BAM file
$gsvar_input3 = data_folder() . "an_somatic_gsvar_in3.GSvar";
$gsvar_ref3 = data_folder() . "an_somatic_gsvar_ref3.GSvar";
$rna_bam = data_folder() . "an_somatic_gsvar_in3.bam";
$gsvar_output3 = output_folder() . "an_somatic_gsvar_out3.GSvar"; 
$rna_counts = data_folder() . "an_somatic_gsvar_in3_counts.tsv"; 
check_exec("php ".src_folder()."/NGS/an_somatic_gsvar.php -gsvar_in $gsvar_input3 -rna_bam $rna_bam -rna_counts $rna_counts -rna_id an_somatic_gsvar_in3 -rna_ref_tissue colon -out $gsvar_output3");
check_file($gsvar_output3,$gsvar_ref3,false);
end_test();
?>
