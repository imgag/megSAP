<?php

require_once("framework.php");

$name = "an_somatic_gsvar";

start_test($name);

//1.test: annotate RNA BAM file no target
$gsvar_input1 = data_folder() . "an_somatic_gsvar_in1.GSvar";
$gsvar_ref1 = data_folder() . "an_somatic_gsvar_ref1.GSvar";
$rna_bam = data_folder() . "an_somatic_gsvar_in1.bam";
$gsvar_output1 = output_folder() . "an_somatic_gsvar_out1.GSvar"; 
$rna_counts = data_folder() . "an_somatic_gsvar_in1_counts.tsv"; 
check_exec("php ".src_folder()."/Tools/an_somatic_gsvar.php -gsvar_in $gsvar_input1 -rna_bam $rna_bam -rna_counts $rna_counts -rna_id an_somatic_gsvar_in1 -rna_ref_tissue colon -out $gsvar_output1");
check_file($gsvar_output1,$gsvar_ref1,false);

//2. test: NCG Oncogene/TSG for annotation
$gsvar_input2 = data_folder() . "an_somatic_gsvar_in2.GSvar";
$gsvar_ref2 = data_folder() . "an_somatic_gsvar_ref2.GSvar";
$gsvar_output2 = output_folder() . "an_somatic_gsvar_out2.GSvar"; 
check_exec("php ".src_folder()."/Tools/an_somatic_gsvar.php -gsvar_in $gsvar_input2 -out $gsvar_output2 -include_ncg");
check_file($gsvar_output2,$gsvar_ref2,true);

//3.test: annotate RNA BAM file with target
$gsvar_ref3 = data_folder() . "an_somatic_gsvar_ref3.GSvar";
$gsvar_output3 = output_folder() . "an_somatic_gsvar_out3.GSvar"; 
$rna_target = data_folder() . "an_somatic_gsvar_in3.bed"; 
check_exec("php ".src_folder()."/Tools/an_somatic_gsvar.php -gsvar_in $gsvar_input1 -rna_bam $rna_bam -rna_counts $rna_counts -rna_id an_somatic_gsvar_in1 -rna_ref_tissue colon -out $gsvar_output3 -rna_target {$rna_target}" );
check_file($gsvar_output3,$gsvar_ref3,false);

end_test();
?>
