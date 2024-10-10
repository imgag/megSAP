<?php

require_once("framework.php");


function check_signature_file($out, $ref)
{
	check_test_started();
	
	$passed = true;
	$logfile = $out."_diff";
	
	exec("numdiff -a 0.01 -D '% \  , \n \t ( )' $ref $out > $logfile 2>&1", $output, $return);
	$passed = ($return==0);
	
	if ($passed)
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED (see $logfile)";
		++$GLOBALS["failed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}

$name = "extract_signatures";
start_test($name);

$seeds = data_folder()."/extract_signatures_seeds.txt";

//snvs
$out = output_folder()."{$name}_out0/";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder()."/extract_signatures_in0.vcf -mode snv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out1.log");

$out_file_SBS = $out."De_Novo_map_to_COSMIC_SBS96.csv";
check_signature_file($out_file_SBS, data_folder().$name."_out0_SBS96.csv", true);
$out_file_ID = $out."De_Novo_map_to_COSMIC_ID83.csv";
check_signature_file($out_file_ID, data_folder().$name."_out0_ID83.csv", true);
$out_file_DBS = $out."De_Novo_map_to_COSMIC_DBS78.csv";
check_signature_file($out_file_DBS, data_folder().$name."_out0_DBS78.csv", true);

//snvs - empty file
$out = output_folder()."{$name}_out3/";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder()."/extract_signatures_in3.vcf -mode snv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out2.log");

//cnvs
$out = output_folder()."{$name}_out1/";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder()."/extract_signatures_clincnv_in1.tsv -mode cnv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out3.log");

$out_file_CNV = $out."De_Novo_map_to_COSMIC_CNV48.csv";
check_signature_file($out_file_CNV, data_folder().$name."_clincnv_out1_CNV48.csv", true);

//cnvs empty file (no cnvs):
$out = output_folder()."{$name}_out2/";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder()."/extract_signatures_clincnv_in2.tsv -mode cnv -reference GRCh38 -threads 4 -out $out -seeds $seeds -nmfRep 10 --log ".output_folder().$name."_out4.log");

end_test();

?>