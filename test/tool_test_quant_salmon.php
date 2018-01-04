<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

//calculate correlation of two columns in two tsv files
function column_correlation($tsv1, $tsv2, $col1, $col2)
{
	$table1 = Matrix::fromTSV($tsv1);
	$table1->removeRow(0);
	$table2 = Matrix::fromTSV($tsv2);
	$table2->removeRow(0);
	$corr = correlation($table1->getCol($col1), $table2->getCol($col2));
	return($corr);
}

$name = "quant_salmon";
start_test($name);

//input files
$in1_file = data_folder()."mapping_star_in1.fastq.gz";
$in2_file = data_folder()."mapping_star_in2.fastq.gz";

//single-end
$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 {$in1_file} -out {$out_file1}");
//check_file($out_file1, data_folder().$name."_out1.tsv");
//check TPM correlation, at least 0.995
check(column_correlation($out_file1, data_folder().$name."_out1.tsv", 3, 3), 1.0, 0.005);

//paired-end
$in_file2 = data_folder().$name."_in2.tsv";
$out_file2 = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 {$in1_file} -in2 {$in2_file} -out {$out_file2}");
//check_file($out_file2, data_folder().$name."_out2.tsv");
//check TPM correlation, at least 0.995
check(column_correlation($out_file2, data_folder().$name."_out2.tsv", 3, 3), 1.0, 0.005);



end_test();