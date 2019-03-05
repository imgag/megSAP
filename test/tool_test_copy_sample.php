<?php

require_once("framework.php");

$name = "copy_sample";
$file = data_folder().$name;

start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
mkdir("Unaligned");



//test 1 - default project handling, projects with 'analyze=fastq'
$out_file = output_folder().$name."_out1Makefile";
$runinfo = data_folder().$name."_in1_RunInfo.xml";
check_exec("php ".src_folder()."/NGS/{$name}.php -samplesheet {$file}_in1.csv -out {$out_file} -db NGSD_TEST -runinfo {$runinfo}");
check_file($out_file, data_folder()."{$name}_out1Makefile");

//test 2 - NextSeq, several lanes per sample, somatic samples
$out_file = output_folder().$name."_out2Makefile";
$runinfo = data_folder().$name."_in2_RunInfo.xml";
check_exec("php ".src_folder()."/NGS/{$name}.php -samplesheet {$file}_in2.csv -out {$out_file} -db NGSD_TEST -runinfo {$runinfo}");
check_file($out_file, data_folder()."{$name}_out2Makefile");

rmdir("Unaligned");
end_test();

?>
