<?php

require_once("framework.php");

$name = "copy_sample";
$file = data_folder().$name;

start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
exec("rm -rf 190228_NB501582_0169_AH5LG5BDXX_0*");
mkdir(output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001");
mkdir(output_folder()."190228_NB501582_0169_AH5LG5BDXX_01489");
mkdir(output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/Unaligned");
mkdir(output_folder()."190228_NB501582_0169_AH5LG5BDXX_01489/Unaligned");
exec("cp ".data_folder().$name."_in1.csv ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/");
exec("cp ".data_folder().$name."_in2.csv ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/");
exec("cp ".data_folder().$name."_in3.csv ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_01489/");
exec("cp ".data_folder().$name."_in_runparameters1.xml ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/RunParameters.xml");
exec("cp ".data_folder().$name."_in_runparameters1.xml ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_01489/RunParameters.xml");
exec("unzip ".data_folder().$name."_in4.zip -d ".output_folder());
exec("cp ".data_folder().$name."_in_runparameters2.xml ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/RunParameters.xml");

//test 1 - default project handling, projects with 'analyze=fastq'
$out_file = output_folder().$name."_out1Makefile";
$runinfo = data_folder().$name."_in1_RunInfo.xml";
check_exec("cd ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/ && php ".src_folder()."/NGS/{$name}.php -test -samplesheet {$name}_in1.csv -out {$out_file} -db NGSD_TEST -runinfo {$runinfo}");
check_file($out_file, data_folder()."{$name}_out1Makefile");

//test 2 - NextSeq, several lanes per sample, somatic samples
$out_file = output_folder().$name."_out2Makefile";
$runinfo = data_folder().$name."_in2_RunInfo.xml";
check_exec("cd ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_00001/ && php ".src_folder()."/NGS/{$name}.php -test -samplesheet {$name}_in2.csv -out {$out_file} -db NGSD_TEST -runinfo {$runinfo}");
check_file($out_file, data_folder()."{$name}_out2Makefile");

//test 3 - Trio handling from GenLab
$out_file = output_folder().$name."_out3Makefile";
$runinfo = data_folder().$name."_in3_RunInfo.xml";
check_exec("cd ".output_folder()."190228_NB501582_0169_AH5LG5BDXX_01489/ && php ".src_folder()."/NGS/{$name}.php -test -samplesheet {$name}_in3.csv -out {$out_file} -db NGSD_TEST -runinfo {$runinfo}");
check_file($out_file, data_folder()."{$name}_out3Makefile");

//test 4 - NovaSeq X with analysis
$out_file = output_folder().$name."_out4Makefile";
check_exec("cd ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/ && php ".src_folder()."/NGS/{$name}.php -test -out {$out_file} -db NGSD_TEST");
check_file($out_file, data_folder()."{$name}_out4Makefile");

//test 5 - NovaSeq X with analysis (CRAM)
$out_file = output_folder().$name."_out5Makefile";
// remove BAMs
exec("rm ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/DX181277_05/enrichment_seq/DX181277_05.bam");
exec("rm ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/DX181277_05/enrichment_seq/DX181277_05.bam.bai");
exec("rm ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/DX181278_05/enrichment_seq/DX181278_05.bam");
exec("rm ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/DX181278_05/enrichment_seq/DX181278_05.bam.bai");
check_exec("cd ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/ && php ".src_folder()."/NGS/{$name}.php -test -out {$out_file} -db NGSD_TEST");
check_file($out_file, data_folder()."{$name}_out5Makefile");

end_test();

?>
