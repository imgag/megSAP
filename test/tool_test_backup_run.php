<?php

require_once("framework.php");

$name = "backup_run";
$file = data_folder().$name;
start_test($name);

//init
// check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
exec("rm -rf ".output_folder()."20230823_LH00240_0003_AFCID0001_00123");
exec("rm -rf ".output_folder()."backup");
exec("rm -rf ".output_folder()."restore");
exec("unzip ".data_folder().$name."_in1.zip -d ".output_folder());
mkdir(output_folder()."backup");
mkdir(output_folder()."restore");

//test 1 - try to backup NovaSeqX run with 2 analysis folder (should fail)
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".output_folder()."20230823_LH00240_0003_AFCID0001_00123 -out_folder ".output_folder()."backup -test --log ".output_folder()."test1.log", FALSE);

//test 2 - remove 2nd analysis folder  backup NovaSeqX run 
exec("rm -d ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/2");
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".output_folder()."20230823_LH00240_0003_AFCID0001_00123 -out_folder ".output_folder()."backup -test --log ".output_folder()."test2.log", TRUE);
check_exec("tar xzf ".output_folder()."backup/230823_LH00240_0003_AFCID0001_00123.tar.gz -C ".output_folder()."restore/");
//check files
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/BCLConvert/fastq/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/BCLConvert/fastq/";
check_file("{$ref_folder}RX123456_03_S5_L006_R1_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L006_R1_001.fastq.gz");
check_file("{$ref_folder}RX123456_03_S5_L006_R2_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L006_R2_001.fastq.gz");
check_file("{$ref_folder}RX123456_03_S5_L006_R3_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L006_R3_001.fastq.gz");
check_file("{$ref_folder}RX123456_03_S5_L008_R1_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L008_R1_001.fastq.gz");
check_file("{$ref_folder}RX123456_03_S5_L008_R2_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L008_R2_001.fastq.gz");
check_file("{$ref_folder}RX123456_03_S5_L008_R3_001.fastq.gz", "{$restored_folder}RX123456_03_S5_L008_R3_001.fastq.gz");
check_file("{$ref_folder}Undetermined_S5_L006_R1_001.fastq.gz", "{$restored_folder}Undetermined_S5_L006_R1_001.fastq.gz");
check_file("{$ref_folder}Undetermined_S5_L006_R2_001.fastq.gz", "{$restored_folder}Undetermined_S5_L006_R2_001.fastq.gz");
check_file("{$ref_folder}Undetermined_S5_L006_R3_001.fastq.gz", "{$restored_folder}Undetermined_S5_L006_R3_001.fastq.gz");
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/ora_fastq/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/ora_fastq/";
check_file("{$ref_folder}DX181277_05_S5_L005_R1_001.fastq.ora", "{$restored_folder}DX181277_05_S5_L005_R1_001.fastq.ora");
check_file("{$ref_folder}DX181277_05_S5_L005_R2_001.fastq.ora", "{$restored_folder}DX181277_05_S5_L005_R2_001.fastq.ora");
check_file("{$ref_folder}DX181278_05_S5_L007_R1_001.fastq.ora", "{$restored_folder}DX181278_05_S5_L007_R1_001.fastq.ora");
check_file("{$ref_folder}DX181278_05_S5_L007_R2_001.fastq.ora", "{$restored_folder}DX181278_05_S5_L007_R2_001.fastq.ora");
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/fastq/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenEnrichment/fastq/";
check_file("{$ref_folder}DX181277_05_S5_L005_R1_001.fastq.gz", "{$restored_folder}DX181277_05_S5_L005_R1_001.fastq.gz");
check_file("{$ref_folder}DX181277_05_S5_L005_R2_001.fastq.gz", "{$restored_folder}DX181277_05_S5_L005_R2_001.fastq.gz");
check_file("{$ref_folder}DX181278_05_S5_L007_R1_001.fastq.gz", "{$restored_folder}DX181278_05_S5_L007_R1_001.fastq.gz");
check_file("{$ref_folder}DX181278_05_S5_L007_R2_001.fastq.gz", "{$restored_folder}DX181278_05_S5_L007_R2_001.fastq.gz");
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenGermline/DX180049_05/germline_seq/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Analysis/3/Data/DragenGermline/DX180049_05/germline_seq/";
check(file_exists("{$restored_folder}DX180049_05.hard-filtered.gvcf.gz"), false);
check(file_exists("{$restored_folder}DX180049_05.hard-filtered.vcf.gz"), false);
check(file_exists("{$restored_folder}DX180049_05.sv.vcf.gz"), false);
check(file_exists("{$restored_folder}DX180049_05.bam"), false);
check(file_exists("{$restored_folder}DX180049_05.bam.bai"), false);
check(file_exists("{$restored_folder}DX180049_05.cram"), false);
check(file_exists("{$restored_folder}DX180049_05.cram.crai"), false);
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Data/Intensities/BaseCalls/L001/";
check(file_exists("{$restored_folder}s_1_1101.filter"), false);
check(file_exists("{$restored_folder}s_1_1102.filter"), false);
check(file_exists("{$restored_folder}C1.1/L001_1.cbcl"), false);
check(file_exists("{$restored_folder}C1.1/L001_2.cbcl"), false);
check(file_exists("{$restored_folder}C2.1/L001_1.cbcl"), false);
check(file_exists("{$restored_folder}C2.1/L001_2.cbcl"), false);

//test 3 - normal NovaSeq 6000 run
exec("rm ".output_folder()."backup/230823_LH00240_0003_AFCID0001_00123.tar.gz");
exec("rm -r ".output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123");
exec("cp ".data_folder().$name."_in_runparameters1.xml ".output_folder()."20230823_LH00240_0003_AFCID0001_00123/RunParameters.xml"); //replace RunParameters.xml
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".output_folder()."20230823_LH00240_0003_AFCID0001_00123 -out_folder ".output_folder()."backup -test --log ".output_folder()."test2.log", TRUE);
check_exec("tar xzf ".output_folder()."backup/230823_LH00240_0003_AFCID0001_00123.tar.gz -C ".output_folder()."restore/");

//check files
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Data/Intensities/BaseCalls/L001/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Data/Intensities/BaseCalls/L001/";
check_file("{$restored_folder}s_1_1101.filter", "{$ref_folder}s_1_1101.filter");
check_file("{$restored_folder}s_1_1102.filter", "{$ref_folder}s_1_1102.filter");
check_file("{$restored_folder}C1.1/L001_1.cbcl", "{$ref_folder}C1.1/L001_1.cbcl");
check_file("{$restored_folder}C1.1/L001_2.cbcl", "{$ref_folder}C1.1/L001_2.cbcl");
check_file("{$restored_folder}C2.1/L001_1.cbcl", "{$ref_folder}C2.1/L001_1.cbcl");
check_file("{$restored_folder}C2.1/L001_2.cbcl", "{$ref_folder}C2.1/L001_2.cbcl");
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Images/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Images/";
check(file_exists("{$restored_folder}pic1.png"), false);
check(file_exists("{$restored_folder}pic2.png"), false);
check(file_exists("{$restored_folder}pic3.png"), false);
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Thumbnail_Images/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Thumbnail_Images/";
check(file_exists("{$restored_folder}pic1.png"), false);
check(file_exists("{$restored_folder}pic2.png"), false);
check(file_exists("{$restored_folder}pic3.png"), false);
$ref_folder = output_folder()."20230823_LH00240_0003_AFCID0001_00123/Unaligned/";
$restored_folder = output_folder()."restore/20230823_LH00240_0003_AFCID0001_00123/Unaligned/";
check(file_exists("{$restored_folder}Sample_DX181277_05/DX181277_05_S5_L005_R1_001.fastq.gz"), false);
check(file_exists("{$restored_folder}Sample_DX181277_05/DX181277_05_S5_L005_R2_001.fastq,gz"), false);
check(file_exists("{$restored_folder}Sample_DX181278_05/DX181278_05_S5_L007_R1_001.fastq.gz"), false);
check(file_exists("{$restored_folder}Sample_DX181278_05/DX181278_05_S5_L007_R2_001.fastq.gz"), false);

end_test();

?>
