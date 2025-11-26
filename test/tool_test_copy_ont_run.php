<?php

require_once("framework.php");

$name = "copy_ont_run";
$file = data_folder().$name;

start_test($name);

//init NGSD
init_ngsd($name);

//extract run data
exec2("rm -rf ".output_folder()."/21073LRa277_01234");
exec2("tar -xvzf ".data_folder().$name."_in1.tar.gz -C ".output_folder());

//test 1: bam
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 --log ".output_folder()."{$name}_1.log");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam", data_folder()."{$name}_out1.bam");
# move folder to don't interfere with next tests
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test1");

//check runQC import
$db_con = DB::getInstance("NGSD_TEST");
check($db_con->getValue("SELECT COUNT(*) FROM runqc_ont"), 1);
check($db_con->getValue("SELECT read_num FROM runqc_ont WHERE sequencing_run_id=1"), 9574413);
check($db_con->getValue("SELECT yield FROM runqc_ont WHERE sequencing_run_id=1"), 89394700000);
check($db_con->getValue("SELECT passing_filter_perc FROM runqc_ont WHERE sequencing_run_id=1"), 96.4236);
check($db_con->getValue("SELECT fraction_skipped FROM runqc_ont WHERE sequencing_run_id=1"), 0.000000313335);
check($db_con->getValue("SELECT q30_perc FROM runqc_ont WHERE sequencing_run_id=1"), 0.103372);
check($db_con->getValue("SELECT q20_perc FROM runqc_ont WHERE sequencing_run_id=1"), 47.327);
check($db_con->getValue("SELECT n50 FROM runqc_ont WHERE sequencing_run_id=1"), 14183);
check($db_con->getValue("SELECT protocol_id FROM runqc_ont WHERE sequencing_run_id=1"), "sequencing/sequencing_PRO114_DNA_e8_2_400K:FLO-PRO114M:SQK-LSK114-XL:400");
check($db_con->getValue("SELECT software_args FROM runqc_ont WHERE sequencing_run_id=1"), "--fast5=off --pod5=on --fastq=off --bam=on --generate_bulk_file=off --active_channel_selection=on --base_calling=on "
																							."--bam_batch_duration=3600 --mux_scan_period=1.5 --pore_reserve=on --guppy_filename=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg "
																							."--read_filtering min_qscore=7 --min_read_length=200 --split_files_by_barcode=off");
check($db_con->getValue("SELECT device_firmware_versions FROM runqc_ont WHERE sequencing_run_id=1"), "Hublett Board:2.1.10;Satellite Board:2.1.9");
check($db_con->getValue("SELECT minknow_version FROM runqc_ont WHERE sequencing_run_id=1"), "5.9.7");

//TODO: recheck if fastq support is back
/*
//test 2: fastq
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -fastq --log ".output_folder()."{$name}_2.log");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
# move folder to don't interfere with next tests
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test2");

//test 3: fastq single file
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -fastq -single_fastq ".output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz --log ".output_folder()."{$name}_3.log");
check_file(output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
exec2("rm ".output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz");
*/

//test 4: fastq+bam and queueing
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -queue_sample --log ".output_folder()."{$name}_4.log");
// check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam", data_folder()."{$name}_out1.bam");
# move folder to don't interfere with next tests
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test4");


check($db_con->getValue("SELECT status FROM sequencing_run WHERE name = '#01234'"), "analysis_started");
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job"), 1);
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job_history"), 1);
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job_sample"), 1);

//test 5: skipped basecalling -> fail
exec2("dd if=/dev/zero of=".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/pod5_skip/trash.pod5  bs=8M  count=1"); //create data in pod5_skip folder
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/IMGAG/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 --log ".output_folder()."{$name}_5.log", FALSE); //should fail
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "'pod5_skip' directory present in ") != FALSE), TRUE); //check error message
check((strpos(implode("", $stderr), "some data has not been basecalled!") != FALSE), TRUE);

//test 6: queue skipped
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -queue_basecalling --log ".output_folder()."{$name}_6.log");
//check log file
$logs = file(output_folder()."{$name}_6.log", FILE_IGNORE_NEW_LINES);
foreach ($logs as $line) 
{
	if (str_contains($line, "NOTICE: 'SGE command:"))
	{
		check(str_contains($line, "src/Tools/basecall_ont_run.php"), TRUE);
		check(str_contains($line, "tool_test_copy_ont_run/21073LRa277_01234"), TRUE);
		check(str_contains($line, "Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam"), TRUE);
		check(str_contains($line, "-basecall_model hac -min_qscore 9 -secondary_output"), TRUE);
		check(str_contains($line, "-skipped_only -rename_skip_folder"), TRUE);
		break;
	}
}
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test6");

//test 7: queue full basecalling (sup)
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model sup -queue_basecalling --log ".output_folder()."{$name}_7.log");
//check log file
$logs = file(output_folder()."{$name}_7.log", FILE_IGNORE_NEW_LINES);
foreach ($logs as $line) 
{
	if (str_contains($line, "NOTICE: 'SGE command:"))
	{
		check(str_contains($line, "src/Tools/basecall_ont_run.php"), TRUE);
		check(str_contains($line, "tool_test_copy_ont_run/21073LRa277_01234"), TRUE);
		check(str_contains($line, "Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam"), TRUE);
		check(str_contains($line, "-basecall_model sup -min_qscore 9 -secondary_output"), TRUE);
		check(str_contains($line, "-skipped_only"), FALSE);
		break;
	}
}
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test7");

// remove pod5s from pod_skip and rename bam folder back
exec2("rm ".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/pod5_skip/trash.pod5");
exec2("mv ".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/bam_pass_hac ".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/bam_pass");

//test 8: test wrong basecall model -> fail
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model sup --log ".output_folder()."{$name}_8.log", false);
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "ERROR: 'Basecall Models doesn't match: On-device basecall model: 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2', requested basecall model: 'sup'. Please check setting or re-do basecalling!'") != FALSE), TRUE); //check error message

// add second flowcell with different basecall model
exec2("mkdir ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/");
exec2("mkdir ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/bam_pass");
exec2("mkdir ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/bam_fail");
exec2("mkdir ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/pod5");
exec2("cp ".data_folder().$name."_in2.bam ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/bam_pass/");
exec2("cp ".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/report_PAW01234_20000101_0000_abcdef12.json ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/report_PAW05678_20000002_0000_abcdef12.json");

//test 9: test mixed flowcell id -> fail
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model sup --log ".output_folder()."{$name}_9.log", false);
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "Flowcell ID 'PAW01234' not found in directory name") != FALSE), TRUE); //check error message

//test 10: test mixed basecall models -> fail
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model sup -ignore_flowcell_id_check --log ".output_folder()."{$name}_10.log", false);
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "Multiple basecall models found (dna_r10.4.1_e8.2_400bps_sup@v4.1.0, dna_r10.4.1_e8.2_400bps_hac@v3.5.2)!") != FALSE), TRUE); //check error message

//test 11: test mixed basecall and uncalled flow cell runs -> fail
exec2("rm ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/bam_pass/{$name}_in2.bam");
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model hac -ignore_flowcell_id_check --log ".output_folder()."{$name}_11.log", false);
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "Multiple basecall models found (undefined, dna_r10.4.1_e8.2_400bps_hac@v3.5.2)!") != FALSE), TRUE); //check error message

//copy bams
exec2("cp ".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/bam_pass/*.bam ".output_folder()."/21073LRa277_01234/20000002_1200_3D_PAW05678_abcdef12/bam_pass/");

//test 12: import QC from 2 runs
check_exec("php ".src_folder()."/IMGAG/{$name}.php -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -basecall_model hac -ignore_flowcell_id_check --log ".output_folder()."{$name}_12.log");
# move folder to don't interfere with next tests
exec2("mv ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01 ".output_folder()."/TEST_Sample_21073LRa277_01_test12");

//check runQC import
check($db_con->getValue("SELECT COUNT(*) FROM runqc_ont"), 1);
check($db_con->getValue("SELECT read_num FROM runqc_ont WHERE sequencing_run_id=1"), 19148826);
check($db_con->getValue("SELECT yield FROM runqc_ont WHERE sequencing_run_id=1"), 178789000000);
check($db_con->getValue("SELECT passing_filter_perc FROM runqc_ont WHERE sequencing_run_id=1"), 96.4236);
check($db_con->getValue("SELECT fraction_skipped FROM runqc_ont WHERE sequencing_run_id=1"), 0.000000313335);
check($db_con->getValue("SELECT q30_perc FROM runqc_ont WHERE sequencing_run_id=1"), 0.103372);
check($db_con->getValue("SELECT q20_perc FROM runqc_ont WHERE sequencing_run_id=1"), 47.327);
check($db_con->getValue("SELECT n50 FROM runqc_ont WHERE sequencing_run_id=1"), 14183);
check($db_con->getValue("SELECT protocol_id FROM runqc_ont WHERE sequencing_run_id=1"), "sequencing/sequencing_PRO114_DNA_e8_2_400K:FLO-PRO114M:SQK-LSK114-XL:400");
check($db_con->getValue("SELECT software_args FROM runqc_ont WHERE sequencing_run_id=1"), "--fast5=off --pod5=on --fastq=off --bam=on --generate_bulk_file=off --active_channel_selection=on --base_calling=on "
																							."--bam_batch_duration=3600 --mux_scan_period=1.5 --pore_reserve=on --guppy_filename=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg "
																							."--read_filtering min_qscore=7 --min_read_length=200 --split_files_by_barcode=off");
check($db_con->getValue("SELECT device_firmware_versions FROM runqc_ont WHERE sequencing_run_id=1"), "Hublett Board:2.1.10;Satellite Board:2.1.9");
check($db_con->getValue("SELECT minknow_version FROM runqc_ont WHERE sequencing_run_id=1"), "5.9.7");



end_test();


?>
