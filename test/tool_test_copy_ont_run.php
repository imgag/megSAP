<?php

require_once("framework.php");

$name = "copy_ont_run";
$file = data_folder().$name;

start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//extract run data
exec2("rm -rf ".output_folder()."/21073LRa277_01234");
exec2("tar -xvzf ".data_folder().$name."_in1.tar.gz -C ".output_folder());

//test 1: bam
check_exec("php ".src_folder()."/NGS/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -bam --log -");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam", data_folder()."{$name}_out1.bam");
exec2("rm -rf ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01");
//check runQC import
$db_con = DB::getInstance("NGSD_TEST");
check($db_con->getValue("SELECT COUNT(*) FROM runqc_ont"), 1);
check($db_con->getValue("SELECT read_num FROM runqc_ont WHERE sequencing_run_id=1"), 9574413);
check($db_con->getValue("SELECT yield FROM runqc_ont WHERE sequencing_run_id=1"), 89394689391);
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

//test 2: fastq
check_exec("php ".src_folder()."/NGS/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -fastq --log -");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
exec2("rm -rf ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01");

//test 3: fastq single file
check_exec("php ".src_folder()."/NGS/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -fastq -single_fastq ".output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz --log -");
check_file(output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
exec2("rm ".output_folder()."/21073LRa277_01234/TEST_21073LRa277.fastq.gz");

//test 4: fastq+bam and queueing
check_exec("php ".src_folder()."/NGS/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -fastq -bam -queue_sample --log -");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.fastq.gz", data_folder()."{$name}_out2.fastq.gz");
check_file(output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01/21073LRa277_01.mod.unmapped.bam", data_folder()."{$name}_out1.bam");
exec2("rm -rf ".output_folder()."/21073LRa277_01234/TEST_Sample_21073LRa277_01");


check($db_con->getValue("SELECT status FROM sequencing_run WHERE name = '#01234'"), "analysis_started");
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job"), 1);
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job_history"), 1);
check($db_con->getValue("SELECT COUNT(*) FROM analysis_job_sample"), 1);

//test 5: skipped basecalling -> fail
exec2("dd if=/dev/zero of=".output_folder()."/21073LRa277_01234/20000101_1200_3D_PAW01234_abcdef12/pod5_skip/trash.pod5  bs=2M  count=1"); //create data in pod5_skip folder
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/NGS/{$name}.php -run_name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -threads 4 -bam --log -", FALSE); //should fail
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "'pod5_skip' directory present in ") != FALSE), TRUE); //check error message
check((strpos(implode("", $stderr), "some data has not been basecalled!") != FALSE), TRUE);

end_test();


?>
