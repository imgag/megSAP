<?php

require_once("framework.php");

$name = "runqc_parser_ont";
$file = data_folder().$name;

start_test($name);

//init NGSD
init_ngsd($name);

$db_con = DB::getInstance("NGSD_TEST");

//extract run data
exec2("rm -rf ".output_folder()."/21073LRa277_01234");
exec2("tar -xvzf ".data_folder().$name."_in1.tar.gz -C ".output_folder());

//test 1: no db import
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -no_db");
check($output[0], "#run\tfc_id\tread_num\tyield\tpassing_filter_perc\tfraction_skipped\tq30_perc\tq20_perc\tN50\tprotocol_id\tsoftware_args\tdevice_firmware_versions\tminknow_version\tsub_folder");
check($output[1], "01234\tPAW01234\t9574413\t89394689391\t96.423576795467\t3.13335136E-7\t0.10337179157905\t47.327042850332\t14183\tsequencing/sequencing_PRO114_DNA_e8_2_400K:FLO-PRO114M:SQK-LSK114-XL:400\t"
					."--fast5=off --pod5=on --fastq=off --bam=on --generate_bulk_file=off --active_channel_selection=on --base_calling=on --bam_batch_duration=3600 --mux_scan_period=1.5 --pore_reserve=on "
					."--guppy_filename=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg --read_filtering min_qscore=7 --min_read_length=200 --split_files_by_barcode=off\t"
					."Hublett Board:2.1.10;Satellite Board:2.1.9\t5.9.7\t20000101_1200_3D_PAW01234_abcdef12");
check($output[2], "01234\tPAW01234\t9574413\t89394689391\t96.423576795467\t3.13335136E-7\t0.10337179157905\t47.327042850332\t14183\tsequencing/sequencing_PRO114_DNA_e8_2_400K:FLO-PRO114M:SQK-LSK114-XL:400\t"
					."--fast5=off --pod5=on --fastq=off --bam=on --generate_bulk_file=off --active_channel_selection=on --base_calling=on --bam_batch_duration=3600 --mux_scan_period=1.5 --pore_reserve=on "
					."--guppy_filename=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg --read_filtering min_qscore=7 --min_read_length=200 --split_files_by_barcode=off\t"
					."Hublett Board:2.1.10;Satellite Board:2.1.9\t5.9.7\tcombined");
check($db_con->getValue("SELECT COUNT(*) FROM runqc_ont"), 0);

//test 2: with db import
check_exec("php ".src_folder()."/Tools/{$name}.php -name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST");
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


//test 3: re-import without '-force'
list($stdout, $stderr, $exit_code) = exec2("php ".src_folder()."/Tools/{$name}.php -name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST", FALSE); //should fail
check(($exit_code == 0), FALSE); // check exit code
check((strpos(implode("", $stderr), "QC data for this run and read was already imported. Use the flag ") != FALSE), TRUE); //check error message

//test 4: re-import with '-force'
check_exec("php ".src_folder()."/Tools/{$name}.php -name '#01234' -run_dir ".output_folder()."/21073LRa277_01234 -db NGSD_TEST -force");
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

end_test();


?>
