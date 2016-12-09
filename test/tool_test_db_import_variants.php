<?php

require_once("framework.php");

$name = "db_import_variants";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//germline - import 3 samples
$log1_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120676_01 -var ".data_db_folder()."GS120676_01_annotated.tsv -db NGSD_TEST --log $log1_file");

$log2_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120677_01 -var ".data_db_folder()."GS120677_01_annotated.tsv -db NGSD_TEST --log $log2_file");

$log3_file = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120678_01 -var ".data_db_folder()."GS120678_01_annotated.tsv -db NGSD_TEST --log $log3_file");

//germline - import third sample again (with -force flag)
$log4_file = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120678_01 -var ".data_db_folder()."GS120678_01_annotated.tsv -db NGSD_TEST --log $log4_file -force");

//germline - import third sample again (with -force flag) and check that comment/report are still present after re-import
$db = DB::getInstance("NGSD_TEST");
$db->executeStmt("UPDATE detected_variant SET comment='BLA', report='1' WHERE processed_sample_id='3' AND variant_id='2'");

$log5_file = output_folder().$name."_out5.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120678_01 -var ".data_db_folder()."GS120678_01_annotated.tsv -db NGSD_TEST --log $log5_file -force");

$result = $db->executeQuery("SELECT comment, report FROM detected_variant WHERE processed_sample_id='3' AND variant_id='2'");
check($result[0]['comment'], 'BLA');
check($result[0]['report'], 1);
$result = $db->executeQuery("SELECT comment, report FROM detected_variant WHERE processed_sample_id='3' AND variant_id='3'");
check($result[0]['comment'], '');
check($result[0]['report'], 0);
$result = $db->executeQuery("SELECT comment, report FROM detected_variant WHERE processed_sample_id='3' AND variant_id='4'");
check($result[0]['comment'], '');
check($result[0]['report'], 0);
$result = $db->executeQuery("SELECT comment, report FROM detected_variant WHERE processed_sample_id='3' AND variant_id='5'");
check($result[0]['comment'], '');
check($result[0]['report'], 0);

//somatic
$log6_file = output_folder().$name."_out6.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS130796_01-GS130797_01 -var ".data_db_folder()."GS130796_01-GS130797_01.GSvar -db NGSD_TEST -mode somatic --log $log6_file");


end_test();

?>