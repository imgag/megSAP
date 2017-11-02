<?php

require_once("framework.php");

$name = "db_import_variants";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//germline - import 3 samples
$log1_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120676_01 -var ".data_db_folder()."GS120676_01_annotated.tsv -db NGSD_TEST --log $log1_file --verbose");

$log2_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120677_01 -var ".data_db_folder()."GS120677_01_annotated.tsv -db NGSD_TEST --log $log2_file --verbose");

$log3_file = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120678_01 -var ".data_db_folder()."GS120678_01_annotated.tsv -db NGSD_TEST --log $log3_file --verbose");

//germline - import third sample again (with -force flag)
$log4_file = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS120678_01 -var ".data_db_folder()."GS120678_01_annotated.tsv -db NGSD_TEST --log $log4_file -force --verbose");

//somatic (tumor-normal)
$log5_file = output_folder().$name."_out5.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS130796_01-GS130797_01 -var ".data_db_folder()."GS130796_01-GS130797_01.GSvar -db NGSD_TEST -mode somatic --log $log5_file --verbose");

//somatic (tumor only)
$log6_file = output_folder().$name."_out6.log";
check_exec("php ".src_folder()."/NGS/db_import_variants.php -id GS130796_01 -var ".data_db_folder()."GS130796_01.GSvar -db NGSD_TEST -mode somatic --log $log6_file --verbose");


end_test();

?>