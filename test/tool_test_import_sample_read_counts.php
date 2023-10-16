<?php

require_once("framework.php");

$name = "import_sample_read_counts";

start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
// get id for read count in qc_terms
$db_conn = DB::getInstance("NGSD_TEST");

//test 1 - default NovaSeq 6000
$stats_file = data_folder().$name."_in1.json";
check_exec("php ".src_folder()."/NGS/{$name}.php -stats ${stats_file} -db NGSD_TEST --log ".output_folder()."test1.log");
//check db entries
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=1001 AND qc_terms_id=5"), 119248104);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=1005 AND qc_terms_id=5"), 143768808);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=2002 AND qc_terms_id=5"), 120769616);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=6005 AND qc_terms_id=5"), 129593696);


//test 2 - NextSeq X
$stats_file = data_folder().$name."_in2.csv";
check_exec("php ".src_folder()."/NGS/{$name}.php -stats ${stats_file} -csv_mode -db NGSD_TEST --log ".output_folder()."test2.log");
//check db entries
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=1001 AND qc_terms_id=5"), 283437964);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=1005 AND qc_terms_id=5"), 345253854);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=2002 AND qc_terms_id=5"), 261370342);
check($db_conn->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id=6005 AND qc_terms_id=5"), 317723122);
end_test();

?>
