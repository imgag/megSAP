<?php

require_once("framework.php");

$name = "db_import_variant_classifications";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//first import (create classification)
check_exec("php ".src_folder()."/NGS/db_import_variant_classifications.php -user klaus -in ".data_folder()."/db_import_variant_classifications_in1.tsv -db NGSD_TEST");

//second import (update classification)
check_exec("php ".src_folder()."/NGS/db_import_variant_classifications.php -user klaus -in ".data_folder()."/db_import_variant_classifications_in1.tsv -db NGSD_TEST");

end_test();

?>