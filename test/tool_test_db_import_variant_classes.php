<?php

require_once("framework.php");

$name = "db_import_variant_classes";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test");

//first import (create variants and set classification)
check_exec("php ".src_folder()."/Tools/{$name}.php -user admin -in ".data_folder()."/{$name}_in1.vcf -db NGSD_TEST");

//second import (update classification)
check_exec("php ".src_folder()."/Tools/{$name}.php -user admin -in ".data_folder()."/{$name}_in1.vcf -db NGSD_TEST");


end_test();

?>