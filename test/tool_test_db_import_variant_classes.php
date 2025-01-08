<?php

require_once("framework.php");

$name = "db_import_variant_classes";
start_test($name);

//init NGSD
init_ngsd($name);

//first import (create variants and set classification)
check_exec("php ".src_folder()."/Auxilary/{$name}.php -user admin -in ".data_folder()."/{$name}_in1.vcf -db NGSD_TEST");

//second import (update classification)
check_exec("php ".src_folder()."/Auxilary/{$name}.php -user admin -in ".data_folder()."/{$name}_in1.vcf -db NGSD_TEST");


end_test();

?>