<?php

require_once("framework.php");

$name = "check_tsv";

start_test($name);

$in_file = data_folder()."/".$name."_in1.tsv";
$output = check_exec("php ".src_folder()."/NGS/".$name.".php -in $in_file");
check($output, array());

end_test();

?>