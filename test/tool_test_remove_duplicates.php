<?php

require_once("framework.php");

$name = "remove_duplicates";
start_test($name);

$in_file = data_folder().$name."_in.bam";
$out_file = output_folder().$name."_out.bam";
$log_file = output_folder().$name."_out.log";

check_exec("php ".src_folder()."/NGS/".$name.".php -in $in_file -out $out_file --log $log_file");

check_file($out_file, data_folder().$name."_out.bam");

end_test();
?>
