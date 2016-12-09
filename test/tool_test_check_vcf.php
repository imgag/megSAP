<?php

require_once("framework.php");

$name = "check_vcf";
$script = realpath(src_folder())."/NGS/".$name.".php";

start_test($name);

//no error
$in_file = data_folder()."/".$name."_in1.vcf";
$output = check_exec("php $script -in $in_file", FALSE);
check($output, array());

end_test();

?>